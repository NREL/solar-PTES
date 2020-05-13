% Reset Load structure
Load = Load0;

% Reset fluid states and stages
gas = reset_fluid(gas);
if any(Load.mode==[3,7])
    steam = reset_fluid(steam);
end

% These storage fluid streams only exist in certain cases
if PBmode == 0 || PBmode == 2
    for ir = 1:length(fluidH)
        fluidH(ir) = reset_fluid(fluidH(ir)); %#ok<*SAGROW>
    end
    for ir = 1:length(fluidC)
        fluidC(ir) = reset_fluid(fluidC(ir));
    end
    
    % Reset cold tanks
    for ir = 1 : Ncld
        CT(ir) = reset_tanks(CT(ir),TC_dis0(ir),p0,MC_dis0(ir),TC_chg0(ir),p0,MC_chg0(ir),T0);
    end
    
    % Reset hot tanks
    for ir = 1 : Nhot
        HT(ir) = reset_tanks(HT(ir),TH_dis0(ir),p0,MH_dis0(ir),TH_chg0(ir),p0,MH_chg0(ir),T0);
    end
end

% Reset atmospheric tanks
AT = reset_tanks(AT,T0,p0,0,T0,p0,0,T0);

% Estimate cycle pressures and temperatures
% Set bottom pressure line based on maximum pressure and pressure ratio
pbot = pmax/PRch;
T1   = TH_dis0(1); % compressor inlet temperature estimate
% Estimate gamma
if strcmp(gas.read,'IDL')
    Gama = gas.IDL.gam ;
else
    Gama = CP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CPMASS',gas.handle)/CP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CVMASS',gas.handle);
end

% Hot temperatures
if setTmax
    % Obtain PRch from maximum temperature and estimated temperature ratio
    if strcmp(gas.read,'IDL')
        Gama = gas.IDL.gam ;
    else
        Gama = RP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CPMASS',gas)/RP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CVMASS',gas);
    end
    PRestim = ((Tmax/T1)^((Gama*eta)/(Gama-1)))^Nc_ch;
    pbot = pmax/PRestim;
else
    Tmax = T1 * PRch ^((Gama-1)/(Gama*eta))^(1/Nc_ch) ;
end

% Cold temperatures
T3 = TC_dis0(1) ;
if setTmax
    TC_chg0(1) = T3 / (PRestim ^((Gama-1)*eta/(Gama))^(1/Ne_ch)) ;
else
    TC_chg0(1) = T3 / (PRch ^((Gama-1)*eta/(Gama))^(1/Ne_ch)) ;
end

switch PBmode
    case 1
        % Set up packed beds
        for ii = 1 : Nhot 
            pbH(ii).TC = TH_chg0 ;
            pbH(ii).TD = TH_dis0 ;
            pbH(ii) = PB_INITIALISE( pbH(ii), gas, pmax, Load ) ;
        end
        
        for ii = 1 : Ncld
            pbC(ii).TC = TC_chg0 ;
            pbC(ii).TD = TC_dis0 ;
            pbC(ii) = PB_INITIALISE( pbC(ii), gas, pbot, Load ) ;
        end
        
        [pbH, pbC] = PB_TIMINGS(pbH, pbC) ; % This routine selects the smallest timesteps to run at
    case 2
        error ('Not implemented')
end

% Set number of compressions and expansions symmetrically between charge
% and discharge
Nc_dis = Ne_ch; % compressions during discharge
Ne_dis = Nc_ch; % expansions during discharge

% Construct compressor and expander classes
switch Load.mode
    case {0,1,2} % Ideal gas Joule-Brayton PTES
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', CCMPmode, eta, Load.num) ; % Charging compressors
        DEXP(1:Nc_ch) = compexp_class('exp', 'poly', DEXPmode, eta, Load.num) ; % Discharging expanders
        
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', CEXPmode, eta, Load.num) ; % Charging expanders
        DCMP(1:Ne_ch) = compexp_class('comp', 'poly', DCMPmode, eta, Load.num) ; % Discharging compressors
    case {3,7} % JB (charge) + Rankine (discharge)
        % Charging components
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', CCMPmode, eta, Load.num) ; % Charging compressors
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', CEXPmode, eta, Load.num) ; % Charging expanders
        
        % Discharging components
        DCMP(1:3) = compexp_class('comp', 'isen', 0, eta, Load.num) ; % Discharging compressors
        DEXP(1:3) = compexp_class('exp', 'isen', 0, eta, Load.num) ; % Discharging expanders
        
    case {4,5} % sCO2-PTES type cycles
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', 7, eta, Load.num) ; % Charging compressors
        DEXP(1:Nc_ch) = compexp_class('exp', 'poly', 17, eta, Load.num) ; % Discharging expanders
        
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', 7, eta, Load.num) ; % Charging expanders
        DCMP(1:Ne_ch) = compexp_class('comp', 'poly', 17, eta, Load.num) ; % Discharging compressors
        
        %Recompressor
        if Lrcmp
            RCMP = compexp_class('comp', 'poly', 7, eta, Load.num) ; % Re-compressors
        end
        
    case {6} % sCO2-PTES type cycles
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', 0, eta, Load.num) ; % Charging compressors
        DEXP(1:Nc_ch) = compexp_class('exp', 'poly', 0, eta, Load.num) ; % Discharging expanders
        
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', 17, eta, Load.num) ; % Charging expanders
        DCMP(1:Ne_ch) = compexp_class('comp', 'poly', 7, eta, Load.num) ; % Discharging compressors
        
        %Recompressor
        if Lrcmp
            RCMP = compexp_class('comp', 'poly', 0, eta, Load.num) ; % Re-compressors
        end
        
end

% Construct heat exchangers
ihx_hot = 1:Nc_ch;
ihx_reg = ihx_hot(end)+1;
ihx_rej = ihx_reg(end)+1;
ihx_cld = ihx_rej(end)+(1:Ne_ch);
switch Load.mode
    case {0,1,2}
        switch PBmode
            case {0,2}
                % Call HX classes for ideal-gas PTES cycle
                HX(ihx_hot)  = hx_class('hot',  'hex',   hotHXmode, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
                HX(ihx_reg)  = hx_class('regen','regen', rcpHXmode, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Recuperator
                HX(ihx_rej)  = hx_class('rej',  'hex',   rejHXmode, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Heat rejection unit
                HX(ihx_cld)  = hx_class('cold', 'hex',   cldHXmode, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Cold heat exchanger
            case 1
                HX(1) = hx_class('rej',  'hex',   rejHXmode, 100, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Heat rejection unit
                HX(2) = hx_class('rej',  'hex',   rejHXmode, 100, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Heat rejection unit
        end
        
        
    case {3,7}
        % Call HX classes for ideal-gas PTES heat pump with Rankine cycle discharge
        ihx_JB  = ihx_cld(end);
        HX(ihx_hot)  = hx_class('hot',  'hex',   hotHXmode, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(ihx_reg)  = hx_class('regen','regen', rcpHXmode, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Recuperator
        HX(ihx_rej)  = hx_class('rej',  'hex',   rejHXmode, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Heat rejection unit
        HX(ihx_cld)  = hx_class('cold', 'hex',   cldHXmode, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Cold heat exchanger
        
        HX(ihx_JB+1) = hx_class('hot',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Reheat
        HX(ihx_JB+2) = hx_class('cold', 'hex',   0, HX_NX, Load.num, Load.num, 'eff', eff, 0.1/100, HX_D1, HX_shape) ; % Condenser
        HX(ihx_JB+3) = hx_class('rej',  'regen', 0, HX_NX, Load.num, Load.num, 'eff', eff, 0.1/100, HX_D1, HX_shape) ; % Air-cooled condenser
        HX(ihx_JB+4) = hx_class('hot',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Boiler
                
    case 4
        iHX = 1 ; % Heat exchanger counter
    for ii = 1 : Nhot
        HX(iHX) = hx_class('hot', 'hex', 25, HX_NX, Load.num, Load.num, 'eff', eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        iHX = iHX + 1 ;
    end
    for ii = 1 : Ncld
        HX(iHX) = hx_class('cold', 'hex', 25, HX_NX, Load.num, Load.num, 'eff', eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        iHX = iHX + 1 ;
    end
    if (Nhot < 2) && (Ncld < 2) && (Nrcp > 0)
        for ii = 1 : Nrcp
            HX(iHX) = hx_class('regen', 'regen', 25, HX_NX, Load.num, Load.num, 'eff', eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
            iHX = iHX + 1 ;
        end
    end
    
    case 5
        % Heat exchangers set up to match Ty's work
        HX(1) = hx_class('hot',   'hex',   25, HX_NX, Load.num, Load.num, 'eff', 0.879, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(2) = hx_class('cold',  'hex',   25, HX_NX, Load.num, Load.num, 'eff',   eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(3) = hx_class('regen', 'regen', 25, HX_NX, Load.num, Load.num, 'eff', 0.968, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(4) = hx_class('regen', 'regen', 25, HX_NX, Load.num, Load.num, 'eff', 0.937, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        
    case 6
        HX(1) = hx_class('hot',   'hex',    0, HX_NX, Load.num, Load.num, 'eff', eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(2) = hx_class('hot',   'hex',   25, HX_NX, Load.num, Load.num, 'eff', eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(3) = hx_class('cold',  'hex',   25, HX_NX, Load.num, Load.num, 'eff', eff, ploss, HX_D1, HX_shape) ; % Cold heat exchanger
        HX(4) = hx_class('regen', 'regen',  0, HX_NX, Load.num, Load.num, 'eff', eff, ploss, HX_D1, HX_shape) ; % Recuperator exchanger
        HX(5) = hx_class('regen', 'regen',  0, HX_NX, Load.num, Load.num, 'eff', eff, ploss, HX_D1, HX_shape) ; % Recuperator heat exchanger    
        
end


% Fans --> NOT SURE WHAT cost_mode should be selected in this case
CFAN(1:10) = compexp_class('comp', 'isen', FANmode, 0.5, Load.num) ;
DFAN(1:10) = compexp_class('comp', 'isen', FANmode, 0.5, Load.num) ;

% Fluid pumps --> NOT SURE WHAT cost_mode should be selected in this case
CPMP(1:10) = compexp_class('pump', 'isen', PMPmode, 0.8, Load.num) ;
DPMP(1:10) = compexp_class('pump', 'isen', PMPmode, 0.8, Load.num) ;


% Mixers may be required (e.g. in Rankine cycle)
MIX(1:2)    = misc_class('mix',Load.num) ;

% Put design case load cycles in load for the first iteration.
if Loffdesign
    Load = Design_Load ; % Later reset Load using Load0
end

design_mode = 1 ; % Logical - currently in design mode
