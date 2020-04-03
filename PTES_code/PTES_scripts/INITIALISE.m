% Reset Load structure
Load = Load0;

% Reset fluid states and stages
gas = reset_fluid(gas);
if Load.mode==3
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
        CT(ir) = reset_tanks(CT(ir),TC_dis0(ir),10.*p0,MC_dis0(ir),TC_chg0(ir),10.*p0,MC_chg0(ir),T0);
    end
    
    % Reset hot tanks
    for ir = 1 : Nhot
        HT(ir) = reset_tanks(HT(ir),TH_dis0(ir),10.*p0,MH_dis0(ir),TH_chg0(ir),10.*p0,MH_chg0(ir),T0);
    end
end

% Reset atmospheric tanks
AT = reset_tanks(AT,T0,p0,huge,T0,p0,huge,T0);

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
    PR_estim = ((Tmax/T1)^((Gama*eta)/(Gama-1)))^Nc_ch;
    pbot = pmax/PR_estim;
else
    Tmax = T1 * PRch ^((Gama-1)/(Gama*eta))^(1/Nc_ch) ;
end

% Cold temperatures
T3 = TC_dis0(1) ;
if setTmax
    TC_chg0(1) = T3 / (PR_estim ^((Gama-1)/(Gama*eta))^(1/Ne_ch)) ;
else
    TC_chg0(1) = T3 / (PR_ch ^((Gama-1)/(Gama*eta))^(1/Ne_ch)) ;
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
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', 4, eta, Load.num) ; % Charging compressors
        DEXP(1:Nc_ch) = compexp_class('exp', 'poly', 13, eta, Load.num) ; % Discharging expanders
        
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', 13, eta, Load.num) ; % Charging expanders
        DCMP(1:Ne_ch) = compexp_class('comp', 'poly', 4, eta, Load.num) ; % Discharging compressors
    case 3 % JB (charge) + Rankine (discharge)
        % Charging components
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', 4, eta, Load.num) ; % Charging compressors
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', 13, eta, Load.num) ; % Charging expanders
        
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

% Fans --> NOT SURE WHAT cost_mode should be selected in this case
CFAN(1:10) = compexp_class('comp', 'isen', 40, 0.5, Load.num) ;
DFAN(1:10) = compexp_class('comp', 'isen', 40, 0.5, Load.num) ;

% Put design case load cycles in load for the first iteration.
if Loffdesign
    Load = Design_Load ; % Later reset Load using Load0
end

design_mode = 1 ; % Logical - currently in design mode