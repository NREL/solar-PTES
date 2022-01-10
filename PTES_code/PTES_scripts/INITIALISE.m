% Reset Load structure
Load = Load0;

switch Load.mode
    case {0,1,2,3,4,5,6,7,20,22,24}
        
        % Reset fluid states and stages
        gas = reset_fluid(gas);
        if any(Load.mode==[3,7])
            steam = reset_fluid(steam);
        end
        if Load.mode == 20
            gas2 = reset_fluid(gas2);
        end
        
        % These storage fluid streams only exist in certain cases
        if PBmode == 0 || PBmode == 2
            for ir = 1:length(fluidH)
                fluidH(ir) = reset_fluid(fluidH(ir)); %#ok<*SAGROW>
            end
            for ir = 1:length(fluidM)
                fluidM(ir) = reset_fluid(fluidM(ir));
            end
            for ir = 1:length(fluidC)
                fluidC(ir) = reset_fluid(fluidC(ir));
            end
            
            % Reset cold tanks
            for ir = 1 : Ncld
                CT(ir) = reset_tanks(CT(ir),TC_dis0(ir),p0,MC_dis0(ir),TC_chg0(ir),p0,MC_chg0(ir),T0);
            end
            
            % Reset medium tanks
            for ir = 1 : Nmid
                MT(ir) = reset_tanks(MT(ir),TM_dis0(ir),p0,MM_dis0(ir),TM_chg0(ir),p0,MM_chg0(ir),T0);
            end
            
            % Reset hot tanks
            for ir = 1 : Nhot
                HT(ir) = reset_tanks(HT(ir),TH_dis0(ir),p0,MH_dis0(ir),TH_chg0(ir),p0,MH_chg0(ir),T0);
            end
            
            if any(Load.mode == [20,22,24])
                % Reset liquid air tanks
                LAT(ir) = reset_single_tank(LAT(ir),TLA_chg0(ir),p0,MLA_chg0(ir),T0);
                
                % Reset atmospheric tank
                AT = reset_single_tank(AT,TAT_chg0,p0,MAT_chg0,T0);
                
            else 
                % Reset atmospheric tanks
                AT = reset_tanks(AT,T0,p0,0,T0,p0,0,T0);
            end
        end
        
    otherwise
        error('not implemented')
end

switch Load.mode
    case {0,1,2,3,4,5,6,7}
        
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
                Gama = RPN('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CPMASS',gas)/RPN('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CVMASS',gas);
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
end

% Construct compressor and expander classes
switch Load.mode
    case {0,1,2} % Ideal gas Joule-Brayton PTES
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', CCMPmode(1), eta, Load.num) ; % Charging compressors
        DEXP(1:Nc_ch) = compexp_class('exp', 'poly', DEXPmode(1), eta, Load.num) ; % Discharging expanders
        
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', CEXPmode(1), eta, Load.num) ; % Charging expanders
        DCMP(1:Ne_ch) = compexp_class('comp', 'poly', DCMPmode(1), eta, Load.num) ; % Discharging compressors
    case {3,7} % JB (charge) + Rankine (discharge)
        % Charging components
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', CCMPmode(1), eta, Load.num) ; % Charging compressors
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', CEXPmode(1), eta, Load.num) ; % Charging expanders
        
        % Discharging components
        if ~Load.options.superRank
            DCMP(1:3) = compexp_class('comp', 'isen', 0, eta, Load.num) ; % Discharging compressors
            DEXP(1:3) = compexp_class('exp', 'isen', 0, eta, Load.num) ; % Discharging expanders
        else
            DCMP(1:4) = compexp_class('comp', 'isen', 0, eta, Load.num) ; % Discharging compressors
            DEXP(1:4) = compexp_class('exp', 'isen', 0, eta, Load.num) ; % Discharging expanders
        end
        
    case {4,5} % sCO2-PTES type cycles
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', CCMPmode(1), eta, Load.num) ; % Charging compressors
        DEXP(1:Nc_ch) = compexp_class('exp', 'poly', DEXPmode(1), eta, Load.num) ; % Discharging expanders
        
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', CEXPmode(1), eta, Load.num) ; % Charging expanders
        DCMP(1:Ne_ch) = compexp_class('comp', 'poly', DCMPmode(1), eta, Load.num) ; % Discharging compressors
        
        %Recompressor
        if Lrcmp
            RCMP = compexp_class('comp', 'poly', RCMPmode(1), eta, Load.num) ; % Re-compressors
        end
        
    case {6} % sCO2-PTES type cycles
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', CCMPmode(1), eta, Load.num) ; % Charging compressors
        DEXP(1:Nc_ch) = compexp_class('exp', 'poly', DEXPmode(1), eta, Load.num) ; % Discharging expanders
        
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', CEXPmode(1), eta, Load.num) ; % Charging expanders
        DCMP(1:Ne_ch) = compexp_class('comp', 'poly', DCMPmode(1), eta, Load.num) ; % Discharging compressors
        
        %Recompressor
        if Lrcmp
            RCMP = compexp_class('comp', 'poly', CCMPmode(1), eta, Load.num) ; % Re-compressors
        end
        
    case 20 % CCES
        CCMP(1:2) = compexp_class('comp', 'poly', CCMPmode(1), eta, Load.num) ; % Charging compressors (LAES)
        CCMP(3:4) = compexp_class('comp', 'poly', CCMPmode(1), eta, Load.num) ; % Charging compressors (PTES)
        CEXP(1)   = compexp_class('exp',  'isen', CEXPmode(1), eta_cryo_exp, Load.num) ; % Cryo expander (LAES side)
        CEXP(2:4) = compexp_class('exp',  'poly', CEXPmode(1), eta, Load.num) ; % Charging expanders (PTES side)
        
        DCMP(1)   = compexp_class('comp', 'isen', DCMPmode(1), eta_cryo_pump, Load.num) ; % Cryo pump (LAES side)
        DCMP(2:5) = compexp_class('comp', 'poly', DCMPmode(1), eta, Load.num) ; % Disch. compressors (PTES side)
        DEXP(1:2) = compexp_class('exp',  'poly', DEXPmode(1), eta, Load.num) ; % Discharging expanders (LAES)
        DEXP(3)   = compexp_class('exp',  'poly', DEXPmode(1), eta, Load.num) ; % Discharging expanders (PTES)
        
    case {22,24} % Integrated PLAES
        % Charging turbomachines:
        % 1 shared compressor, 1 LAES compressor,
        % 1 cryogenic (liquid) expander, and 2 PTES turbines
        CCMP(1:2) = compexp_class('comp', 'poly', CCMPmode(1), eta, Load.num) ;
        CEXP(1)   = compexp_class('exp',  'isen', CEXPmode(1), eta_cryo_exp, Load.num) ;
        CEXP(2:3) = compexp_class('exp',  'poly', CEXPmode(1), eta, Load.num) ;
        
        % Discharging turbomachines:
        % 1 cryogenic pump; 2 PTES compressors and 1 auxiliary compressor;
        % and % 1 LAES turbine and 1 shared turbine.
        DCMP(1)   = compexp_class('comp', 'isen', DCMPmode(1), eta_cryo_pump, Load.num) ;
        DCMP(2:4) = compexp_class('comp', 'poly', DCMPmode(1), eta, Load.num) ;
        DEXP(1:2) = compexp_class('exp',  'poly', DEXPmode(1), eta, Load.num) ;
end

% Construct heat exchangers
switch Load.mode
    case 0
        % Call HX classes for ideal-gas PTES cycle
        % First, set heat exchanger indeces
        ihx_hot  = 1:Nc_ch;
        ihx_reg  = ihx_hot(end)+1;        
        ihx_cld  = ihx_reg(end)+(1:Ne_ch);
        ihx_rejc = ihx_cld(end)+1;
        ihx_rejd = ihx_rejc(end)+(1:Ne_ch);
        switch PBmode
            case {0,2}
                HX(ihx_hot)  = hx_class('hot',  'hex',   hotHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Hot heat exchanger
                HX(ihx_reg)  = hx_class('regen','regen', rcpHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Recuperator
                HX(ihx_cld)  = hx_class('cold', 'hex',   cldHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Cold heat exchanger
                HX(ihx_rejc) = hx_class('rej',  'hex',   rejHXmode(1), HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shapeX) ; % Heat rejection unit (charge)
                HX(ihx_rejd) = hx_class('rej',  'hex',   rejHXmode(1), HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shapeX) ; % Heat rejection unit (discharge)
                
            case 1
                HX(1) = hx_class('rej',  'hex',   rejHXmode(1), 100, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shapeX) ; % Heat rejection unit
                HX(2) = hx_class('rej',  'hex',   rejHXmode(1), 100, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shapeX) ; % Heat rejection unit
        end
        
    case 1
        % Call HX classes for ideal-gas PTES heat pump
        % First, set heat exchanger indeces
        ihx_hot  = 1:Nc_ch;
        ihx_reg  = ihx_hot(end)+1;
        ihx_cld  = ihx_reg(end)+(1:Ne_ch);
        ihx_rejc = ihx_cld(end)+1;
        HX(ihx_hot)  = hx_class('hot',  'hex',   hotHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Hot heat exchanger
        HX(ihx_reg)  = hx_class('regen','regen', rcpHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Recuperator
        HX(ihx_cld)  = hx_class('cold', 'hex',   cldHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Cold heat exchanger
        HX(ihx_rejc) = hx_class('rej',  'hex',   rejHXmode(1), HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shapeX) ; % Heat rejection unit (charge)
        
    case 2
        % Call HX classes for ideal-gas PTES heat engine
        % First, set heat exchanger indeces
        ihx_hot  = 1:Nc_ch;
        ihx_reg  = ihx_hot(end)+1;
        ihx_rejd = ihx_reg(end)+Ne_ch;
        HX(ihx_hot)  = hx_class('hot',  'hex',   hotHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Hot heat exchanger
        HX(ihx_reg)  = hx_class('regen','regen', rcpHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Recuperator
        HX(ihx_rejd) = hx_class('rej',  'hex',   rejHXmode(1), HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shapeX) ; % Heat rejection unit (discharge)
        
    case {3,7}
        % Call HX classes for ideal-gas PTES heat pump with Rankine cycle discharge
        ihx_hot = 1:Nc_ch;
        ihx_reg = ihx_hot(end)+1;
        %ihx_rejc= ihx_reg(end)+1;
        %ihx_hin = ihx_rejc(end)+(1:Ne_ch);
        ihx_hin = ihx_reg(end)+(1:Ne_ch);
        ihx_cld = ihx_hin(end)+(1:Ne_ch);
        ihx_htf = ihx_cld(end)+(1:Ne_ch);
        ihx_JB  = ihx_htf(end);
        HX(ihx_hot) = hx_class('hot',  'hex',   hotHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Hot heat exchanger
        HX(ihx_reg) = hx_class('regen','regen', rcpHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Recuperator
        %HX(ihx_rejc)= hx_class('rej',  'hex',   rejHXmode(1), HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shapeX) ; % Heat rejection unit (charge)
        HX(ihx_hin) = hx_class('hin',  'hex',   rejHXmode(1), HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shape)  ; % Heat intake unit      (no cold tanks scenario)
        HX(ihx_cld) = hx_class('cold', 'hex',   cldHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Cold heat exchanger   (cold tanks)
        HX(ihx_htf) = hx_class('htf',  'hex',   cldHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Intermediate HTF loop (cold tanks)
        
        
        
        % Normal Rankine
        if ~Load.options.superRank
            HX(ihx_JB+1) = hx_class('hot',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Reheat
            HX(ihx_JB+2) = hx_class('cold', 'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Condenser
            HX(ihx_JB+3) = hx_class('rej',  'regen', 0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shapeX) ; % Air-cooled condenser
            HX(ihx_JB+4) = hx_class('hot',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Boiler
        else
            % Supercritical Rankine
            HX(ihx_JB+1) = hx_class('hot',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % First Reheat
            HX(ihx_JB+2) = hx_class('hot',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Second Reheat
            HX(ihx_JB+3) = hx_class('cold', 'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Condenser
            HX(ihx_JB+4) = hx_class('rej',  'regen', 0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shapeX) ; % Air-cooled condenser
            HX(ihx_JB+5) = hx_class('hot',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape)  ; % Boiler
        end
        
    case 4
        
        % Should implement a scheme similar to JB_PTES with ihx indices
        iHX = 1 ; % Heat exchanger counter
        for ii = 1 : Nhot
            HX(iHX) = hx_class('hot', 'hex', hotHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
            iHX = iHX + 1 ;
        end
        for ii = 1 : Ncld
            HX(iHX) = hx_class('cold', 'hex', cldHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
            iHX = iHX + 1 ;
        end
        if (Nhot < 2) && (Ncld < 2) && (Nrcp > 0)
            for ii = 1 : Nrcp
                HX(iHX) = hx_class('regen', 'regen', rcpHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
                iHX = iHX + 1 ;
            end
        end
        for ii = 1 : 3 % How many heat rejection units?
            HX(iHX) = hx_class('rej',  'hex',   rejHXmode(1), HX_NX, Load.num, Load.num, 'eff', eff, 0.001, HX_D1, HX_shape) ; % Heat rejection unit 
            iHX     = iHX + 1 ;
        end
    
    case 5
        % Heat exchangers set up to match Ty's work
        %HX(1) = hx_class('hot',   'hex',   25, HX_NX, Load.num, Load.num, 'eff', 0.879, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        %HX(2) = hx_class('cold',  'hex',   25, HX_NX, Load.num, Load.num, 'eff',   eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        %HX(3) = hx_class('regen', 'regen', 25, HX_NX, Load.num, Load.num, 'eff', 0.968, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        %HX(4) = hx_class('regen', 'regen', 25, HX_NX, Load.num, Load.num, 'eff', 0.937, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        
        HX(1) = hx_class('hot',   'hex',   25, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(2) = hx_class('cold',  'hex',   25, HX_NX, Load.num, Load.num, HX_model,   eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(3) = hx_class('regen', 'regen', 25, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(4) = hx_class('regen', 'regen', 25, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        
        
    case 6
        HX(1) = hx_class('hot',   'hex',    0, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(2) = hx_class('hot',   'hex',   hotHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Hot heat exchanger
        HX(3) = hx_class('cold',  'hex',   cldHXmode(1), HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Cold heat exchanger
        HX(4) = hx_class('regen', 'regen',  0, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Recuperator exchanger
        HX(5) = hx_class('regen', 'regen',  0, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Recuperator heat exchanger    
        
    case 20
        % HEXs for the LAES side
        HX(1) = hx_class('hot', 'hex', 0, HX_NX, Load.num, Load.num, HX_model, 1-(1-eff)/2,  ploss,  HX_D1, HX_shape) ; % Hot heat exchanger 1
        HX(2) = hx_class('hot', 'hex', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Hot heat exchanger 2
        HX(3) = hx_class('rej', 'hex', 0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shape) ; % Heat rejection unit
        
        % Couplers
        HX(4) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Coupler A (air-neon)
        HX(5) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Coupler B (air-neon)
        HX(6) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff, ploss, HX_D1, HX_shape) ; % Coupler C (air-neon)
        
        % HEXs for the PTES side
        HX(7) = hx_class('hot',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Hot heat exchanger
        HX(8) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Regen A
        HX(9) = hx_class('rej',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shape) ; % Heat rejection unit
        HX(10)= hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Regen B
        HX(11)= hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Regen C
        
    case {22}
        % Shared hot HEX and Heat Rejection Unit
        % EXTRA EFFECTIVE (AND DOUBLE PRESSURE LOSS) AS IT SIMULATES TWO HEXS IN SERIES!
        HX(1) = hx_class('hot', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, 1-(1-eff)/2,  2*ploss,  HX_D1, HX_shape) ; % Hot heat exchanger A
        %HX(1) = hx_class('hot', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Hot heat exchanger A
        HX(2) = hx_class('rej',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shape) ; % Heat rejection unit
        
        % LAES hot HEX and Heat Rejection Unit
        % EXTRA EFFECTIVE (AND DOUBLE PRESSURE LOSS) AS IT SIMULATES TWO HEXS IN SERIES!
        HX(3) = hx_class('hot', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, 1-(1-eff)/2,  2*ploss,  HX_D1, HX_shape) ; % Hot heat exchanger B
        %HX(3) = hx_class('hot', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Hot heat exchanger B
        HX(4) = hx_class('rej', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shape) ; % Heat rejection unit
        
        % Couplers
        HX(5) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Coupler A
        HX(6) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Coupler B
        
        % Regen for the PTES side and discharge Heat Rejection Unit
        HX(7) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Regen
        HX(8) = hx_class('rej', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shape) ; % Heat rejection unit
        
        % Discharge coupler
        HX(9) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, 1-(1-eff)/2,  2*ploss,  HX_D1, HX_shape) ; % Coupler B
        
    case {24}
        % Shared hot HEX and Heat Rejection Unit
        % EXTRA EFFECTIVE (AND DOUBLE PRESSURE LOSS) AS IT SIMULATES TWO HEXS IN SERIES!
        HX(1) = hx_class('hot', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, 1-(1-eff)/2,  2*ploss,  HX_D1, HX_shape) ; % Hot heat exchanger A
        HX(2) = hx_class('rej',  'hex',   0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shape) ; % Heat rejection unit
        
        % LAES hot HEX and Heat Rejection Unit
        HX(3) = hx_class('hot', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Hot heat exchanger B
        HX(4) = hx_class('rej', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shape) ; % Heat rejection unit
        
        % Cold HEX
        HX(5) = hx_class('cold','hex', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Coupler A
        
        % Coupler B
        HX(6) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Coupler B
        
        % Regen for the PTES side and discharge Heat Rejection Unit
        HX(7) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, eff,  ploss,  HX_D1, HX_shape) ; % Regen
        HX(8) = hx_class('rej', 'hex',    0, HX_NX, Load.num, Load.num, HX_model, effX, plossX, HX_D1, HX_shape) ; % Heat rejection unit
        
        % Discharge coupler
        HX(9) = hx_class('regen','regen', 0, HX_NX, Load.num, Load.num, HX_model, 1-(1-eff)/2,  2*ploss,  HX_D1, HX_shape) ; % Coupler B
end

% Motor-generator
if Load.mode == 3 || Load.mode == 6
    % Probably have a separate generator and motor for solar-PTES
    GEN(1)     = gen_class('mot',GENmode(1)) ;
    GEN(2)     = gen_class('gen',0) ; % Generator costs nothing as part of Rankine cycle already
else
    GEN        = gen_class('mot-gen',GENmode(1)) ;
end
    

% Fans 
CFAN(1:10) = compexp_class('comp', 'poly', FANmode(1), 0.75, Load.num) ;
DFAN(1:10) = compexp_class('comp', 'poly', FANmode(1), 0.75, Load.num) ;

% Fluid pumps 
CPMP(1:10) = compexp_class('pump', 'isen', PMPmode(1), 0.8, Load.num) ;
DPMP(1:10) = compexp_class('pump', 'isen', PMPmode(1), 0.8, Load.num) ;

% Mixers may be required (e.g. in Rankine cycle)
MIX(1:6)   = misc_class('mix',Load.num) ;

% Put design case load cycles in load for the first iteration.
if Loffdesign
    Load = Design_Load ; % Later reset Load using Load0
end

design_mode = 1 ; % Logical - currently in design mode
