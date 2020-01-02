% Reset Load structure
Load = Load0;

% Reset fluid states and stages
gas = reset_fluid(gas);
if Load.mode==3
    steam = reset_fluid(steam);
end
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

% Reset atmospheric tanks
AT = reset_tanks(AT,T0,p0,huge,T0,p0,huge,T0);

% Set bottom pressure line based on maximum pressure and pressure ratio
pbot = pmax/PRch;
T1   = TH_dis0(1); % compressor inlet temperature estimate
if setTmax
    % Obtain PRch from maximum temperature and estimated temperature ratio
    if strcmp(gas.read,'IDL')
        Gama = gas.IDL.gam ;
    else
        Gama = CP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CPMASS',gas.handle)/CP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CVMASS',gas.handle);
    end
    PR_estim = ((Tmax/T1)^((Gama*eta)/(Gama-1)))^Nc_ch;
    pbot = pmax/PR_estim;
end

% Set number of compressions and expansions symmetrically between charge
% and discharge
Nc_dis = Ne_ch; % compressions during discharge
Ne_dis = Nc_ch; % expansions during discharge

% Construct compressor and expander classes
switch Load.mode
    case {0,1,2} % Ideal gas Joule-Brayton PTES
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', 2, eta, Load.num) ; % Charging compressors
        DEXP(1:Nc_ch) = compexp_class('exp', 'poly', 11, eta, Load.num) ; % Discharging expanders
        
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', 2, eta, Load.num) ; % Charging expanders
        DCMP(1:Ne_ch) = compexp_class('comp', 'poly', 11, eta, Load.num) ; % Discharging compressors
    case 3 % JB (charge) + Rankine (discharge)
        % Charging components
        CCMP(1:Nc_ch) = compexp_class('comp', 'poly', 2, eta, Load.num) ; % Charging compressors
        CEXP(1:Ne_ch) = compexp_class('exp', 'poly', 11, eta, Load.num) ; % Charging expanders
        
        % Discharging components
        DCMP(1:3) = compexp_class('comp', 'isen', 22, eta, Load.num) ; % Discharging compressors
        DEXP(1:3) = compexp_class('exp', 'isen', 11, eta, Load.num) ; % Discharging expanders
        
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
