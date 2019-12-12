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
% Reset tanks
%HT = reset_tanks(HT,TH_dis0,p0,MH_dis0,TH_chg0,p0,MH_chg0,T0);
%CT = reset_tanks(CT,TC_dis0,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0);

% NEW CODE JDM 8 Nov. Reset tanks when there are multiple hot/cold tanks in
% series. Pau to double check this works.
% Reset cold tanks
for ir = 1 : Ncld
    CT(ir) = reset_tanks(CT(ir),TC_dis0(ir),p0,MC_dis0(ir),TC_chg0(ir),p0,MC_chg0(ir),T0);
end

% Reset hot tanks
for ir = 1 : Nhot
    HT(ir) = reset_tanks(HT(ir),TH_dis0(ir),p0,MH_dis0(ir),TH_chg0(ir),p0,MH_chg0(ir),T0);
end

% Set bottom pressure line based on maximum pressure and pressure ratio
pbot = pmax/PRch;
T1   = TH_dis0(1); % compressor inlet temperature estimate
if setTmax
    % Obtain PRch from maximum temperature and estimated temperature ratio
    Gama = CP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CPMASS',gas.handle)/CP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CVMASS',gas.handle);
    PR_estim = ((Tmax/T1)^((Gama*eta)/(Gama-1)))^Nc_ch;
    pbot = pmax/PR_estim;
end

% Set number of compressions and expansions symmetrically between charge
% and discharge
Nc_dis = Ne_ch; % compressions during discharge
Ne_dis = Nc_ch; % expansions during discharge

% Construct compressor and expander classes
switch Load.mode
    case {0,1,2} % Ideal gas Joule-Bratyon PTES
        for i = 1 : Nc_ch
            CCMP(i) = compexp_class('comp', 'poly', 2, eta, Load.num) ; % Charging compressors
            DEXP(i) = compexp_class('exp', 'poly', 11, eta, Load.num) ; % Discharging expanders
        end
        
        for i = 1 : Ne_ch
            CEXP(i) = compexp_class('exp', 'poly', 2, eta, Load.num) ; % Charging expanders
            DCMP(i) = compexp_class('comp', 'poly', 11, eta, Load.num) ; % Discharging compressors
        end
    case 3 % JB (charge) + Rankine (discharge)
        % Charging components
        for i = 1 : Nc_ch
            CCMP(i) = compexp_class('comp', 'poly', 2, eta, Load.num) ; % Charging compressors
        end
        
        for i = 1 : Ne_ch
            CEXP(i) = compexp_class('exp', 'poly', 11, eta, Load.num) ; % Charging expanders
        end
        
        % Discharging components
        for i = 1 : 3
            DCMP(i) = compexp_class('comp', 'isen', 22, eta, Load.num) ; % Discharging compressors
            DEXP(i) = compexp_class('exp', 'isen', 11, eta, Load.num) ; % Discharging expanders
        end
        
    case 4 % sCO2-PTES type cycles
        for i = 1 : Nc_ch
            CCMP(i) = compexp_class('comp', 'poly', 7, eta, Load.num) ; % Charging compressors
            DEXP(i) = compexp_class('exp', 'poly', 17, eta, Load.num) ; % Discharging expanders
        end
        
        for i = 1 : Ne_ch
            CEXP(i) = compexp_class('exp', 'poly', 7, eta, Load.num) ; % Charging expanders
            DCMP(i) = compexp_class('comp', 'poly', 17, eta, Load.num) ; % Discharging compressors
        end
        
        %Recompressor
        if Lrcmp
            RCMP = compexp_class('comp', 'poly', 7, eta, Load.num) ; % Re-compressors
        end
end
