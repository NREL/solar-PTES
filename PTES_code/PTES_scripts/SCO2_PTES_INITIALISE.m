% Reset fluid states and stages
gas = reset_fluid(gas);
for ir = 1:length(fluidH)
    fluidH(ir) = reset_fluid(fluidH(ir)); %#ok<*SAGROW>
end
for ir = 1:length(fluidC)
    fluidC(ir) = reset_fluid(fluidC(ir));
end
% Reset tanks

CT = reset_tanks(CT,TC_dis0,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0);

switch Ncld
    case 1
        CT = reset_tanks(CT,TC_dis0,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0);
    case 2
        CT  = reset_tanks(CT,TC_int,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0);
        CT2 = reset_tanks(CT,TC_dis0,p0,MC_dis0,TC_int,p0,MC_chg0,T0);
    case 3
        error('Not implemented')
end
% Hot tanks
switch Nhot
    case 1
        HT = reset_tanks(HT,TH_dis0,p0,MH_dis0,TH_chg0,p0,MH_chg0,T0);
    case 2
        HT  = reset_tanks(HT,TH_int,p0,MH_dis0,TH_chg0,p0,MH_chg0,T0);
        HT2 = reset_tanks(HT,TH_dis0,p0,MH_dis0,TH_int,p0,MH_chg0,T0);
    case 3
        error('Not implemented')
end


% Set bottom pressure line based on maximum pressure and pressure ratio
pbot = pmax/PRch;
T1   = TH_chg0; % compressor inlet temperature estimate
if setTmax
    % Obtain PRch from maximum temperature and estimated temperature ratio
    Gama = CP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CPMASS',gas.handle)/CP1('PT_INPUTS',pmax/2,0.5*(T1+Tmax),'CVMASS',gas.handle);
    PR_estim = ((Tmax/T1)^((Gama*eta)/(Gama-1)))^Nc_ch;
    pbot = pmax/PR_estim;
end

% Set number of compressions and expansions simmetrically between charge
% and discharge
Nc_dis = Ne_ch; % compressions during discharge
Ne_dis = Nc_ch; % expansions during discharge