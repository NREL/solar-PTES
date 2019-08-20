% Reset fluid states and stages
gas = reset_fluid(gas);
for ir = 1:length(fluidH)
    fluidH(ir) = reset_fluid(fluidH(ir)); %#ok<*SAGROW>
end
for ir = 1:length(fluidC)
    fluidC(ir) = reset_fluid(fluidC(ir));
end
% Reset tanks
HT = reset_tanks(HT,TH_0,p0);
CT = reset_tanks(CT,TC_0,p0);


% Set bottom pressure line based on maximum pressure and pressure ratio
pbot = pmax/PRch;
T1   = TH_0; % compressor inlet temperature estimate
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