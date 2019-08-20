% Reset fluid states and stages
gas1  = reset_fluid(gas1);
gas2  = reset_fluid(gas2);
gas3  = reset_fluid(gas3);
gasHP = reset_fluid(gasHP);
for ir = 1:length(fluidH)
    fluidH(ir) = reset_fluid(fluidH(ir)); %#ok<*SAGROW>
end
for ir = 1:length(fluidC1)
    fluidC1(ir) = reset_fluid(fluidC1(ir));
end
for ir = 1:length(fluidC2)
    fluidC2(ir) = reset_fluid(fluidC2(ir));
end
% Reset tanks
HT  = reset_tanks(HT,TH_0,p0);
CT1 = reset_tanks(CT1,TC1_0,p0);
CT2 = reset_tanks(CT2,TC2_0,p0);


% % Set bottom pressure line based on maximum pressure and pressure ratio
% pbot = pmax/PR;
% T1   = TH_0; % compressor inlet temperature estimate
% if setTmax
%     % Obtain PR from maximum temperature and estimated temperature ratio
%     Gamma = CoolProp.PropsSI('CPMASS','P',pmax/2,'T',0.5*(T1+Tmax),gas.name)/CoolProp.PropsSI('CVMASS','P',pmax/2,'T',0.5*(T1+Tmax),gas.name);
%     PR_estim = ((Tmax/T1)^((Gamma*eta)/(Gamma-1)))^Nc_ch;
%     pbot = pmax/PR_estim;
% end