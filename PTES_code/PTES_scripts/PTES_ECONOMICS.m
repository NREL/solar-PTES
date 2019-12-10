% Script that calculates the cost of each component and subsequently the
% total cost of the system, as well as other economic metrics (possibly
% ...)
cap_cost = 0 ;
Nsens    = 10000 ; % How many points to take from distribution for sensitivity analysis
cap_sens = zeros(Nsens,1) ;

% Compressors and expanders
for ii = 1 : Nc_ch %% THIS DOESN'T WORK FOR RANKINE DISCHARGE
    CCMP(ii) = compexp_econ(CCMP(ii))  ;
    DEXP(ii) = compexp_econ(DEXP(ii))  ;
    cap_cost = cap_cost + CCMP(ii).cmpexp_cost.COST + DEXP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(CCMP(ii).cmpexp_cost, Nsens) + cost_sens(DEXP(ii).cmpexp_cost, Nsens) ;
end

for ii = 1 : Ne_ch
    CEXP(ii) = compexp_econ(CEXP(ii))  ;
    DCMP(ii) = compexp_econ(DCMP(ii))  ;
    cap_cost = cap_cost + CEXP(ii).cmpexp_cost.COST + DCMP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(CEXP(ii).cmpexp_cost, Nsens) + cost_sens(DCMP(ii).cmpexp_cost, Nsens) ;
end

% Hot tank cost and hot fluid cost
fluidH(1).cost = 1.0 ; % Specify in input file in class constructor
fluidC(1).cost = 1.0 ; % Specify in input file in class constructor
for ii = 1 : Nhot
   HT(ii) = tank_cost(HT(ii)) ;  
   HT(ii) = fld_cost(HT(ii),fluidH(ii).cost) ; 
   cap_cost = cap_cost + HT(ii).tankA_cost.COST + HT(ii).tankB_cost.COST + HT(ii).fluid_cost.COST ;
   cap_sens = cap_sens + cost_sens(HT(ii).tankA_cost, Nsens) + cost_sens(HT(ii).tankB_cost, Nsens) + cost_sens(HT(ii).fluid_cost, Nsens) ;
end

% Cold tank cost and cold fluid cost
for ii = 1 : Ncld
   CT(ii) = tank_cost(CT(ii)) ;  
   CT(ii) = fld_cost(CT(ii),fluidC(ii).cost) ;  
   cap_cost = cap_cost + CT(ii).tankA_cost.COST + CT(ii).tankB_cost.COST + CT(ii).fluid_cost.COST ;
   cap_sens = cap_sens + cost_sens(CT(ii).tankA_cost, Nsens) + cost_sens(CT(ii).tankB_cost, Nsens) + cost_sens(CT(ii).fluid_cost, Nsens) ;
end

cap_cost_std = std(cap_sens) ;
cap_cost_lo  = cap_cost - cap_cost_std ;
cap_cost_hi  = cap_cost + cap_cost_std ;

% Move this kind of stuff into an 'ECONOMICS' class
pow = W_out_dis/t_dis/1e3 ;
en  = W_out_dis/(1e3*3600) ;

cap_cost_pow = cap_cost / pow ; % Cost, $ / kW
cap_cost_en  = cap_cost / en ;  % Cost, $ / kWh