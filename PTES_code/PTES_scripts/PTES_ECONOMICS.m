% Script that calculates the cost of each component and subsequently the
% total cost of the system, as well as other economic metrics (possibly
% ...)
cap_cost = 0 ;

% Compressors and expanders
for ii = 1 : Nc_ch
    CCMP(ii) = compexp_econ(CCMP(ii))  ;
    DEXP(ii) = compexp_econ(DEXP(ii))  ;
    cap_cost = cap_cost + CCMP(ii).COST + DEXP(ii).COST ;
end

for ii = 1 : Ne_ch
    CEXP(ii) = compexp_econ(CEXP(ii))  ;
    DCMP(ii) = compexp_econ(DCMP(ii))  ;
    cap_cost = cap_cost + CEXP(ii).COST + DCMP(ii).COST ;
end

% Hot tank cost and hot fluid cost
fluidH(1).cost = 1.0 ; % Specify in input file in class constructor
fluidC(1).cost = 1.0 ; % Specify in input file in class constructor
for ii = 1 : Nhot
   HT(ii) = tank_cost(HT(ii),1) ;  
   HT(ii) = fld_cost(HT(ii),fluidH(ii).cost) ; 
   cap_cost = cap_cost + HT(ii).tank_costA + HT(ii).tank_costB + HT(ii).fluid_cost ;
end

% Cold tank cost and cold fluid cost
for ii = 1 : Ncld
   CT(ii) = tank_cost(CT(ii),1) ;  
   CT(ii) = fld_cost(CT(ii),fluidC(ii).cost) ;  
   cap_cost = cap_cost + CT(ii).tank_costA + CT(ii).tank_costB + CT(ii).fluid_cost ;
end

% Move this kind of stuff into an 'ECONOMICS' class
pow = W_out_dis/t_dis/1e3 ;
en  = W_out_dis/(1e3*3600) ;

cap_cost_pow = cap_cost / pow ; % Cost, $ / kW
cap_cost_en  = cap_cost / en ;  % Cost, $ / kWh