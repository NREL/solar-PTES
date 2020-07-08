% Script that calculates the cost of each component and subsequently the
% total cost of the system, as well as other economic metrics such as the
% levelized cost of storage.

% FALSE: Calculate the cost using one set of cost correlation, then
% calculate the sensitivity assuming each cost is normally distributed
% TRUE: Calculate the cost numerous times using different combinations of
% different cost correlations
Lsuper = 1 ;

% Some input variables - move these to an input file?
price = [0.033,0.025,0.06] ;
life  = [25,30,35];
OnM   = [0.0225,0.01,0.05];
cont  = [0,0.1,0.5] ;

% Have a structure called Cdata
Cdata.lifetime    = life(1) ;      % Lifetime
Cdata.price       = price(1) ;      % Electricity price - dollars per kWhe. Lazard uses 0.033, ARPA-E uses 0.025. 0.06 is a value that I've used in the past.
Cdata.inflation   = 0.025;      % Inflation
Cdata.irr         = 0.10 ;      % Internal Rate of Return
Cdata.debt_frac   = 0.60 ;      % Project debt fraction - SAM is 0.60
Cdata.debt_IR     = 0.08 ;      % Debt interest rate
Cdata.tax_rate    = 0.40 ;      % Tax rate
Cdata.deprec      = [0.20 0.32 0.20 0.14 0.14]  ; % Depreciation[0.20 0.32 0.192 0.1152 0.1152 0.0576] ;
Cdata.annual_cost = [1.0 0.0 0.]; % Capital cost incurred in which years [0.80 0.10 0.10] ;
Cdata.construc_IR = 0.0 ;         % Construction interest rate <- new assumption 25/1/17 to make CFF =1. SAM value -> % 0.08 ;
Cdata.OnM         = OnM(1) ;      % Operations and maintenance cost as a fraction of total capital cost - see Georgiou et al 2018
Cdata.conting     = cont(1);%0.07 ;        % Contingency
Cdata.indirect    = 0;%0.25 ;        % Indirect costs


if Lsuper
    Nsens    = 1 ;      % How many points to take from distribution for sensitivity analysis
    Ncomb    = 1000 ;   % How many combinations of cost correlations?
    
    costMAT      = zeros(Ncomb,1) ;
    cost_enMAT   = zeros(Ncomb,1) ;
    cost_powMAT  = zeros(Ncomb,1) ;
    lcosMAT      = zeros(Ncomb,1) ;
    
    compMAT      = zeros(Ncomb,13) ; % Each column is for a different component
    
    rng('shuffle') % Shuffle the random number generator
else
    Nsens    = 10000 ;  % How many points to take from distribution for sensitivity analysis
    Ncomb    = 1 ;      % How many combinations of cost correlations?
end


for jj = 1 : Ncomb
    
% For each component select a random cost mode
if Lsuper
    for ii = 1 : length(CCMP)
        CCMP(ii).cmpexp_cost.cost_mode = CCMPmode(randi(length(CCMPmode))) ;    
    end
    for ii = 1 : length(DCMP)
        DCMP(ii).cmpexp_cost.cost_mode = DCMPmode(randi(length(DCMPmode))) ;    
    end
    for ii = 1 : length(CEXP)
        CEXP(ii).cmpexp_cost.cost_mode = CEXPmode(randi(length(CEXPmode))) ;    
    end
    for ii = 1 : length(DEXP)
        DEXP(ii).cmpexp_cost.cost_mode = DEXPmode(randi(length(DEXPmode))) ;    
    end
    for ii = 1 : length(CPMP)
        CPMP(ii).cmpexp_cost.cost_mode = PMPmode(randi(length(PMPmode))) ;    
    end
    for ii = 1 : length(DPMP)
        DPMP(ii).cmpexp_cost.cost_mode = PMPmode(randi(length(PMPmode))) ;    
    end
    for ii = 1 : length(CFAN)
        CFAN(ii).cmpexp_cost.cost_mode = FANmode(randi(length(FANmode))) ;    
    end
    for ii = 1 : length(DFAN)
        DFAN(ii).cmpexp_cost.cost_mode = FANmode(randi(length(FANmode))) ;    
    end
    for ii = 1 : numel(HX)
        if strcmp(HX(ii).name,'hot')
            HX(ii).hx_cost.cost_mode = hotHXmode(randi(length(hotHXmode))) ;
        elseif strcmp(HX(ii).name,'cold')
            HX(ii).hx_cost.cost_mode = cldHXmode(randi(length(cldHXmode))) ;
        elseif strcmp(HX(ii).name,'regen')
            HX(ii).hx_cost.cost_mode = rcpHXmode(randi(length(rcpHXmode))) ;
        elseif strcmp(HX(ii).name,'rej')
            HX(ii).hx_cost.cost_mode = rejHXmode(randi(length(rejHXmode))) ;
        end
    end
    for ii = 1 : numel(HT)
        nH = length(HTmode.tankmode) ;
        HT(ii).tankA_cost.cost_mode = HTmode.tankmode(randi(nH)) ;
        HT(ii).tankB_cost.cost_mode = HT(ii).tankA_cost.cost_mode ;
        
        nF = length(HTmode.fld_cost) ;
        HT(ii).fluid_cost.cost_mode = HTmode.fld_cost(randi(nF)) ;
        
        nI = length(HTmode.ins_cost) ;
        HT(ii).insA_cost.cost_mode = HTmode.ins_cost(randi(nI)) ;
        HT(ii).insB_cost.cost_mode = HTmode.ins_cost(randi(nI)) ;
    end
    for ii = 1 : numel(CT)
        nC = length(CTmode.tankmode) ;
        CT(ii).tankA_cost.cost_mode = CTmode.tankmode(randi(nC)) ;
        CT(ii).tankB_cost.cost_mode = CT(ii).tankA_cost.cost_mode ;
        
        nF = length(CTmode.fld_cost) ;
        CT(ii).fluid_cost.cost_mode = CTmode.fld_cost(randi(nF)) ;
        
        nI = length(CTmode.ins_cost) ;
        CT(ii).insA_cost.cost_mode = CTmode.ins_cost(randi(nI)) ;
        CT(ii).insB_cost.cost_mode = CTmode.ins_cost(randi(nI)) ;
    end
    
    GEN.gen_cost.cost_mode = GENmode(randi(length(GENmode))) ;
    
    Cdata.lifetime = life(randi(length(life))) ;
    Cdata.price    = price(randi(length(price))) ;
    Cdata.OnM      = OnM(randi(length(OnM))) ;
    Cdata.conting  = cont(randi(length(cont))) ;
    
end

    


% Make array of chemical engineering cost indices
CEind = create_CEindex() ;

cap_cost = 0 ;
cap_sens = zeros(Nsens,1) ;

% Retrofit? Then don't pay for steam turbine or hot storage system        
Lretro = true ; 

% Compressors and expanders
for ii = 1 : length(CCMP)
    if any(Load.mode == [2,7])
        CCMP(ii).cmpexp_cost.COST = 0.01 ;
    else
        CCMP(ii) = compexp_econ(CCMP(ii), CEind, gas)  ;
    end
    cap_cost = cap_cost + CCMP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(CCMP(ii).cmpexp_cost, Nsens) ;
end
for ii = 1 : length(DEXP)
    if Load.mode == 1 || (Lretro && Load.mode == 3)
        DEXP(ii).cmpexp_cost.COST = 0.01 ;
    else
        DEXP(ii) = compexp_econ(DEXP(ii), CEind, gas)  ;
    end
    cap_cost = cap_cost + DEXP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(DEXP(ii).cmpexp_cost, Nsens) ;
end
for ii = 1 : length(CEXP)
    if any(Load.mode == [2,7])
        CEXP(ii).cmpexp_cost.COST = 0.01 ;
    else
        CEXP(ii) = compexp_econ(CEXP(ii), CEind, gas)  ;
    end
    cap_cost = cap_cost + CEXP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(CEXP(ii).cmpexp_cost, Nsens) ;
end
for ii = 1 : length(DCMP)
    if Load.mode == 1 || (Lretro && Load.mode == 3)
        DCMP(ii).cmpexp_cost.COST = 0.01 ;
    else
        DCMP(ii) = compexp_econ(DCMP(ii), CEind, gas)  ;
    end
    cap_cost = cap_cost + DCMP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(DCMP(ii).cmpexp_cost, Nsens) ;
end

if Load.mode == 4 || Load.mode == 5 || Load.mode == 6
    %Recompressor
    if Lrcmp
        RCMP     = compexp_econ(RCMP, CEind, gas) ;
        cap_cost = cap_cost + RCMP.cmpexp_cost.COST ;
        cap_sens = cap_sens + cost_sens(RCMP.cmpexp_cost, Nsens) ;
    end
end

% Pumps
for ii = 1 : length(CPMP)
    if CPMP(ii).W0 == 0
        CPMP(ii).cmpexp_cost.COST = 0.01 ;
    else
        CPMP(ii) = compexp_econ(CPMP(ii), CEind, gas)  ;
    end
    cap_cost = cap_cost + CPMP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(CPMP(ii).cmpexp_cost, Nsens) ;
end
for ii = 1 : length(DPMP)
    if DPMP(ii).W0 == 0
        DPMP(ii).cmpexp_cost.COST = 0.01 ;
    else
        DPMP(ii) = compexp_econ(DPMP(ii), CEind, gas)  ;
    end
    cap_cost = cap_cost + DPMP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(DPMP(ii).cmpexp_cost, Nsens) ;
end

% FANS
for ii = 1 : length(CFAN)
    if CFAN(ii).W0 == 0
        CFAN(ii).cmpexp_cost.COST = 0.01 ;
    else
        CFAN(ii) = compexp_econ(CFAN(ii), CEind, gas)  ;
    end
    cap_cost = cap_cost + CFAN(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(CFAN(ii).cmpexp_cost, Nsens) ;
end
for ii = 1 : length(DFAN)
    if DFAN(ii).W0 == 0 || (Lretro && Load.mode == 3)
        DFAN(ii).cmpexp_cost.COST = 0.01 ;
    else
        DFAN(ii) = compexp_econ(DFAN(ii), CEind, gas)  ;
    end
    cap_cost = cap_cost + DFAN(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(DFAN(ii).cmpexp_cost, Nsens) ;
end


% Motor-generator. Assume this is just to provide the net work (i.e. don't
% have a motor on the compressor and a separate generator on the expander)
% This needs to be based on design values of the compressors and expanders
if Load.mode == 7
    GEN.gen_cost.cost_mode = 0 ;
end
GEN = gen_econ(GEN,CEind) ;
cap_cost = cap_cost + GEN.gen_cost.COST ;
cap_sens = cap_sens + cost_sens(GEN.gen_cost, Nsens) ;
    
   
% Electric heater
if any(Load.mode == [2,7])
   Leh = true ; % Is there an electric heater to charge the hot tanks?
   if Leh
      % Heat from hot tanks
      QH = ((HT.A(2).H - HT.A(end).H) + (HT.B(2).H - HT.B(end).H))/(t_dis*1e6) ; % This equals heat added by electric heater
      EH.eh_cost = econ_class(1, 0.2, 5, 0.2) ;
      EH.eh_cost.COST = 75e3 * QH ^ 0.9 ; % Correlation from Benato et al. 2017
      %EH.eh_cost.COST = 475e3 * QH ; 
      cap_cost = cap_cost + EH.eh_cost.COST ;
      cap_sens = cap_sens + cost_sens(EH.eh_cost, Nsens) ;
   end
end

% Heat exchangers
% If the heat exchanger was employed with the 'eff' or 'DT' modes, the
% required geometry is computed now
for ii = 1 : numel(HX)
    if any(strcmp(HX(ii).model,{'eff','DT'})) && (~HX(ii).Lgeom_set)
        HX(ii)   = hex_set_geom(HX(ii)); %#ok<*SAGROW>
    end
    HX(ii)   = HX_cost(HX(ii), CEind) ;
    cap_cost = cap_cost + HX(ii).hx_cost.COST ;
    cap_sens = cap_sens + cost_sens(HX(ii).hx_cost, Nsens) ;
end

% Hot tank cost and hot fluid cost
for ii = 1 : Nhot
    if Lretro && any(Load.mode == [3,7])
        HT(ii).tankA_cost.COST = 0.01 ;
        HT(ii).tankB_cost.COST = 0.01 ;
        HT(ii).fluid_cost.COST = 0.01 ;
        HT(ii).insA_cost.COST = 0.01 ;
        HT(ii).insB_cost.COST = 0.01 ;
    elseif Lretro && Load.mode == 6 && ii == 1
        HT(ii).tankA_cost.COST = 0.01 ;
        HT(ii).tankB_cost.COST = 0.01 ;
        HT(ii).fluid_cost.COST = 0.01 ;
    else
        HT(ii) = tank_cost(HT(ii), CEind) ;
        HT(ii) = ins_cost(HT(ii), CEind) ;
        HT(ii) = fld_cost(HT(ii), CEind) ;
    end
    cap_cost = cap_cost + HT(ii).tankA_cost.COST + HT(ii).tankB_cost.COST + HT(ii).fluid_cost.COST + HT(ii).insA_cost.COST + HT(ii).insB_cost.COST;
    cap_sens = cap_sens + cost_sens(HT(ii).tankA_cost, Nsens) + cost_sens(HT(ii).tankB_cost, Nsens) + cost_sens(HT(ii).fluid_cost, Nsens) + cost_sens(HT(ii).insA_cost, Nsens) + cost_sens(HT(ii).insB_cost, Nsens);
end

% Cold tank cost and cold fluid cost
for ii = 1 : Ncld
   if any(Load.mode == [2,7])
       CT(ii).tankA_cost.COST = 0.01 ;
       CT(ii).tankB_cost.COST = 0.01 ;
       CT(ii).fluid_cost.COST = 0.01 ;
   else
       CT(ii) = tank_cost(CT(ii), CEind) ;
       CT(ii) = ins_cost(CT(ii), CEind) ;
       CT(ii) = fld_cost(CT(ii), CEind) ;
   end
   cap_cost = cap_cost + CT(ii).tankA_cost.COST + CT(ii).tankB_cost.COST + CT(ii).fluid_cost.COST + CT(ii).insA_cost.COST + CT(ii).insB_cost.COST ;
   cap_sens = cap_sens + cost_sens(CT(ii).tankA_cost, Nsens) + cost_sens(CT(ii).tankB_cost, Nsens) + cost_sens(CT(ii).fluid_cost, Nsens) + cost_sens(HT(ii).insA_cost, Nsens) + cost_sens(HT(ii).insB_cost, Nsens);
end

cap_cost           = cap_cost * (1 + Cdata.conting) * (1 + Cdata.indirect) ; 
Cdata.cap_sens     = cap_sens .* (1 + Cdata.conting) .* (1 + Cdata.indirect) ; 
Cdata.cap_cost     = cap_cost ;
Cdata.cap_costM    = mean(Cdata.cap_sens) ;
Cdata.cap_costSD   = std(Cdata.cap_sens) ;
Cdata.cap_cost_lo  = Cdata.cap_costM - Cdata.cap_costSD ;
Cdata.cap_cost_hi  = Cdata.cap_costM + Cdata.cap_costSD ;

switch Load.mode
    case {0,3,4,6}
        pow  = W_out_dis/t_dis/1e3 ;
        Wout = W_out_dis/(1e3*3600) ;
        Win  = -W_in_chg/(1e3*3600) ;

        Cdata.cap_cost_pow = Cdata.cap_costM / pow ; % Cost, $ / kW
        Cdata.cap_cost_en  = Cdata.cap_costM / Wout ;  % Cost, $ / kWh

        % Use FCR method - see Short and Packey
        Cdata = calc_fcr(Cdata) ;

        % Calculate the LCOS
        Cdata = calc_lcos(Cdata, Win, Wout, Load.time, Nsens) ;

    case 1
        pow  = -W_in_chg/t_chg/1e3 ;
        Win  = -W_in_chg/(1e3*3600) ;
        
        Cdata.cap_cost_pow = Cdata.cap_costM / pow ; % Cost, $ / kW
        Cdata.cap_cost_en  = Cdata.cap_costM / Win ;  % Cost, $ / kWh
        
    case {2,5,7}
        pow  = W_out_dis/t_dis/1e3 ;
        Wout = W_out_dis/(1e3*3600) ;
        
        Cdata.cap_cost_pow = Cdata.cap_costM / pow ; % Cost, $ / kW
        Cdata.cap_cost_en  = Cdata.cap_costM / Wout ;  % Cost, $ / kWh

end

if Lsuper
   costMAT(jj) = Cdata.cap_cost ;
   cost_enMAT(jj) = Cdata.cap_cost_en ;
   cost_powMAT(jj) = Cdata.cap_cost_pow ;
   lcosMAT(jj) = Cdata.lcosM ;
   
   % Assign component costs to the cost matrix
   % Compressors
   for ii = 1 : length(CCMP)
       compMAT(jj,1) = compMAT(jj,1) + CCMP(ii).cmpexp_cost.COST ;
   end
   for ii = 1 : length(DCMP)
       compMAT(jj,1) = compMAT(jj,1) + DCMP(ii).cmpexp_cost.COST ;
   end
   
   % Expanders
   for ii = 1 : length(CEXP)
       compMAT(jj,2) = compMAT(jj,2) + CEXP(ii).cmpexp_cost.COST ;
   end
   for ii = 1 : length(DEXP)
       compMAT(jj,2) = compMAT(jj,2) + DEXP(ii).cmpexp_cost.COST ;
   end
   
   % Pumps
   for ii = 1 : length(CPMP)
       compMAT(jj,3) = compMAT(jj,3) + CPMP(ii).cmpexp_cost.COST ;
   end
   for ii = 1 : length(DPMP)
       compMAT(jj,3) = compMAT(jj,3) + DPMP(ii).cmpexp_cost.COST ;
   end
   
   % Fans
   for ii = 1 : length(CFAN)
       compMAT(jj,4) = compMAT(jj,4) + CFAN(ii).cmpexp_cost.COST ;
   end
   for ii = 1 : length(DFAN)
       compMAT(jj,4) = compMAT(jj,4) + DFAN(ii).cmpexp_cost.COST ;
   end
   
   
   % Heat exchangers and rejection
   for ii = 1 : length(HX)
       switch HX(ii).name
           case 'hot'
               compMAT(jj,5) = compMAT(jj,5) + HX(ii).hx_cost.COST ;
           case 'cold'
               compMAT(jj,6) = compMAT(jj,6) + HX(ii).hx_cost.COST ;
           case 'regen'
               compMAT(jj,7) = compMAT(jj,7) + HX(ii).hx_cost.COST ;
           case 'rej'
               compMAT(jj,8) = compMAT(jj,8) + HX(ii).hx_cost.COST ;
       end
   end
   
   % Tanks, fluid, insulation
   for ii = 1 : length(HT)
       compMAT(jj,9)  = compMAT(jj,9)  + HT(ii).tankA_cost.COST + HT(ii).tankB_cost.COST ;
       compMAT(jj,10) = compMAT(jj,10) + HT(ii).insA_cost.COST  + HT(ii).insB_cost.COST ;
       compMAT(jj,11) = compMAT(jj,11) + HT(ii).fluid_cost.COST;
   end
   
   for ii = 1 : length(CT)
       compMAT(jj,9)  = compMAT(jj,9)  + CT(ii).tankA_cost.COST + CT(ii).tankB_cost.COST ;
       compMAT(jj,10) = compMAT(jj,10) + CT(ii).insA_cost.COST  + CT(ii).insB_cost.COST ;
       compMAT(jj,12) = compMAT(jj,12) + CT(ii).fluid_cost.COST;
   end
   
   % Generator/motor
   compMAT(jj,13) = compMAT(jj,13) + GEN.gen_cost.COST;
end

end

if Lsuper
    % Plot box and whisker plot of component costs
    PLOT_BOX(compMAT) ;
    
    Cdata.cap_costM    = mean(costMAT) ;
    Cdata.cap_costSD   = std(costMAT) ;
    Cdata.cap_cost_pow = mean(cost_powMAT) ;
    Cdata.cap_cost_en  = mean(cost_enMAT) ;
    Cdata.lcosM        = mean(lcosMAT) ;
    Cdata.lcosSD       = std(lcosMAT) ;

end
% Write out some results
switch Load.mode
    case {0,3,4,6}
        fprintf(1,'ECONOMIC RESULTS:\n');
        fprintf(1,'Capital cost:                  %8.1f M$\n',Cdata.cap_costM/1e6);
        fprintf(1,'     Standard deviation:       %8.1f M$\n\n',Cdata.cap_costSD/1e6);
        fprintf(1,'Cost per unit power:           %8.1f $/kW-e\n',Cdata.cap_cost_pow);
        fprintf(1,'Cost per unit energy:          %8.1f $/kWh-e\n\n',Cdata.cap_cost_en);
        fprintf(1,'Levelised cost of storage:     %8.3f $/kWh-e\n',Cdata.lcosM);
        fprintf(1,'     Standard deviation:       %8.3f $/kWh-e\n\n',Cdata.lcosSD);
    case {1,2,5,7}
        fprintf(1,'ECONOMIC RESULTS:\n');
        fprintf(1,'Capital cost:                  %8.1f M$\n',Cdata.cap_costM/1e6);
        fprintf(1,'     Standard deviation:       %8.1f M$\n\n',Cdata.cap_costSD/1e6);
        fprintf(1,'Cost per unit power:           %8.1f $/kW-e\n',Cdata.cap_cost_pow);
        fprintf(1,'Cost per unit energy:          %8.1f $/kWh-e\n\n',Cdata.cap_cost_en);
end

% Calculate the fixed charge rate and other economic factors
% Calculations based on the model in SAM
function obj = calc_fcr(obj)
    obj.RROE = (1 + obj.irr) / (1 + obj.inflation) - 1.0 ;
    obj.RINT = (1 + obj.debt_IR) / (1 + obj.inflation) - 1.0 ;
    obj.WACC = ((1 + obj.inflation) * (1 + obj.RROE) - 1) * (1 - obj.debt_frac) + 1 ;
    obj.WACC = obj.WACC + obj.debt_frac * ((1 + obj.RINT) * (1 + obj.inflation) - 1) * (1 - obj.tax_rate) ;
    obj.WACC = obj.WACC / (1 + obj.inflation) - 1.0 ;

    obj.CRF  = obj.WACC / (1.0 - (1.0 / (1 + obj.WACC)^obj.lifetime) ) ;

    obj.PVDEP = 0.0 ;
    for i = 1:numel(obj.deprec)
        obj.PVDEP = obj.PVDEP + obj.deprec(i) / (((1 + obj.WACC) * (1 + obj.inflation))^i );
    end

    obj.PFF   = (1 - obj.tax_rate * obj.PVDEP) / (1.0 - obj.tax_rate) ;

    obj.CFF   = 0.0 ;
    for i = 1:numel(obj.annual_cost)
        obj.CFF = obj.CFF + (((1 + obj.construc_IR)^(i-0.5) - 1) * (1 - obj.tax_rate) + 1) * obj.annual_cost(i) ;
    end

    obj.FCR   = obj.CRF * obj.PFF * obj.CFF ;
end


% Calculate the levelised cost of energy
function obj = calc_lcos(obj, Win, Wout, times, Nsens)

    % Calculate the total time that elapsed in this run
    totT = sum(times) ; % seconds 
    Ncyc = 8760 * 3600 / totT ; % How many cycles in one year
    
    totWin  = Win * Ncyc ;   % Total work input in one year
    totWout = Wout * Ncyc ;  % Total work output in one year
    
    % Make distributions of OnM, FCR, elec price
    % Create classes first
    OnMclass  = econ_class(0, 0.2, 2, 0.3) ;
    FCRclass  = econ_class(0, 0.2, 2, 0.3) ;
    ELECclass = econ_class(0, 0.2, 2, 0.3) ;
    
    OnMclass.COST = obj.OnM ;
    FCRclass.COST = obj.FCR ;
    ELECclass.COST = obj.price ;
    
    OnM_sens   = cost_sens(OnMclass, Nsens) ;
    FCR        = cost_sens(FCRclass, Nsens) ;
    price_sens = cost_sens(ELECclass, Nsens) ;
    
    % Total O&M cost
    OnM = OnM_sens .* obj.cap_sens ;
    
    % Total cost of electricity for charging
    elec = price_sens .* totWin ;
    
    % Levelised cost of storage
    obj.lcos_sens = (obj.cap_sens .* FCR + OnM + elec) ./ totWout ;
    
    obj.lcosM  = mean(obj.lcos_sens) ;
    obj.lcosSD = std(obj.lcos_sens) ;
    obj.lcos_lo = obj.lcosM - obj.lcosSD ;
    obj.lcos_hi = obj.lcosM + obj.lcosSD ;

end

function PLOT_BOX(cost)

figure(77)

xlab = {'Compressors','Expanders','Pumps','Fans','Hot HXs','Cold HXs','Recuperators','Heat rejection','Storage tanks','Insulation','Hot fluid','Cold fluid','Motor-generator'} ;
bplot(cost,'nomean','whisker',1) ;
set(gca, 'XTick', 1:numel(cost(1,:)), 'XTickLabel', xlab, 'TickLabelInterpreter', 'latex')
xlim([0 numel(cost(1,:))+1]) ;
xtickangle(45)
ylabel(strcat('Capital cost, \$'))
box on

end

% Chemical engineering cost indices from
% chemengonline.com
function CEindex = create_CEindex()

    CEindex = zeros(2019,1) ;
    CEindex(1947) = 	64.8  ;
    CEindex(1948) = 	70.2  ;
    CEindex(1949) = 	71.4  ;
    CEindex(1950) = 	73.9  ;
    CEindex(1951) = 	80.4  ;
    CEindex(1952) = 	81.3  ;
    CEindex(1953) = 	84.7  ;
    CEindex(1954) = 	86.1  ;
    CEindex(1955) = 	88.3  ;
    CEindex(1956) = 	93.9  ;
    CEindex(1957) = 	98.5  ;
    CEindex(1958) = 	99.7  ;
    CEindex(1959) = 	101.8 ;
    CEindex(1960) = 	102	;
    CEindex(1961) = 	101.5 ;
    CEindex(1962) = 	102 ;
    CEindex(1963) = 	102.4 ;
    CEindex(1964) = 	103.3 ;
    CEindex(1965) = 	104.2 ;
    CEindex(1966) = 	107.2 ;
    CEindex(1967) = 	109.7 ;
    CEindex(1968) = 	113.7 ;
    CEindex(1969) = 	119 ;
    CEindex(1970) = 	125.7 ;
    CEindex(1971) = 	132.2 ;
    CEindex(1972) = 	137.2 ;
    CEindex(1973) = 	144.1 ;
    CEindex(1974) = 	165.4 ;
    CEindex(1975) = 	182.3 ;
    CEindex(1976) = 	192 ;
    CEindex(1977) = 	204.1 ;
    CEindex(1978) = 	218.8 ;
    CEindex(1979) = 	238.7 ;
    CEindex(1980) = 	261.1 ;
    CEindex(1981) = 	297 ;
    CEindex(1982) = 	314 ;
    CEindex(1983) = 	317 ;
    CEindex(1984) = 	322.6 ;
    CEindex(1985) = 	325.3 ;
    CEindex(1986) = 	318.3 ;
    CEindex(1987) = 	323.7 ;
    CEindex(1988) = 	342.4 ;
    CEindex(1989) = 	355.5 ;
    CEindex(1990) = 	357.6 ;
    CEindex(1991) = 	361.3 ;
    CEindex(1992) = 	358.2 ;
    CEindex(1993) = 	359.2 ;
    CEindex(1994) = 	368.1 ;
    CEindex(1995) = 	381.1 ;
    CEindex(1996) = 	381.7 ;
    CEindex(1997) = 	386.5 ;
    CEindex(1998) = 	389.5 ;
    CEindex(1999) = 	390.6 ;
    CEindex(2000) = 	394.1 ;
    CEindex(2001) = 	394.3 ;
    CEindex(2002) = 	395.6 ;
    CEindex(2003) = 	402   ;
    CEindex(2004) = 	444.2 ;
    CEindex(2005) = 	468.2 ;
    CEindex(2006) = 	499.6 ;
    CEindex(2007) = 	525.4 ;
    CEindex(2008) = 	575.4 ;
    CEindex(2009) = 	521.9 ;
    CEindex(2010) = 	550.8 ;
    CEindex(2011) = 	585.7 ;
    CEindex(2012) = 	584.6 ;
    CEindex(2013) = 	567.3 ;
    CEindex(2014) = 	576.1 ;
    CEindex(2015) = 	556.8 ;
    CEindex(2016) = 	541.7 ;
    CEindex(2017) = 	567.5 ;
    CEindex(2018) = 	603.1 ;
    CEindex(2019) = 	619.2 ;

end

