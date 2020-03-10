% Script that calculates the cost of each component and subsequently the
% total cost of the system, as well as other economic metrics (possibly
% ...)

% Some input variables - move these to an input file?
% Have a structure called Cdata
Cdata.lifetime    = 25.0 ;      % Lifetime
Cdata.price       = 0.06 ;      % Electricity price - dollars per kWhe. Lazard uses 0.033
Cdata.inflation   = 0.025;      % Inflation
Cdata.irr         = 0.10 ;      % Internal Rate of Return
Cdata.debt_frac   = 0.60 ;      % Project debt fraction - SAM is 0.60
Cdata.debt_IR     = 0.08 ;      % Debt interest rate
Cdata.tax_rate    = 0.40 ;      % Tax rate
Cdata.deprec      = [0.20 0.32 0.20 0.14 0.14]  ; % Depreciation[0.20 0.32 0.192 0.1152 0.1152 0.0576] ;
Cdata.annual_cost = [1.0 0.0 0.]; % Capital cost incurred in which years [0.80 0.10 0.10] ;
Cdata.construc_IR = 0.0 ;         % Construction interest rate <- new assumption 25/1/17 to make CFF =1. SAM value -> % 0.08 ;
Cdata.OnM         = 0.0225 ;      % Operations and maintenance cost as a fraction of total capital cost - see Georgiou et al 2018
Cdata.conting     = 0.07 ;        % Contingency
Cdata.indirect    = 0.25 ;        % Indirect costs

% Make array of chemical engineering cost indices
CEind = create_CEindex() ;

cap_cost = 0 ;
Nsens    = 10000 ; % How many points to take from distribution for sensitivity analysis
cap_sens = zeros(Nsens,1) ;

% Retrofit? Then don't pay for steam turbine or hot storage system        
Lretro = true ; 

% Compressors and expanders
for ii = 1 : length(CCMP)
    CCMP(ii) = compexp_econ(CCMP(ii), CEind, false, 0)  ;
    cap_cost = cap_cost + CCMP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(CCMP(ii).cmpexp_cost, Nsens) ;
end
for ii = 1 : length(DEXP)
    DEXP(ii) = compexp_econ(DEXP(ii), CEind, false, 0)  ;
    cap_cost = cap_cost + DEXP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(DEXP(ii).cmpexp_cost, Nsens) ;
end
for ii = 1 : length(CEXP)
    CEXP(ii) = compexp_econ(CEXP(ii), CEind, false, 0)  ;
    cap_cost = cap_cost + CEXP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(CEXP(ii).cmpexp_cost, Nsens) ;
end
for ii = 1 : length(DCMP)
    DCMP(ii) = compexp_econ(DCMP(ii), CEind, false, 0)  ;
    cap_cost = cap_cost + DCMP(ii).cmpexp_cost.COST ;
    cap_sens = cap_sens + cost_sens(DCMP(ii).cmpexp_cost, Nsens) ;
end

if Load.mode == 4 || Load.mode == 5 || Load.mode == 6
    %Recompressor
    if Lrcmp
        RCMP     = compexp_econ(RCMP, CEind, false, 0) ;
        cap_cost = cap_cost + RCMP.cmpexp_cost.COST ;
        cap_sens = cap_sens + cost_sens(RCMP.cmpexp_cost, Nsens) ;
    end
end

% Motor-generator. Assume this is just to provide the net work (i.e. don't
% have a motor on the compressor and a separate generator on the expander)
powIN  = W_in_chg/t_chg/1e3 ;
GEN.gen_cost = econ_class(1, 0.2, 5, 0.2) ;
GEN.gen_cost.COST = 1.85e6 * (powIN / 1.18e4)^0.94 ; % Really need to make a class .... just for this
cap_cost = cap_cost + GEN.gen_cost.COST ;
cap_sens = cap_sens + cost_sens(GEN.gen_cost, Nsens) ;

% Heat exchangers
for ii = 1 : numel(HX)
   HX(ii)   = HX_cost(HX(ii), CEind) ;
   cap_cost = cap_cost + HX(ii).hx_cost.COST ;
   cap_sens = cap_sens + cost_sens(HX(ii).hx_cost, Nsens) ;
end

% Hot tank cost and hot fluid cost
for ii = 1 : Nhot
    fluidH(ii).cost = 1.0 ; % Specify in input file in class constructor
    if Lretro && Load.mode == 3
        HT(ii).tankA_cost.COST = 0.01 ;
        HT(ii).tankB_cost.COST = 0.01 ;
        HT(ii).fluid_cost.COST = 0.01 ;
    elseif Lretro && Load.mode == 6 && ii == 1
        HT(ii).tankA_cost.COST = 0.01 ;
        HT(ii).tankB_cost.COST = 0.01 ;
        HT(ii).fluid_cost.COST = 0.01 ;
    else
        HT(ii) = tank_cost(HT(ii), CEind) ;
        HT(ii) = fld_cost(HT(ii),fluidH(ii).cost, CEind) ;
    end
    cap_cost = cap_cost + HT(ii).tankA_cost.COST + HT(ii).tankB_cost.COST + HT(ii).fluid_cost.COST ;
    cap_sens = cap_sens + cost_sens(HT(ii).tankA_cost, Nsens) + cost_sens(HT(ii).tankB_cost, Nsens) + cost_sens(HT(ii).fluid_cost, Nsens) ;
end

% Cold tank cost and cold fluid cost
for ii = 1 : Ncld
   fluidC(ii).cost = 1.0 ; % Specify in input file in class constructor
   CT(ii) = tank_cost(CT(ii), CEind) ;  
   CT(ii) = fld_cost(CT(ii),fluidC(ii).cost, CEind) ;  
   cap_cost = cap_cost + CT(ii).tankA_cost.COST + CT(ii).tankB_cost.COST + CT(ii).fluid_cost.COST ;
   cap_sens = cap_sens + cost_sens(CT(ii).tankA_cost, Nsens) + cost_sens(CT(ii).tankB_cost, Nsens) + cost_sens(CT(ii).fluid_cost, Nsens) ;
end

cap_cost           = cap_cost * (1 + Cdata.conting) * (1 + Cdata.indirect) ; 
Cdata.cap_sens     = cap_sens .* (1 + Cdata.conting) .* (1 + Cdata.indirect) ; 
Cdata.cap_cost     = cap_cost ;
Cdata.cap_costM    = mean(Cdata.cap_sens) ;
Cdata.cap_costSD   = std(Cdata.cap_sens) ;
Cdata.cap_cost_lo  = Cdata.cap_costM - Cdata.cap_costSD ;
Cdata.cap_cost_hi  = Cdata.cap_costM + Cdata.cap_costSD ;

pow  = W_out_dis/t_dis/1e3 ;
Wout = W_out_dis/(1e3*3600) ;
Win  = W_in_chg/(1e3*3600) ;

Cdata.cap_cost_pow = Cdata.cap_costM / pow ; % Cost, $ / kW
Cdata.cap_cost_en  = Cdata.cap_costM / Wout ;  % Cost, $ / kWh

% ** Now calculate metrics such as the Levelized Cost of Storage **
% Use FCR method - see Short and Packey
Cdata = calc_fcr(Cdata) ;

% Calculate the LCOS
Cdata = calc_lcos(Cdata, Win, Wout, Load.time, Nsens) ;

% Write out some results
switch Load.mode
    case {0,1,3,4, 5, 6}
        fprintf(1,'ECONOMIC RESULTS:\n');
        fprintf(1,'Capital cost:                  %8.1f M$\n',Cdata.cap_costM/1e6);
        fprintf(1,'     Standard deviation:       %8.1f M$\n\n',Cdata.cap_costSD/1e6);
        fprintf(1,'Cost per unit power:           %8.1f $/kW-e\n',Cdata.cap_cost_pow);
        fprintf(1,'Cost per unit energy:          %8.1f $/kWh-e\n\n',Cdata.cap_cost_en);
        fprintf(1,'Levelised cost of storage:     %8.3f $/kWh-e\n',Cdata.lcosM);
        fprintf(1,'     Standard deviation:       %8.3f $/kWh-e\n\n',Cdata.lcosSD);
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

