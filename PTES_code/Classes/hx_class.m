classdef hx_class
   properties
       name
       mode
       
       model
       stage_type
       
       eff      % Effectiveness
       ploss    % Pressure loss
       
       UA0      % Conductance, W/K - (design)
       NTU0     % Design NTU
       LMTD0    % Design LMTD
       Hname    % Hot fluid name
       Cname    % Cold fluid name
       
       Qmax     % Max heat transfer - design point
       Qact     % Actual heat transfer - design point
       
       Nsave    % How many load cycles to save data for. All of them, or just two (first charge and first discharge)?
       NX       % Number of points along HX
       
       % These are structures containing T, P, mdot along the HX. 
       % These structures are arrays of length Nsave
       H = stream ; % Hot side
       C = stream ; % Cold side
       
       QS       % Cumulative heat transfer
       AS       % Heat transfer area
       LMTD     % Why not
       Cmin 
       UA
       NTU
       DppH
       DppC
       
       % **
       
       % Geometry
       shape
       L        % Tube length, m
       AfT      % Total flow area, m2
       D1       % Tube diameter, m
       AfR      % Ratio of flow areas (shell/tube)
       
       % More geometry
       t1       % Tube wall thickness
       N1       % Number of tubes
       D2       % Shell-side hydraulic diameter
       G1       % Mass flux, tube-side
       G2       % Mass flux, shell-side
       Af1      % Flow area
       Af2      % Flow area
       A1       % Heat transfer area
       A2       % Heat transfer area
       Vm       % Volume of metal
       
       % Loss factors - one for each time period
       w       % specific work transfer, J/kg
       q       % specific heat transfer, J/kg
       Dh      % specific enthalpy change, J/kg
       sirr    % specific entropy generation, J/kgK
       
       W        % work transfer, W
       Q        % heat transfer, W
       DH       % enthalpy change, W
       Sirr     % entropy generation, W/K
       
       % Costs
       LestA   % Logical - estimate the area (if using effectivenesss/UA method) to enable cost calculations
       hx_cost = econ_class(0,0,0,0) ;
       
   end
   
   methods
       function obj = hx_class(name, stage_type, model, eff, ploss, cost_mode, Ngrid, Nsave, numPeriods)
            obj.name       = name ;
            obj.stage_type = stage_type ;
            obj.model      = model ;
            obj.eff        = eff ;
            obj.ploss      = ploss ;
            
            obj.NX         = Ngrid ;
            
            % Loss data          
            % Two columns. First for hot side. Second for cold side.
            obj.w    = zeros(numPeriods,2) ;
            obj.q    = zeros(numPeriods,2) ;
            obj.Dh   = zeros(numPeriods,2) ;
            obj.sirr = zeros(numPeriods,2) ;
                        
            obj.W    = zeros(numPeriods,2) ;
            obj.Q    = zeros(numPeriods,2) ;
            obj.DH   = zeros(numPeriods,2) ;
            obj.Sirr = zeros(numPeriods,2) ;
            
            % Property data
            %obj.H    = zeros(Nsave, 1) ;
            %obj.C    = zeros(Nsave, 1) ;
                        
            obj.QS    = zeros(Ngrid+1,Nsave) ;
            
            obj.hx_cost = econ_class(cost_mode, 0.2, 5, 0.2) ;
            
       end
        
       
       % Calculate the HX cost
       % Costs come from Q2 report unless otherwise stated
       function [obj] = HX_cost(obj, CEind)
           
           curr = 2019 ; % Current year
           switch obj.hx_cost.cost_mode
               case 0
                   COST = 0 ;
                   
               % ** 
               case 1
                   % Shell and tube. Stainless steel. Eq. HX3.2 
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   COST = 2768 * A ^ 0.573 ;
                   COST = COST * CEind(curr) / CEind(1986) ;
                   
               case 2
                   % Shell and tube. Carbon steel. 10 bar. Eq. HX5
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   COST = 9583.8 + 251.5 * A ;
                   COST = COST * CEind(curr) / CEind(1998) ;
                   
               case 3
                   % Shell and tube from Couper et al 3rd ed 2012 (2009
                   % costs). Have to choose material and pressure range.
                   % For 150 < A < 12e3 sqft
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   % Convert area to square feet
                   A = A * 3.28084^2 ;
                   
                   COST = exp(8.821 - 0.30863 * log(A) + 0.0681 * (log(A))^2) ;
                   
                   type = 'fixed head' ;
                   
                   switch type
                       case 'fixed head'
                           f1   = exp(-1.1156 + 0.0906 * log(A)) ; % For fixed head
                       case 'kettle'
                           f1 = 1.35 ; % For kettle reboiler
                       case 'U-tube'
                           f1   = exp(-0.9816 + 0.083 * log(A)) ; % For U-tube
                   end
                   
                   P = 0.000145038 * max(obj.H(1).pin,obj.C(1).pin) ; % Convert max pressure to psi
                   if P < 100 
                       f2 = 1.0 ;
                   elseif P < 300
                       f2 = 0.7771 + 0.04981 * log(A) ;
                   elseif P < 600
                       f2 = 1.0305 + 0.0714 * log(A) ;
                   else
                       f2 = 1.14 + 0.12088 * log(A) ;
                   end
                   
                   mat = 'SS316' ;
                   
                   switch mat
                       case 'SS316'
                           f3 = 0.86030 + 0.23296 * log(A) ;                            
                       case 'Nickel 200'
                           f3 = 1.5092 + 0.60859 * log(A) ;
                       case 'Inconel 600'
                           f3 = 1.204 + 0.50764 * log(A) ;
                       case 'Monel 400'
                           f3 = 1.2989 + 0.43377 * log(A) ;
                       case 'Ti'
                           f3 = 1.542 + 0.42913 * log(A) ;
                   end
                    
                   COST = COST * f1 * f2 * f3 ;
                   COST = COST * CEind(curr) / CEind(2009) ;
                   
                   
               case 4
                   % Shell and tube from Sieder et al 3rd ed 2017 
                   % Have to choose material and pressure.
                   % For 150 < A < 12e3 sqft
                   % Assumed to use 3/4 or 1 inch O.D. carbon steel tubes,
                   % 20 ft long, on square or triangular pitch, at
                   % pressures up to 100 psig
                 
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   % Convert area to square feet
                   A = A * 3.28084^2 ;
                   
                   type = 'floating head';
                   
                   switch type
                       case 'floating head'
                           a = 12.031 ;
                           b = -0.8709 ;
                           c = 0.09005 ;
                       case 'fixed head'
                           a = 11.4185 ; 
                           b = -0.9228 ;
                           c = 0.09861 ;
                       case 'U-tube'
                           a = 11.551 ; 
                           b = -0.9186 ; 
                           c = 0.0979 ;
                       case 'kettle'
                           a = 12.331 ;
                           b = -0.8709 ; 
                           c = 0.09005 ;
                   end
                   
                   COST = exp( a + b * log(A) + c * (log(A))^2) ;
                   
                   P = 0.000145038 * max(obj.H(1).pin,obj.C(1).pin) ; % Convert max pressure to psi
                   if P < 100 
                       f2 = 1.0 ;
                   else
                       f2 = 0.9803 + 0.018 * (P/100) + 0.0017 * (P/100)^2 ;
                   end
                   
                   mat = 'CS' ;
                   
                   switch mat
                       case 'CS'
                           d = 0.0 ;
                           e = 0.0 ;
                       case 'SS'
                           d = 2.7 ;
                           e = 0.07 ;
                       case 'Monel'
                           d = 3.3 ;
                           e = 0.08 ;
                       case 'Ti'
                           d = 9.6 ;
                           e = 0.06 ;
                   end
                    
                   f3 = d + (A/100)^e ;
                   
                   COST = COST * f2 * f3 ;
                   COST = COST * CEind(curr) / CEind(2017) ;    
               
               case 10
                   % Plate-frame heat exchanger. Stainless and carbon steel
                   % Eq. HX3.3
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   COST = 635.14 * A ^ 0.778 ;
                   COST = COST * CEind(curr) / CEind(1986) ;
                   
               case 11
                   % Fin-plate heat exchanger. Up to 150 bar
                   % Eq. HX4.1
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   COST = 5e3 + 450 * A ^ 0.82 ;
                   COST = COST * CEind(curr) / CEind(2009) ;
                   
               case 12
                   % Coil-wound heat exchanger. Over 150 bar
                   % Eq. HX4.2
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   COST = 10e3 + 900 * A ^ 0.82 ;
                   COST = COST * CEind(curr) / CEind(2009) ;
                   
               % ** EXCHANGERS SUGGESTED FOR sCO2 CYCLES
               case 20
                   % Primary heat exchanger. Eq. HX 7
                   if (obj.UA0 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   COST = 3.5 * obj.UA0 * CEind(curr) / CEind(2017);
                   
               case 21
                   % Primary heat exchanger. Eq. HX 8
                   if (obj.UA0 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   COST = 5. * obj.UA0 * CEind(curr) / CEind(2016);
                   
               case 22
                   % Recuperator. Eq. HX1
                   if (obj.UA0 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   COST = 1.25 * obj.UA0 * CEind(curr) / CEind(2017);
               
               case 23
                   % Recuperator. Eq. HX2
                   if (obj.UA0 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   COST = 2.5 * obj.UA0 * CEind(curr) / CEind(2016);
                   
               case 24
                   % Natural gas fired primary heat exchanger
                   % From Weiland et al 2019. Valid up to 715 C and 230-275 bar, 10-50 MWth
                   COST = 6.329e5 * (obj.QS(end,1)/1e6) ^ 0.6 ;
                   if max(obj.H(1).T) > 550 + 273.15
                       COST = COST * (1 + 5.4e-5 * (max(obj.H(1).T) - 550 - 273.15)^2) ;
                   end
                   COST = COST * CEind(curr) / CEind(2017) ;
                   
               case 25
                   % Recuperator
                   % From Weiland et al 2019. UA in range 1.6e5 - 2.15e8 W/mK
                   if (obj.UA0 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   COST = 49.45 * obj.UA0 ^ 0.7544 ;
                   if max(obj.H(1).T) > 550 + 273.15
                       COST = COST * (1 + 0.02141 * (max(obj.H(1).T) - 550 - 273.15)) ;
                   end
                   COST = COST * CEind(curr) / CEind(2017) ;
                   
               % ** AIR COOLERS
               case 30
                   % Air cooler for sCO2 system. Eq. HX9
                   if (obj.UA0 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   COST = 2.3 * obj.UA0 * CEind(curr) / CEind(2017);
                   
               case 31
                   % For ideal-gas PTES system. Eq. HX12
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   COST = (5.3e5/3563) * A ^ 0.9 ;
                   COST = COST * CEind(curr) / CEind(2017) ;
                   
               case 32
                   % Air cooler for sCO2 system from Weiland et al 2019
                   % UA: 8.6e5-7.5e7, T: 50-170C, p: 54-100bar
                   if (obj.UA0 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   COST = 32.88 * obj.UA0^0.75 * CEind(curr) / CEind(2017);
                   
               case 33
                   % From Couper et al 3rd ed 2012 (2009 costs). 
                   % For 0.05 < A < 500 Ksqft
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   % Convert area to kilo-square feet
                   A = (A * 3.28084^2)/1e3 ;
                   
                   COST = 30 * A ^ 0.4 ;
                   COST = COST * 1e3 * CEind(curr) / CEind(2009) ;
                   
           end
           
           obj.hx_cost.COST = COST ;
           obj.Qact         = obj.QS(end,1) ; % Temporary
           obj.hx_cost.cost = COST / obj.Qact ;
           
           % Temporary
           if isempty(obj.hx_cost.COST)
               obj.hx_cost.COST = 0.001 ;
           end
           
       end
   end
end