classdef hx_class
   properties
       name
       mode
       
       model
       stage_type
       
       eff      % Effectiveness
       ploss    % Pressure loss
       DT       % Pinch-point temperature difference
       
       plossH0  % Pressure loss, specific to hot side
       plossC0  % Pressure loss, specific to cold side
       UA0      % Conductance, W/K - (design)
       NTU0     % Design NTU
       LMTD0    % Design LMTD
       iL0      % First load period in which hx is called
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
       AS       % Cumulative heat transfer area
       dL       % Array of heat exchanger sections (length)
       Ul       % Local overall heat transfer coefficient
       LMTD     % Why not
       Cmin 
       DTmin
       effDT
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
       Lgeom_set % Has the geometry been set
       hx_cost = econ_class(0,0,0,0) ;
       mat_fac % Additional cost factor, depending on material used (temperature dependent).
       
   end
   
   methods
       function obj = hx_class(name, stage_type, cost_mode, Ngrid, Nsave, numPeriods, model, varargin)
           % There are four possible ways to construct a hx, depending on
           % the selected hx model.
           %
           % If model is 'eff',
           % eff = varargin{1}, ploss = varargin{2} and D1 = varargin{3}.
           %
           % If model is 'UA',
           % UA  = varargin{1}, ploss = varargin{2} and D1 = varargin{3}.
           %
           % If model is 'DT',
           % DT  = varargin{1}, ploss = varargin{2} and D1 = varargin{3}.
           %
           % If model is 'geom',
           % DT  = varargin{1}, ploss = varargin{2}, D1 = varargin{3} and
           % shape = varargin{4}.
           
           switch model
               case 'eff'
                   if length(varargin)~=4
                       error('incorrect number of inputs');
                   end
                   obj.eff   = varargin{1};
               case 'UA'
                   if length(varargin)~=4
                       error('incorrect number of inputs');
                   end
                   obj.UA    = varargin{1};
               case 'DT'
                   if length(varargin)~=4
                       error('incorrect number of inputs');
                   end
                   obj.DT    = varargin{1};
               case 'geom'
                   if length(varargin)~=4
                       error('incorrect number of inputs');
                   end
                   obj.eff   = varargin{1};
               otherwise
                   error('not implemented')
           end
           obj.ploss = varargin{2};
           obj.D1    = varargin{3};
           obj.shape = varargin{4};
           
           obj.name       = name ;
           obj.stage_type = stage_type ;
           obj.model      = model ;
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
           
           % Performance data
           obj.Cmin = zeros(numPeriods,1);
           obj.NTU  = zeros(numPeriods,1);
           obj.DppH = zeros(numPeriods,1);
           obj.DppC = zeros(numPeriods,1);
           obj.UA   = zeros(numPeriods,1);
           obj.LMTD = zeros(numPeriods,1);
           obj.DTmin= zeros(numPeriods,1);
           obj.effDT= zeros(numPeriods,1);
           
           % Property data
           obj.H(1:numPeriods) = stream;
           obj.C(1:numPeriods) = stream;
           
           obj.QS   = zeros(Nsave,Ngrid+1) ;
           obj.AS   = zeros(Nsave,Ngrid+1) ;
           obj.Ul   = zeros(Nsave,Ngrid+1) ;
           obj.dL   = zeros(Nsave,Ngrid) ;
           
           obj.Lgeom_set = false ;
           
           obj.hx_cost = econ_class(cost_mode, 0.2, 5, 0.2) ;
           
       end
       
       % Calculate the energy and entropy terms for the heat exchanger
       function obj = hx_energy(obj , T)  
           % T is the duration of the load cycle in seconds
           
           % Iterate through each load cycle
           for i = 1 : numel(obj.H)
               
               % Only evaluate if mdot is defined
               if ~isempty(obj.H(i).mdot)
                   % Hot streams
                   obj.W(i,1)    = obj.w(i,1)    * obj.H(i).mdot * T(i) ;
                   obj.Q(i,1)    = obj.q(i,1)    * obj.H(i).mdot * T(i) ;
                   obj.DH(i,1)   = obj.Dh(i,1)   * obj.H(i).mdot * T(i) ;
                   obj.Sirr(i,1) = obj.sirr(i,1) * obj.H(i).mdot * T(i) ;
                   
                   % Cold streams
                   obj.W(i,2)    = obj.w(i,2)    * obj.C(i).mdot * T(i) ;
                   obj.Q(i,2)    = obj.q(i,2)    * obj.C(i).mdot * T(i) ;
                   obj.DH(i,2)   = obj.Dh(i,2)   * obj.C(i).mdot * T(i) ;
                   obj.Sirr(i,2) = obj.sirr(i,2) * obj.C(i).mdot * T(i) ;
               end
               
           end
       end
       
       
       % Calculate the HX cost
       % Costs come from Q2 report unless otherwise stated
       function [obj] = HX_cost(obj, CEind)
           
           curr = 2019 ; % Current year
           switch obj.hx_cost.cost_mode
               case 0
                   COST = 0.01 ;
                   
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
                   nHX = 1 ; % Number of heat exchangers
                   
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   % Convert area to square feet
                   A = A * 3.28084^2 ;
                   
                   maxA = 12e3 ;
                   if A > maxA
                       nHX = A / maxA ;
                       A = maxA ;
                   end                       
                                      
                   COST = 1.218 * nHX * exp(8.821 - 0.30863 * log(A) + 0.0681 * (log(A))^2) ;
                   
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
                   
                   mat = 'CS' ;
                   
                   switch mat
                       case 'CS'
                           f3 = 1.0 ;
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
                   nHX = 1 ; % Number of heat exchangers
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   % Convert area to square feet
                   A = A * 3.28084^2 ;
                   
                   maxA = 12e3 ;
                   if A > maxA
                       nHX = A / maxA ;
                       A = maxA ;
                   end                       
                   
                   type = 'fixed head';
                   
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
                   
                   COST = COST * f2 * f3 * nHX;
                   COST = COST * CEind(curr) / CEind(2017) ;  
                   
               case 5
                   % Shell and tube. Data from Hall 1990, based on his 1986
                   % work.
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   
                   mat = 'CS' ;
                   
                   switch mat
                       case 'CS'
                           fm = 1.0 ;
                       case 'SS'
                           fm = 1.7853 ;
                       case 'Ti'
                           fm = 5.876 ;
                   end
                   
                   P = max(obj.H(1).pin,obj.C(1).pin) ; % Max pressure
                   if P <= 10e5 
                       fp = 1.0 ;
                   else
                       fp = 1. + 0.018347 * (P - 10e5) / 1e5 ;
                   end
                   
                   COST = 30800 + 750 * fm * fp * A ^ 0.81 ; 
                   COST = COST * CEind(curr) / CEind(1986) ;  
               
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
                   COST = 6.329e5 * (obj.QS(1,end)/1e6) ^ 0.6 ;
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
                   nHX  = 1 ;
                   maxA = 500 ;
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   % Convert area to kilo-square feet
                   A = (A * 3.28084^2)/1e3 ;
                   if A > maxA
                       nHX = A / maxA ;
                       A   = maxA ;
                   end
                   
                   COST = nHX * 30 * A ^ 0.4 ;
                   COST = COST * 1e3 * CEind(curr) / CEind(2009) ;
                   
               case 34
                   % From NETL 1998. Gives very large results.
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   
                   COST = 21126 + 210.28 * A ;
                   COST = COST * CEind(curr) / CEind(1998) ;
                   
           end
           
           % Find the maximum temperature in the heat exchanger
           len  = length(obj.H) ;
           maxT = 0 ;
           for i = 1 : len
               temp = max(obj.H(i).T) ;
               if ~isempty(temp)
                   maxT = max(maxT,temp) ;
               end
           end
           
           % Multiply the cost by a factor according to the temperature
           if maxT <= 400 + 273.15
               obj.mat_fac = 1.0 ; % Carbon steel
           elseif maxT > 400 + 273.15 && maxT <= 600 + 273.15
               obj.mat_fac = 2.0 ; % Stainless steel
           elseif maxT > 600 + 273.15
               obj.mat_fac = 5.0 ; % Some kind of nickel alloy
           end
           
           switch obj.hx_cost.cost_mode
               case {24,25} % These correlations already seem to take material considerations into account
                   obj.mat_fac = 1.0 ;
           end
                      
           COST = COST * obj.mat_fac ;
           
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