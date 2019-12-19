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
       function [obj] = HX_cost(obj, CEind)
           
           curr = 2019 ; % Current year
           switch obj.hx_cost.cost_mode
               case 0
                   if (obj.UA0 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   COST = 3.500 * obj.UA0 * CEind(curr) / CEind(2017);
               case 1
                   if (obj.A1 == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   A = 0.5 * (obj.A1 + obj.A2) ;
                   COST = 9583.8 + 251.5 * A ;
                   COST = COST * CEind(curr) / CEind(2009) ;
               case 2
                   COST = 0 ;
               case 3
                   error('Not implemented')
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