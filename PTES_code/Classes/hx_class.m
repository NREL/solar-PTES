classdef hx_class
   properties
       name
       mode
       
       eff      % Effectiveness
       ploss    % Pressure loss
       Hname    % Hot fluid name
       Cname    % Cold fluid name
       
       Qmax     % Max heat transfer - design point
       Qact     % Actual heat transfer - design point
       
       Nsave    % How many load cycles to save data for. All of them, or just two (first charge and first discharge)?
       Ngrid    % Number of points along HX
       
       % **
       % The following variables are arrays - for each Nsave (which may equal number of time periods)
       HT       % Hot fluid temp
       Hh       % Hot fluid enthalpy
       Hpin     % Hot fluid inlet temp
       Hmdot    % Hot fluid mass flow
       
       % These are duplicated for hot and cold rather than have a bazillion
       % nested structures (currently have a bazillion-1)
       CT       % Hot fluid temp
       Ch       % Hot fluid enthalpy
       Cpin     % Hot fluid inlet temp
       Cmdot    % Hot fluid mass flow
       
       QS      % Cumulative heat transfer
       % **
       
       % Loss factors - one for each time period
       w       % specific work transfer, J/kg
       q       % specific heat transfer, J/kg
       Dh      % specific enthalpy change, J/kg
       sirr    % specific entropy generation, J/kgK
       
       W        % work transfer, W
       Q        % heat transfer, W
       DH       % enthalpy change, W
       Sirr     % entropy generation, W/K
       
       % Geometry
       UA       % Conductance, W/K - (design)
       A        % Heat transfer area, m2 - (design)
       dh       % Hydraulic diameter, m - (design)
       
       % Costs
       cost_mode
       COST     % Total capital cost, $
       cost     % Capital cost / Heat transfer, $/kW
       
   end
   
   methods
       function obj = hx_class(name, mode, eff, ploss, cost_mode, Nsave, Ngrid, numPeriods)
            obj.name      = name ;
            obj.mode      = mode ;
            obj.eff       = eff ;
            obj.ploss     = ploss ;
            obj.cost_mode = cost_mode ;
            
            % Loss data          
            obj.w    = zeros(numPeriods,1) ;
            obj.q    = zeros(numPeriods,1) ;
            obj.Dh   = zeros(numPeriods,1) ;
            obj.sirr = zeros(numPeriods,1) ;
                        
            obj.W    = zeros(numPeriods,1) ;
            obj.Q    = zeros(numPeriods,1) ;
            obj.DH   = zeros(numPeriods,1) ;
            obj.Sirr = zeros(numPeriods,1) ;
            
            % Property data
            obj.HT    = zeros(Nsave,Ngrid) ;
            obj.Hh    = zeros(Nsave,Ngrid) ;
            obj.Hpin  = zeros(Nsave,Ngrid) ;
            obj.Hmdot = zeros(Nsave,Ngrid) ;
            
            obj.CT    = zeros(Nsave,Ngrid) ;
            obj.Ch    = zeros(Nsave,Ngrid) ;
            obj.Cpin  = zeros(Nsave,Ngrid) ;
            obj.Cmdot = zeros(Nsave,Ngrid) ;
            
            obj.QS    = zeros(Nsave,Ngrid) ;
            
       end
        
       
       % Calculate the HX cost
       function [obj] = HX_cost(obj)
           
           switch obj.cost_mode
               case 0
                   if (obj.UA == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   obj.COST = 3500 * obj.UA ;
               case 1
                   if (obj.A == 0)
                       error('Have picked an unsuitable HX cost mode')
                   end
                   obj.COST = 9583.8 + 251.5 * obj.A ;
               case 2
                   obj.COST = 0 ;
               case 3
                   error('Not implemented')
           end
           
           obj.cost = obj.COST / obj.Qact ;
                      
       end
   end
end