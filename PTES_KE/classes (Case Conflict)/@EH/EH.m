classdef EH
    % EH is a construct for electric heaters
    
    properties
        name    % name
        note    % Any other details
        mdot    % Mass flow rate through the electric heater
        in      % FLOW class - inlet flow properties
        out     % FLOW class - outlet flow properties
        ne      % Heater efficiency
        dP      % Pressure loss factor
        tau     % Temperature ratio
        q       % Specific heat transfer into the system (positive), kJ/kg
        qh      % Electric heater heat
        cost    % Cost of Electric Heater based on Benato 2017 equation ($)
        type    % 'Joule' for all electric heaters
    end
    
    % Constructor function
    methods
        function obj = EH(EH_dat) 
           if nargin == 1 
                                            
               if isfield(EH_dat , 'type')
                   obj.type = EH_dat.type ;
               else
                   error('MACHINE class is not provided with the machine type (CMP or EXP)');
               end
               
               if isfield(EH_dat , 'ne')
                   obj.ne = EH_dat.ne ;
               else
                   obj.ne = 1.0 ;
               end
               
               if isfield(EH_dat , 'mdot')
                   obj.mdot = EH_dat.mdot ;
               else
                   obj.mdot = 1 ;
               end
              
               % Set up flow classes
               flow_dat.nobv = 1 ;
               obj.in  = FLOW(flow_dat) ;
               obj.out = FLOW(flow_dat) ;

               % Inlet pressure and temp may not be specified
               if isfield(EH_dat , 'Pin')
                   obj.in.P = EH_dat.Pin ;
                   if obj.in.P > 0
                   else
                       obj.in.P = 101325 ;
                   end
               end
               
               if isfield(EH_dat , 'Tin')
                   obj.in.T = EH_dat.Tin ;
                   if obj.in.T > 0
                   else
                       obj.in.T = 25 + 273.15 ;
                   end
               end
               
           end
        end
    end
    
    % Other functions
    methods
        % This function calculates the electric heating process -
        % assumes all inlet flow properties are known
        function obj = calc_EH(obj,fld)
            
            % Mass Flow rate
            obj.in.mdot = obj.mdot ;
            obj.out.mdot = obj.mdot ;
            
            % Calculate the outlet pressure, outlet temp is an input
            obj.out.P = obj.in.P * (1 - obj.dP) ;

            % Calculate remaining flow properties at the outlet
            obj.out = flow_state(obj.out,'PT',fld) ;
            dH      = obj.out.h - obj.in.h ;

            % Calculate the useful heat (q) and the overall heat (qh)
            obj.q  = dH ;
            obj.qh = dH/obj.ne ;

        end
        
    end
        
    
    
end