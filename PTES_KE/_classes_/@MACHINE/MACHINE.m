classdef MACHINE 
    % MACHINE is a construct for compressors and expanders
    
    properties
        name    % name
        note    % Any other details
        in      % FLOW class - inlet flow properties
        out     % FLOW class - outlet flow properties
        eta     % Isentropic efficiency
        alp     % Heat loss factor - currently unused
        beta    % Pressure ratio
        tau     % Temperature ratio
        w       % Specific work out of the system (positive), kJ/kg
        q       % Specific heat transfer into the system (positive), kJ/kg
        W       % Work out of the system (positive), MW
        Q       % Heat transfer into the system (positive), kJ
        WLTH    % Work lost due to thermal irreversibilities, kJ/kg
        WLQL    % Work lost due to heat leakage, kJ/kg
        WLT     % Work lost total, kJ/kg   
        type    % 'CMP' for compressor, 'EXP' for expander
    end
    
    % Constructor function
    methods
        function obj = MACHINE(machine_dat) 
           if nargin == 1 
                                            
               if isfield(machine_dat , 'type')
                   obj.type = machine_dat.type ;
               else
                   error('MACHINE class is not provided with the machine type (CMP or EXP)');
               end
               
               if isfield(machine_dat , 'beta')
                   obj.beta = machine_dat.beta ;
               else
                   error('MACHINE class is not provided with the pressure ratio');
               end
               
               if isfield(machine_dat , 'eta')
                   obj.eta = machine_dat.eta ;
               else
                   obj.eta = 1.0 ;
               end
               
               % Set up flow classes
               flow_dat.nobv = 1 ;
               obj.in  = FLOW(flow_dat) ;
               obj.out = FLOW(flow_dat) ;
               
               if isfield(machine_dat , 'Tin')
                   obj.in.T = machine_dat.Tin ;
               end
               
               % Inlet pressure may not be specified
               if isfield(machine_dat , 'Pin')
                   obj.in.P = machine_dat.Pin ;
               end
                              
           end
        end
    end
    
    % Other functions
    methods
        % This function calculates the compression/expansion process -
        % assumes all inlet flow properties are known
        function obj = calc_machine(obj,fld)
            
            % Calculate the outlet pressure
            if (obj.type == 'CMP')
                obj.out.P = obj.beta * obj.in.P ; 
                n         = 1.0 ;
            else
                obj.out.P = obj.in.P / obj.beta ;
                n         = -1.0 ;
            end
            
            % Calculate isentropic behaviour
            obj.out.s = obj.in.s ;
           
            
            % Calculate remaining flow properties at the outlet

            obj.out = flow_state(obj.out,'Ps',fld) ;
            dH      = obj.in.h - obj.out.h ;

            % Calcaulte the correct outlet enthalpy
            obj.out.h = obj.in.h + (obj.out.h - obj.in.h) / (obj.eta ^ n) ;
            
            % Re-calculate remaining flow properties at the outlet
            obj.out = flow_state(obj.out,'Ph',fld) ;
            obj.tau = (obj.out.T / obj.in.T )^n ;
            
            % Calculate work terms    
            obj.w    = obj.in.h - obj.out.h ;
            obj.WLTH = abs(dH - obj.w) ;
            
        end
        
    end
        
    
    
end