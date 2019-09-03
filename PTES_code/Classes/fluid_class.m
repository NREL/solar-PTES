classdef fluid_class
    properties
        name   % substance (helium, air, argon, etc.)
        job    % 'WF' (working fluid) or 'SM' (storage media)
        read   % 'CP' (CoolProp) or 'TAB' (table)
        handle % integer to identify CoolProp AbstractState
        TAB    % matrix of thermophysical properties
        state = state_class;
        stage = stage_class;
    end
    methods
        function obj = fluid_class(name,job,read,backend,numPeriods,numStates)
            %GAS_CLASS Construct an instance of this class
            obj.name = name;
            obj.job  = job;
            obj.read = read;
            switch read
                case 'CP'
                    % Real properties will be defined by Coolprop at each stage
                    % Declaring error variables, selecting backend, and obtaining
                    % CoolProp handle
                    ierr = 0; buffer_size = 10;
                    herr= char((1:1:buffer_size));
                    obj.handle = calllib('coolprop','AbstractState_factory',backend,name,ierr,herr,buffer_size);
                case 'TAB'
                    % Read thermophysical properties from table
                    obj.TAB = create_table(name);                    
                otherwise
                    error('not implemented')
            end
            obj.state(1:numPeriods,1:numStates)   = state_class;
            obj.stage(1:numPeriods,1:numStates-1) = stage_class;
        end
        function obj = reset_fluid(obj)
            obj.state(:,:) = state_class;
            obj.stage(:,:) = stage_class;
        end
    end
end
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
