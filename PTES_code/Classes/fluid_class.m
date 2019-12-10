classdef fluid_class
    properties
        name   % substance (helium, air, argon, etc.)
        job    % 'WF' (working fluid) or 'SM' (storage media)
        read   % 'CP' (CoolProp) or 'TAB' (table)
        handle % integer to identify CoolProp AbstractState
        TAB    % matrix of thermophysical properties
        state = state_class;
        stage = stage_class;
        Nstg   % number of stages
        cost   % Cost per kg
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
            obj.Nstg(1:numPeriods,1) = 0;
            obj.cost = 0.0; % Should specify this properly in the class constructor call
        end
        
        function obj = reset_fluid(obj)
            obj.state(:,:) = state_class;
            obj.stage(:,:) = stage_class;
            obj.Nstg(:)    = 0;
        end
        
        function obj = count_Nstg(obj)
            for iL=1:numel(obj.Nstg)
                a = find(~strcmp({obj.stage(iL,:).type},'0'),1,'last')-1;
                if isscalar(a)
                    obj.Nstg(iL) = a;
                else
                    obj.Nstg(iL) = 0;
                end
            end
        end
        
        function print_states(obj,iL,v0,varargin)
            if nargin == 4 % fourth output must be the Load structure
                Load = varargin{1};
                fprintf(1,'%11s ','T [K]','p [bar]','h [MJ/kg]','s [kJ/kg/K]','mdot [kg/s]','Q [-]','Inlet of','Position','Cycle'); fprintf(1,'\n');
                for i0=v0
                    fprintf(1,'%11.1f %11.3f %11.3f %11.2f %11.2f %11.2f %11s %11d %11s\n',...
                        obj.state(iL,i0).T, obj.state(iL,i0).p/1e5, obj.state(iL,i0).h/1e6,...
                        obj.state(iL,i0).s/1e3, obj.state(iL,i0).mdot, obj.state(iL,i0).Q,...
                        obj.stage(iL,i0).type, i0, Load.type(iL));
                end
                fprintf(1,'\n');
            else
                fprintf(1,'%11s ','T [K]','p [bar]','h [MJ/kg]','s [kJ/kg/K]','mdot [kg/s]','Q [-]','Inlet of','Position','Cycle'); fprintf(1,'\n');
                for i0=v0
                    fprintf(1,'%11.1f %11.3f %11.3f %11.2f %11.2f %11.2f %11s %11d %11d\n',...
                        obj.state(iL,i0).T, obj.state(iL,i0).p/1e5, obj.state(iL,i0).h/1e6,...
                        obj.state(iL,i0).s/1e3, obj.state(iL,i0).mdot, obj.state(iL,i0).Q,...
                        obj.stage(iL,i0).type, i0, iL);
                end
                fprintf(1,'\n');
            end
        end
        
    end
end