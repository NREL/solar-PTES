classdef fluid_class
    properties
        name   % substance (helium, air, argon, etc.)
        job    % 'WF' (working fluid), 'SM' (storage media) or 'ENV' (environment)
        read   % 'CP' (CoolProp) or 'TAB' (table)
        handle % integer to identify CoolProp AbstractState
        HEOS   % integer to identify CoolProp AbstractState (with HEOS backend)
        TAB    % matrix of thermophysical properties
        IDL    % Structure that contains ideal gas properties - e.g. cp, cv, R, T0, P0, h0, s0
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
                    obj.HEOS   = calllib('coolprop','AbstractState_factory','HEOS',name,ierr,herr,buffer_size);
                case 'TAB'
                    
                    if computer() == "PCWIN64"
                       pe = pyenv;
                       if ~any(pe.Version==["2.7" "3.7" "3.8"])
                          error("CPython is required to run some of the code in create_table. Matlab only supports python 2.7, 3.7, 3.8"); 
                       end
                    end
                    
                    % Read thermophysical properties from table
                    obj.TAB = create_table(name); 
                case 'IDL'
                    % In this case backend is used to provide a structure of some useful properties
                    obj.IDL.T0 = backend.T0 ;
                    obj.IDL.P0 = backend.P0 ;
                    obj.IDL.cp = backend.cp ;
                    obj.IDL.cv = backend.cv ;
                                      
                    obj.IDL.mu0   = backend.mu0 ;
                    obj.IDL.TVref = backend.TVref ;
                    obj.IDL.S     = backend.S ;
                    obj.IDL.k     = backend.k ;
                    
                    % Now calculate some other reference points
                    obj.IDL.R   = obj.IDL.cp - obj.IDL.cv ;
                    obj.IDL.gam = obj.IDL.cp / obj.IDL.cv ;
                    obj.IDL.h0  = obj.IDL.cp * obj.IDL.T0 ;
                    obj.IDL.s0  = obj.IDL.cp * log(obj.IDL.T0) - obj.IDL.R * log(obj.IDL.P0) ;
                                        
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
        
        function [ Mdot ] = total_mdot( fluid, iL, ind )
            % TOTAL_MDOT Compute total mass flow rate of a number of fluid
            % streams
            %
            %       USAGE:  TOTAL_MDOT( fluid, iL, ind )
            %
            %       fluid   is an instance of the fluid_class
            %       iL      is the index of the Load period
            %       ind     is an array
            
            Mdot  = 0; % Mass flow rate into sink tank
            for i = ind
                Mdot  = Mdot  + fluid.state(iL,i).mdot;
            end
            
        end
        
%         % Calculate the viscosity of an ideal gas using Sutherland's law
%         function mu = fvisc(fluid,T)
%             
%             mu0   = fluid.IDL.mu0 ;
%             TVref = fluid.IDL.TVref ;
%             S     = fluid.IDL.S ;
%             
%             mu    = mu0 .* ((TVref + S) ./ (T + S)) .* (T./TVref).^1.5 ;
%             
%         end
        
    end
end