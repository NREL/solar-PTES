classdef compexp_class
    properties
        type      % "CMP" or "EXP"
        aim_mode  % "Paim" or "Taim"
        eff_mode  % "isen" or "poly"
        eta0  = 0 % Design value of isen/polytropic efficiency
        mdot0 = 0 % Design mass flow rate, kg/s
        pr0   = 0 % Design pressure ratio
                
        % Following variables are arrays (one for each time period)
        pr      % Pressure ratio
        eta     % Efficiency that is used in calculations 
        
        w    = 0 % specific work transfer, J/kg
        q    = 0 % specific heat transfer, J/kg
        Dh   = 0 % specific enthalpy change, J/kg
        sirr = 0 % specific entropy generation, J/kgK
        
        mdot  = 0 % Mass flow rate, kg/s
        
        W    = 0 % work transfer, W
        Q    = 0 % heat transfer, W
        DH   = 0 % enthalpy change, W
        Sirr = 0 % entropy generation, W/K
        
        % Economics
        COST    % Cost in USD
        cost    % 'Specific cost' =  COST / WORK
        cost_mode
        rhoP    % design point power density
        sd      % Standard deviation (for uncertainty analysis)
        
    end
    
    methods
        function obj = compexp_class(type, aim_mode, eff_mode, eta0, numPeriods)
            obj.type     = type ;
            obj.aim_mode = aim_mode ;
            obj.eff_mode = eff_mode ;
            obj.eta0     = eta0 ;
            
            obj.eta = zeros(numPeriods,1) ;
            obj.pr  = zeros(numPeriods,1) ;
            
            obj.w    = zeros(numPeriods,1) ;
            obj.q    = zeros(numPeriods,1) ;
            obj.Dh   = zeros(numPeriods,1) ;
            obj.sirr = zeros(numPeriods,1) ;
            
            obj.mdot  = zeros(numPeriods,1) ;
            
            obj.W    = zeros(numPeriods,1) ;
            obj.Q    = zeros(numPeriods,1) ;
            obj.DH   = zeros(numPeriods,1) ;
            obj.Sirr = zeros(numPeriods,1) ;
            
        end
            
        function [obj,fluid,i] = compexp_func (obj,fluid,ind,aim)
            
            keyboard

            % Import fluid.state and fluid.stage
            state = fluid.state(ind(1),ind(2));
            stage = fluid.stage(ind(1),ind(2)); % delete this eventually
            
            % Extract initial conditions
            T1   = state.T;
            p1   = state.p;
            h1   = state.h;
            s1   = state.s;
            rho1 = state.rho;
            obj.mdot(ind(1)) = state.mdot ;
            
            % Extract eta for this time period
            obj = compexp_offdesign (obj , ind(1) , 1) ;
            etaI = obj.eta(ind(1)) ;
            
            
            % EXPECT ALL OF THIS COULD BE MOVED INTO CLASS CONSTRUCTOR AND
            % THEREFORE ONLY CALCULATED ONCE >>>>
            if strcmp(obj.type,'CMP') && strcmp(obj.aim_mode,'Taim') && strcmp(obj.eff_mode,'poly') 
                mode = 0 ;
            elseif strcmp(obj.type,'EXP') && strcmp(obj.aim_mode,'Paim') && strcmp(obj.eff_mode,'poly') 
                mode = 1 ;
            elseif strcmp(obj.type,'EXP') && strcmp(obj.aim_mode,'Taim') && strcmp(obj.eff_mode,'poly')  
                mode = 2 ;
            elseif strcmp(obj.type,'CMP') && strcmp(obj.aim_mode,'Paim') && strcmp(obj.eff_mode,'poly') 
                mode = 3;
            elseif strcmp(obj.type,'EXP') && strcmp(obj.aim_mode,'Paim') && strcmp(obj.eff_mode,'isen')
                mode = 4 ;
            elseif strcmp(obj.type,'CMP') && strcmp(obj.aim_mode,'Paim') && strcmp(obj.eff_mode,'isen')
                mode = 5 ;
            else
                error('Not implemented')
            end
            % <<<< MOVE TO CONSTRUCTOR?
            
            % DELETE THIS EVENTUALLY? >>>
            % Determine whether it is a compression or an expansion process. Legend:
            % mode = 0: compressor. final temperature is specified
            % mode = 1: expander.   final pressure is specified
            % mode = 2: expander.   final temperature is specified
            % mode = 3: compressor. final pressure is specified
            % mode = 4: expander.   final pressure is specified (isentropic efficiency)
            % mode = 5: compressor. final pressure is specified (isentropic efficiency)
            if (mode==0 || mode==3 || mode==5) %compressor
                n = 1;
                stage.type = 'comp';
            elseif (mode==1 || mode==2 || mode==4) %expander
                n = -1;
                stage.type = 'exp';
            else
                error('***Mode not implemented***')
            end
            % <<<< DELETE THIS EVENTUALLY?
            
            % Set number of sections
            num = 100;
            
            % Compressor/expander. Final temperature is specified (aim = T_final)
            if (mode==0 || mode==2)
                % Estimate polytropic index to estimate final pressure
                T2   = aim;
                TAV  = 0.5*(T1 + T2);
                Gama = CP1('PT_INPUTS',p1,TAV,'CPMASS',fluid.handle)/CP1('PT_INPUTS',p1,TAV,'CVMASS',fluid.handle);
                phi  = (Gama/(Gama-1))*eff^n;
                p2   = p1*(T2/T1)^phi;
                
                % Compute compression/expansion for estimated final pressure
                h2   = nested_compexp(fluid,p1,h1,s1,rho1,p2,etaI,n,num);
                Tnew = CP1('HmassP_INPUTS',h2,p2,'T',fluid.handle);
                
                % Re-compute polytropic index and adapt final pressure
                phi  = log(p2/p1)/log(Tnew/T1);
                p2   = p1*(T2/T1)^phi;
                
                % Re-compute compression/expansion process
                h2   = nested_compexp(fluid,p1,h1,s1,rho1,p2,etaI,n,num);
            end
            
            % Compressor/expander. Final pressure is specified (aim = p_final)
            if (mode==3 || mode==1)
                p2   = aim;
                h2   = nested_compexp(fluid,p1,h1,s1,rho1,p2,etaI,n,num);
            end
            
            % Compressor/expander. Final pressure is specified (isentropic efficiency)
            if (mode==4 || mode==5)
                p2   = aim;
                h2   = nested_compexp_is(fluid,h1,s1,p2,etaI,n);
            end
            
            % Update state
            state.h = h2;
            state.p = p2;
            state   = update_state(state,fluid.handle,fluid.read,fluid.TAB,2);
            
            obj.Dh   = state.h - h1;
            obj.q    = 0;
            obj.w    = -stage.Dh;
            obj.sirr = state.s - s1;
            
            % DELETE THIS EVENTUALLY >>
            % Compute energy flows along stage
            stage.Dh   = state.h - h1;
            stage.q    = 0;
            stage.w    = -stage.Dh;
            stage.sirr = state.s - s1;
            
            % Export computed state and stage back into fluid
            fluid.state(ind(1),ind(2)+1) = state; % Result goes into next state
            fluid.stage(ind(1),ind(2))   = stage; % Result stays in current stage
            
            %Increase stage counter
            i = ind(2)+1;
            
            % << DELETE THIS EVENTUALLY 
            
            function h2 = nested_compexp(fluid,p1,h1,s1,rho1,p2,eta,n,num)
                % Use isentropic efficiency as an initial guess to speed things up
                % Pressure array (to solve integral)
                pv = logspace(log10(p1),log10(p2),num);
                Dp = pv(2:num) - pv(1:(num-1));
                % Initial guess
                h2_is = CP1('PSmass_INPUTS',p2,s1,'H',fluid.handle);
                h2 = h1 + eta^n*(h2_is - h1);
                % Update until convergence
                err = zeros(1,20);
                for i1 = 1:20
                    h2_0  = h2;
                    rho2  = CP1('HmassP_INPUTS',h2,p2,'D',fluid.handle);
                    xi    = log(rho2/rho1)/log(p2/p1); %assumes rho = K*p^xi along polytropic compression/expansion
                    rhov  = rho1*(pv/p1).^xi;  %density array (estimate)
                    rhoAV = 0.5*(rhov(1:(num-1))+rhov(2:num));
                    dh    = Dp./(rhoAV*eta^n); %apply polytropic efficiency definition
                    h2    = h1 + sum(dh);
                    err(i1) = abs((h2_0-h2)/(h2-h1));
                    if err(i1)<1e-6
                        break
                    else
                    end
                end
                if err(i1)>1e-6
                    error('***Convergence not found***')
                end
                %     % Plot convergence
                %     figure(5)
                %     semilogy(1:length(err),err*100)
                %     xlabel('Iterations')
                %     ylabel('Error [$\%$]')
                %     grid on
            end
            
            function h2 = nested_compexp_is(fluid,h1,s1,p2,eta,n)
                % Use isentropic efficiency
                h2_is = CP1('PSmass_INPUTS',p2,s1,'H',fluid.handle);
                if n == 1 %compressor
                    h2 = h1 + (h2_is - h1)/eta;
                elseif n==-1 %expander
                    h2 = h1 - eta*(h1 - h2_is);
                else
                    error('n must be either 1 or -1')
                end
            end
            
        end
        
        % Calculate the off-design performance of the compressor given how
        % the mass flow or pressure ratio deviate from the design point
        function obj = compexp_offdesign(obj,i,mode)
            % i is the current timestep
            % Mode 1: Specify current mass flow rate
            % Mode 2:
            if obj.type == "CMP"
                switch mode
                    case 1
                        % Currently do not have an off-design map
                        obj.eta(i) = obj.eta0 * obj.mdot(i) / obj.mdot0 ;
                    case 2
                        error('Not implemented')
                end
            elseif obj.type == "EXP"
                switch mode
                    case 1
                        % Currently do not have an off-design map
                        obj.eta(i) = obj.eta0 * obj.mdot(i) / obj.mdot0 ;
                    case 2
                        error('Not implemented')
                end
            end
        end
        
        
        % Calculate energy totals for each load cycle
        function obj = compexp_energy(obj,T)
            % T is the duration of the load cycle in seconds
            obj.W    = obj.w * obj.mdot * T ;
            obj.Q    = obj.q * obj.mdot * T ;
            obj.DH   = obj.dh * obj.mdot * T ;
            obj.Sirr = obj.sirr * obj.mdot * T ;
        end
               
        % Calculate the economic cost of the compressor or expander
        function obj = compexp_econ(obj)
            mode = obj.cost_mode ;
            
            % Expect many, many correlations may be added to this
            switch mode
                case 0
                    % "1 $/W" approximation
                    obj.COST = 1 * obj.W ;
                    obj.cost = obj.COST / obj.W ;   % This looks silly now but will be informative when more detailed correlations are used
                case 1
                    % Cost for an sCO2 compressor developed by Carlson et al. 2017
                    % See equation C6 in Q2 solar-PTES report
                    obj.COST = 6998 * obj.W ^ 0.7865 ;
                    obj.cost = obj.COST / obj.W ;
                case 2
                    error('Not implemented')
            end
            
        end
        
        
    end
end