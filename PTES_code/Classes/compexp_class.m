classdef compexp_class
    properties
        type      % "comp" or "exp"
        aim_mode  % "Paim" or "Taim"
        eff_mode  % "isen" or "poly"
        
        % Design values
        eta0  = 0 % Design value of isen/polytropic efficiency
        mdot0 = 0 % Design mass flow rate, kg/s
        W0    = 0 % Design point work, W
        pr0   = 0 % Design pressure ratio
        Tin   = 0 % Design inlet temperature, K
        Pin   = 0 % Design inlet pressure, Pa
        N0    = 0 % Design rotational speed, rad/s
                
        % Following variables are arrays (one for each time period)
        pr      % Pressure ratio
        eta     % Efficiency that is used in calculations 
        mdot    % Mass flow rate
        N       % Speed
        
        w       % specific work transfer, J/kg
        q       % specific heat transfer, J/kg
        Dh      % specific enthalpy change, J/kg
        sirr    % specific entropy generation, J/kgK
        
        W        % work transfer, W
        Q        % heat transfer, W
        DH       % enthalpy change, W
        Sirr     % entropy generation, W/K
        
        rhoP    % design point power density
                
        % Economics
        cmpexp_cost = econ_class(0,0,0,0) ;
        
    end
    
    methods
        function obj = compexp_class(type, eff_mode, cost_mode, eta0, numPeriods)
            obj.type      = type ;
            obj.eff_mode  = eff_mode ;
            obj.eta0      = eta0 ;
            
            obj.N0  = 2 * pi * 3600 ; % Design rotational speed
                        
            obj.eta  = zeros(numPeriods,1) ;
            obj.pr   = zeros(numPeriods,1) ;
            obj.mdot = zeros(numPeriods,1) ;
            obj.N    = zeros(numPeriods,1) ;
            
            obj.w    = zeros(numPeriods,1) ;
            obj.q    = zeros(numPeriods,1) ;
            obj.Dh   = zeros(numPeriods,1) ;
            obj.sirr = zeros(numPeriods,1) ;
                        
            obj.W    = zeros(numPeriods,1) ;
            obj.Q    = zeros(numPeriods,1) ;
            obj.DH   = zeros(numPeriods,1) ;
            obj.Sirr = zeros(numPeriods,1) ;
            
            obj.cmpexp_cost = econ_class(cost_mode, 0.2, 5, 0.2) ;
            
        end
            
        function [obj,fluid,i] = compexp_func (obj,iL,fluid,i,aim_mode,aim, design_mode)
            % iL is the Load index
            % i is the index of the state point within the cycle
            
            obj.aim_mode = aim_mode ;
            % Import fluid.state and fluid.stage
            state = fluid.state(iL,i);
            stage = fluid.stage(iL,i); % delete this eventually
            
            % Extract initial conditions
            T1   = state.T;
            p1   = state.p;
            h1   = state.h;
            s1   = state.s;
            rho1 = state.rho;
            
            % Extract eta, pr for this time period
            if design_mode
                obj.eta(iL)  = obj.eta0 ;
                obj.mdot(iL) = obj.mdot0 ;
                obj.pr(iL)   = 1 ;
            else
                if strcmp(fluid.name,'Water')
                    obj = compexp_offdesign (obj , state, iL , aim, 2) ;
                elseif strcmp(obj.type,'pump')
                    obj.eta(iL)  = obj.eta0 ;
                    obj.mdot(iL) = state.mdot ;
                    obj.pr(iL)   = 1 ;
                else
                    obj = compexp_offdesign (obj , state, iL , aim, 1) ;
                    %obj.eta(iL)  = obj.eta0 ;
                    %obj.mdot(iL) = state.mdot ;
                    %obj.pr(iL)   = 1 ;
                end
            end
            etaI = obj.eta(iL) ;     
            
            % If not at the design pressure ratio then calculate the
            % required pressure ratio
            %if obj.pr(iL) ~= 1
            if ~design_mode
               obj.aim_mode = 'Paim' ;
               if strcmp(obj.type,'comp')
                   aim = fluid.state(iL,i).p * obj.pr0 * obj.pr(iL) ;
               elseif strcmp(obj.type,'exp')
                   aim = fluid.state(iL,i).p / (obj.pr0 * obj.pr(iL)) ;
               end
            end
            
            
            % EXPECT ALL OF THIS COULD BE MOVED INTO CLASS CONSTRUCTOR AND
            % THEREFORE ONLY CALCULATED ONCE >>>>
            if strcmp(obj.type,'comp') && strcmp(obj.aim_mode,'Taim') && strcmp(obj.eff_mode,'poly') 
                mode = 0 ;
            elseif strcmp(obj.type,'exp') && strcmp(obj.aim_mode,'Paim') && strcmp(obj.eff_mode,'poly') 
                mode = 1 ;
            elseif strcmp(obj.type,'exp') && strcmp(obj.aim_mode,'Taim') && strcmp(obj.eff_mode,'poly')  
                mode = 2 ;
            elseif strcmp(obj.type,'comp') && strcmp(obj.aim_mode,'Paim') && strcmp(obj.eff_mode,'poly') 
                mode = 3;
            elseif strcmp(obj.type,'exp') && strcmp(obj.aim_mode,'Paim') && strcmp(obj.eff_mode,'isen')
                mode = 4 ;
            elseif strcmp(obj.type,'comp') && strcmp(obj.aim_mode,'Paim') && strcmp(obj.eff_mode,'isen')
                mode = 5 ;
            elseif strcmp(obj.type,'pump') && strcmp(obj.aim_mode,'Paim') && strcmp(obj.eff_mode,'isen')
                mode = 6 ;
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
            elseif (mode==6) %pump
                n = 1;
                stage.type = 'pump';
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
                Gama = RP1('PT_INPUTS',p1,TAV,'CPMASS',fluid)/RP1('PT_INPUTS',p1,TAV,'CVMASS',fluid);
                phi  = (Gama/(Gama-1))*etaI^n;
                p2   = p1*(T2/T1)^phi;
                
                % Compute compression/expansion for estimated final pressure
                h2   = nested_compexp(fluid,p1,h1,s1,rho1,p2,etaI,n,num);
                Tnew = RP1('HmassP_INPUTS',h2,p2,'T',fluid);
                
                % Re-compute polytropic index and adapt final pressure
                phi  = log(p2/p1)/log(Tnew/T1);
                p2   = p1*(T2/T1)^phi;
                
                % Re-compute compression/expansion process
                h2   = nested_compexp(fluid,p1,h1,s1,rho1,p2,etaI,n,num);
            end
            
            % Compressor/expander. Final pressure is specified (aim = p_final)
            if (mode==3 || mode==1)
                p2   = aim ;
                h2   = nested_compexp(fluid,p1,h1,s1,rho1,p2,etaI,n,num);
            end
            
            % Compressor/expander. Final pressure is specified (isentropic efficiency)
            if (mode==4 || mode==5)
                p2   = aim ;
                h2   = nested_compexp_is(fluid,h1,s1,p2,etaI,n);
            end
            
            % Pump. Final pressure specified (isentropic efficiency)
            if mode==6
               p2 = aim ;
               h2 = pump_is(obj.mdot(iL), h1, rho1, p1, p2, etaI) ;
            end
            
            % Update state
            state.h = h2;
            state.p = p2;
            state   = update_state(state,fluid,2);
            
            obj.Dh(iL)   = state.h - h1;
            obj.q(iL)    = 0;
            obj.w(iL)    = -obj.Dh(iL);
            obj.sirr(iL) = state.s - s1;
            
            % DELETE THIS EVENTUALLY >>
            % Compute energy flows along stage
            stage.Dh   = state.h - h1;
            stage.q    = 0;
            stage.w    = -stage.Dh;
            stage.sirr = state.s - s1;
            
            % Export computed state and stage back into fluid
            fluid.state(iL,i+1) = state; % Result goes into next state
            fluid.stage(iL,i)   = stage; % Result stays in current stage
            % << DELETE THIS EVENTUALLY 
            %Increase stage counter
            i = i+1;
            
            % Save the main design parameters
            if design_mode
                obj.mdot0 = state.mdot ;
                obj.Tin   = T1 ;
                obj.Pin   = p1 ;
                
                % Design pressure ratio
                if strcmp(obj.type,'comp') || strcmp(obj.type,'pump')
                    obj.pr0 = p2 / p1 ;
                else
                    obj.pr0 = p1 / p2 ;
                end
                
                obj.W0 = abs(obj.w(iL) * obj.mdot0) ;
            end
                        
            function h2 = nested_compexp(fluid,p1,h1,s1,rho1,p2,eta,n,num)
                % Use isentropic efficiency as an initial guess to speed things up
                % Pressure array (to solve integral)
                pv = logspace(log10(p1),log10(p2),num);
                Dp = pv(2:num) - pv(1:(num-1));
                % Initial guess
                h2_is = RP1('PSmass_INPUTS',p2,s1,'H',fluid);
                h2 = h1 + eta^n*(h2_is - h1);
                % Update until convergence
                err = zeros(1,20);
                for i1 = 1:50
                    h2_0  = h2;
                    rho2  = RP1('HmassP_INPUTS',h2,p2,'D',fluid);
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
                    %error('***Convergence not found***')
                    warning('***Convergence not found***')
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
                h2_is = RP1('PSmass_INPUTS',p2,s1,'H',fluid);
                if n == 1 %compressor
                    h2 = h1 + (h2_is - h1)/eta;
                elseif n==-1 %expander
                    h2 = h1 - eta*(h1 - h2_is);
                else
                    error('n must be either 1 or -1')
                end
            end
            
            % Calculate the work input to a pump Q = dP * vol flow rate / eta
            function h2 = pump_is(mdot, h1, rho1, p1, p2, eta) 
                work = (p2 - p1) * mdot / (rho1 * eta) ;
                h2   = h1 + work ;
            end
            
        end
        
        % Calculate the off-design performance of the compressor given how
        % the mass flow or pressure ratio deviate from the design point
        function [obj, Nr, mr] = compexp_offdesign(obj,state,iL,paim,mode)
            % iL is the current timestep
            % Mode 1: Map based on Zhang and Cai 2002
            
            % Unpack
            T1   = state.T;
            P1   = state.p;
            obj.mdot(iL) = state.mdot ; 
            Nr = 0.;
            mr = 0.;
            
            if obj.N(iL) == 0
                obj.N(iL) = obj.N0 ;
            end
                        
            % Check whether at design point
            % (Running the off-design model at the design point should
            % return the same value, running this way requires fewer calcs
            if P1==obj.Pin && T1==obj.Tin && obj.mdot(iL)==obj.mdot0 && obj.N(iL) == obj.N0
                obj.eta(iL) = obj.eta0 ;
                obj.pr(iL)  = 1.;
                mr = 1. ;
                Nr = 1. ;
            else
            
                % Off-design compressor
                if obj.type == "comp"
                    switch mode
                        case 1
                            mr = (obj.mdot(iL) / obj.mdot0) * (T1/obj.Tin)^0.5 * (obj.Pin/P1) ; % Reduced mass
                            Nr = (obj.N(iL) / obj.N0) * (T1/obj.Tin)^0.5 ; % Reduced speed
                            % Constants
                            c1 = 1.8 ;
                            c2 = 1.6;%1.8 ; % These values seem to be variable
                            c3 = c1*(1-c2/Nr) + Nr*(Nr-c2)^2 ;
                            c4 = Nr / c3 ;
                            c5 = (c1 - 2*c2*Nr^2) / c3 ;
                            c6 = -(c1*c2*Nr - c2^2*Nr^3) / c3 ;
                            c7 = 0.3 ;
                            % Pressure ratio and isentropic efficiency
                            obj.pr(iL)  = c4 * mr^2 + c5 * mr + c6 ;
                            obj.eta(iL) = obj.eta0 *((1-c7*(1-Nr)^2)*(Nr/mr)*(2-Nr/mr)) ;
                            
                        case 2
                            error('Not implemented')
                    end
                elseif obj.type == "exp"
                    switch mode
                        case 1 % Axial gas machines
                            mr = (obj.mdot(iL) / obj.mdot0) * (T1/obj.Tin)^0.5 * (obj.Pin/P1) ; % Reduced mass
                            Nr = (obj.N(iL) / obj.N0) * (T1/obj.Tin)^0.5 ; % Reduced speed
                            
                            t1 = (1.4 - 0.4 * obj.N(iL) / obj.N0)^0.5 ;
                            t2 = 0.3 ;
                            
                            % PR equations in references seem to be wrong -
                            % based on Stodola's ellipse - use this!
                            %obj.pr(iL)  = (1+(obj.pr0^2 - 1)*(T1/obj.Tin)*(mr/t1)^2)^0.5;
                            mrat = (obj.mdot(iL) / obj.mdot0) * (T1/obj.Tin)^0.5 ;
                            obj.pr(iL)  = (1./((1.-(1-1./obj.pr0^2)*(mr/t1)^2)))^0.5;
                            obj.pr(iL) = obj.pr(iL) / obj.pr0 ;
                            if imag(obj.pr(iL)) > 0
                                warning('Imaginary pressure ratios achieved in off-design expander. If this occurs in the final converged solution you have a problem')
                                obj.pr(iL) = 1.;
                            end
                            obj.eta(iL) = obj.eta0 *(1-t2*(1-Nr)^2)*(Nr/mr)*(2-Nr/mr) ;
                            
                        case 2 % Axial steam machines
                            % Estimate what exit pressure should be based on Stodola's ellipse
                            pout0 = obj.Pin / obj.pr0 ;
                            pout = P1 * sqrt(1.-(1. - (pout0/obj.Pin)^2)*(obj.mdot(iL) / obj.mdot0)^2);
                            
                            obj.eta(iL) = obj.eta0 ; % Assume this by default
                            
                            % How does outlet presssure compare to paim?
                            if paim < pout
                                pout = paim ;
                                obj.pr(iL) = (P1 / paim) / obj.pr0 ; % Reduce pout and pr accordingly
                                
                                % Calculate mass from Stodola
                                mout = obj.mdot0 * sqrt ( (1. - (pout/P1)^2)/(1. - (pout0/obj.Pin)^2)) ;
                                if mout > obj.mdot0 % If an increased mass is required, specify this
                                    obj.mdot(iL) = mout ; 
                                elseif mout < obj.mdot0 %If a reduced mass is required, use the original mass but reduce the efficiency
                                    obj.eta(iL) = obj.eta0 * (1. - 0.191 + 0.409*mout/obj.mdot0 - 0.218*(mout/obj.mdot0)^2) ; % Correlation from Patnode thesis, p.68
                                end
                            end
                            
                            if pout < paim
                                pout = paim ; % paim is minimum outlet pressure
                                obj.pr(iL) = (P1 / paim) / obj.pr0 ; % Reduce pout and pr accordingly
                                
                                % Calculate mass from Stodola
                                mout = obj.mdot0 * sqrt ( (1. - (pout/P1)^2)/(1. - (pout0/obj.Pin)^2)) ;
                                if mout > obj.mdot0 % If an increased mass is required, specify this
                                    obj.mdot(iL) = mout ; 
                                elseif mout < obj.mdot0 %If a reduced mass is required, use the original mass but reduce the efficiency
                                    obj.eta(iL) = obj.eta0 * (1. - 0.191 + 0.409*mout/obj.mdot0 - 0.218*(mout/obj.mdot0)^2) ; % Correlation from Patnode thesis, p.68
                                end
                            end
                        case 3
                            error('Not implemented')
                    end
                    
                end
            end
        end
        
        
        % Calculate energy totals for each load cycle
        function obj = compexp_energy(obj,T)
            % T is the duration of the load cycle in seconds
            obj.W    = obj.w    .* obj.mdot .* T ;
            obj.Q    = obj.q    .* obj.mdot .* T ;
            obj.DH   = obj.Dh   .* obj.mdot .* T ;
            obj.Sirr = obj.sirr .* obj.mdot .* T ;
        end
               
        % Calculate the economic cost of the compressor or expander
        function obj = compexp_econ(obj, CEind, fld)
            mode = obj.cmpexp_cost.cost_mode ;
            curr = 2019 ; % Current year
            
            % Many, many correlations exist
            % If an equation number is given, then that corresponds to an
            % equation given in the solar-PTES Q2 report
            switch mode
                case 0
                    COST = 0.01 ;
                
                case 1
                    % Reciprocating compressor cost from Georgiou et al 2019 based on Turton
                    % ** NOT SURE IF power should be in kW or horsepower?
                    COST = 10 ^ (2.2897 + 1.3604 * log10(obj.W0/1e3) - 0.1027 * (log10(obj.W0/1e3))^2) ;
                    COST = COST * 3.3 * CEind(curr) / CEind(2013)  ;    
                
                case 2
                    % Reciprocating compressor cost from Georgiou et al 2019 based on Seider
                    COST = exp(7.9661 + 0.8 * log(obj.W0 * 0.00134102)) ;
                    COST = COST * CEind(curr) / CEind(2010)  ;        
                
                case 3
                    % Reciprocating compressor cost from Georgiou et al 2019 based on Couper
                    COST = 7190 * (obj.W0 * 0.00134102)^0.61 ;
                    COST = COST * CEind(curr) / CEind(2010)  ; 
                    
                    
                case 4
                    % This calculates the average cost from Turton, Sieder, and Couper
                
                    % Reciprocating compressor cost from Georgiou et al 2019 based on Seider
                    SIEDER = exp(7.9661 + 0.8 * log(obj.W0 * 0.00134102)) ;
                    SIEDER = SIEDER * 2.15 * CEind(curr) / CEind(2010)  ; 
                    
                    % Turton
                    TURT = 10 ^ (2.2897 + 1.3604 * log10(obj.W0/1e3) - 0.1027 * (log10(obj.W0/1e3))^2) ;
                    TURT = TURT * 3.3 * CEind(curr) / CEind(2001)  ;
                    
                    % Couper
                    COUPER = 7190 * (obj.W0 * 0.00134102)^0.61 ;
                    COUPER = COUPER * 1.3 * CEind(curr) / CEind(2010)  ; 
                    
                    COST = (SIEDER + TURT + COUPER) / 3.;
                
                    
                case 10
                    % Black and Veatch correlation, eq. C1 and E1
                    COST = 250. * obj.W0 / 1e3 ;
                    COST = COST * CEind(curr) / CEind(2018) ;
                    
                case 11
                    % Turbomachinery compressor cost from Valero et al 1994
                    % Have changed eta_max from 0.9 to 0.92
                    COST = 39.5 * obj.mdot0 * obj.pr0* log(obj.pr0) / (0.92 - obj.eta0) ;
                    COST = COST * CEind(curr) / CEind(1994) ;
                    
                case 12
                    % Turbomachinery compressor from Agazzini 1996. Based
                    % on Valero, but coefficients updated against costs
                    % from Gas Turbine World 1995 Handbook. Again have
                    % changed eta_max from 0.9 to 0.92.
                    
                    COST = 1.051 * 39.5 * obj.mdot0 * obj.pr0* log(obj.pr0) / (0.92 - obj.eta0) ;
                    COST = COST * CEind(curr) / CEind(1995) ;
                    
                case 13
                    % Air compressor costs based on quotations for CAES.
                    % From CEC-500-2008-069 CAES scoping report.
                    
                    COST = 30 * obj.mdot0 * 3600. ;
                    COST = COST * CEind(curr) / CEind(2008) ;
                    
                case 14
                    % Centrifugual compressor from Couper
                    % Convert power to horsepower
                    P = obj.W0 * 0.00134102 ;
                    
                    nC = 1 ;
                    maxP = 30e3 ;
                    if P > maxP
                        nC = P / maxP ;
                        P  = maxP ;
                    end
                    
                    COST = nC * 7.9e3 * P ^ 0.62 ;
                    COST = COST * CEind(curr) / CEind(2010) ;
                                        
                    
                case 15 
                    % Use Agazzini's (Valero's) correlation, but scale it by the density
                    % of the inlet air (if higher density, then more
                    % compact compressor)
                    
                    COST = 1.051 * 39.5 * obj.mdot0 * obj.pr0* log(obj.pr0) / (0.92 - obj.eta0) ;
                    COST = COST * CEind(curr) / CEind(1995) ;
                    
                    RHOin = RP1('PT_INPUTS',obj.Pin,obj.Tin,'D',fld) ;
                    RHO0  = 1.225 ; % Density of air at standard conditions
                    n     = 1 ; % This allows you to modify how strongly the cost is reduced
                    scale = (RHOin / RHO0)^n ;
                    
                    COST = COST / scale ;  
                    
                case 16
                    % Scale Agazzini's correlation by the maximum compressor temperature
                    % This way can capture increasing costs
                    % A typical air compressor at a pressure ratio of 17
                    % compresses air to around 450-500 C (750 K)
                    % Therefore the cost should increase above this. Assume
                    % a compressor at 1000 K is double the cost. One
                    % sensible looking factor is Tfac = (1. + exp(0.01 .* Tmax - 10)) ;
                    
                    % MAYBE COMBINE THIS CORRELATION AND THE PREVIOUS ONE?
                    
                    COST = 1.051 * 39.5 * obj.mdot0 * obj.pr0* log(obj.pr0) / (0.92 - obj.eta0) ;
                    COST = COST * CEind(curr) / CEind(1995) ;
                    
                    % Estimate gamma
                    if strcmp(fld.read,'IDL')
                        Gama = fld.IDL.gam ;
                    else
                        Gama = CP1('PT_INPUTS',obj.Pin,obj.Tin,'CPMASS',fld.handle)/CP1('PT_INPUTS',obj.Pin,obj.Tin,'CVMASS',fld.handle);
                    end
                                        
                    Tmax = obj.Tin * obj.pr0 ^ ((Gama-1.)/(obj.eta0*Gama)) ;
                    Tfac = (1. + exp(0.01 .* Tmax - 10)) ;
                    
                    COST = COST * Tfac ; 
                case 20
                    % Cost for an sCO2 compressor developed by Carlson et al. 2017
                    % See equation C6 in Q2 solar-PTES report
                    COST = 6998 * (obj.W0/1e3) ^ 0.7865 ;
                    COST = COST * CEind(curr) / CEind(2017) ;
                    
                case 21 
                    % sCO2 compressor from Morandin 2013 - eq. C7
                    COST = 20e3 + 9e3 * (obj.W0/1e3)^0.6 ;
                    COST = COST * CEind(curr) / CEind(2013) ;
                    
                case 22
                    % sCO2 compressor by Weiland et al 2019
                    % Integrally geared, centrifugal, 1.5-200 MW
                    COST = 1.23e6 * (obj.W0/1e6)^0.3992 ;
                    COST = COST * CEind(curr) / CEind(2019) ;
                    
                case 30
                    % Reciprocating Expander cost from Georgiou et al 2019 based on Turton
                    % ** NOT SURE IF power should be in kW or horsepower?
                    COST = 10 ^ (2.7051 + 1.4398 * log10(obj.W0/1e3) - 0.1776 * (log10(obj.W0/1e3))^2) ;
                    COST = COST * 3.5 * CEind(curr) / CEind(2013) ;
                
                case 31
                    % Reciprocating Expander cost from Georgiou et al 2019 based on Sieder
                    COST = 530 * (obj.W0 * 0.00134102)^0.81 ;
                    COST = COST * CEind(curr) / CEind(2010)  ;   
                    
                case 32 
                    % Reciprocating expander cost from Georgiou et al 2019 based on Couper
                    COST = 378 * (obj.W0 * 0.00134102)^0.81 ;
                    COST = COST * CEind(curr) / CEind(2010)  ; 
                    
                case 33
                    % Average reciprocating expander cost from Turton, Sieder and Couper
                    
                    SIEDER = 530 * (obj.W0 * 0.00134102)^0.81 ;
                    SIEDER = SIEDER * 3.5 * CEind(curr) / CEind(2010)  ; 
                    
                    TURT = 10 ^ (2.7051 + 1.4398 * log10(obj.W0/1e3) - 0.1776 * (log10(obj.W0/1e3))^2) ;
                    TURT = TURT * 3.5 * CEind(curr) / CEind(2001) ;
                    
                    COUPER = 378 * (obj.W0 * 0.00134102)^0.81 ;
                    COUPER = COUPER * 1.5 * CEind(curr) / CEind(2010)  ; 
                    
                    COST = (TURT + SIEDER + COUPER) / 3. ;
                
                    
                case 40
                    % Expander cost from Valero et al 1994
                    % Eq. E3.1                    
                    COST = 266.3 * obj.mdot0 * log(obj.pr0) / (0.92 - obj.eta0) ;
                    COST = COST * (1. + exp(0.018 * obj.Tin - 26.4)) ;
                    COST = COST * CEind(curr) / CEind(1994) ;
                    
                case 41
                    % Turbomachinery expander from Agazzini 1996. Based
                    % on Valero, but coefficients updated against costs
                    % from Gas Turbine World 1995 Handbook. 
                    
                    COST = 1.051 * 266.3 * obj.mdot0 * log(obj.pr0) / (0.92 - obj.eta0) ;
                    COST = COST * (1. + exp(0.018 * obj.Tin - 1.207 * 26.4)) ;
                    COST = COST * CEind(curr) / CEind(1995) ;
                    
                case 42
                    % Turboexpander costs based on quotations for CAES.
                    % From CEC-500-2008-069 CAES scoping report.
                    
                    COST = 160 * obj.W0/1e3  ;
                    COST = COST * CEind(curr) / CEind(2008) ;
                    
                case 43
                    % 'Turbine' from Couper
                    % Convert power to horsepower
                    P = obj.W0 * 0.00134102 ;
                    
                    nC = 1 ;
                    maxP = 8e3 ;
                    if P > maxP
                        nC = P / maxP ;
                        P  = maxP ;
                    end
                    
                    COST = nC * 1.1e3 * P ^ 0.81 ;
                    COST = COST * CEind(curr) / CEind(2010) ;
                
                case 44
                    % Use Agazzini's (Valero's) correlation, but scale it by the density
                    % of the outlet air (if higher density, then more
                    % compact expander)
                    
                    COST = 1.051 * 266.3 * obj.mdot0 * log(obj.pr0) / (0.92 - obj.eta0) ;
                    COST = COST * (1. + exp(0.018 * obj.Tin - 1.207 * 26.4)) ;
                    COST = COST * CEind(curr) / CEind(1995) ;
                    
                    % Estimate gamma
                    if strcmp(fld.read,'IDL')
                        Gama = fld.IDL.gam ;
                    else
                        Gama = CP1('PT_INPUTS',obj.Pin,obj.Tin,'CPMASS',fld.handle)/CP1('PT_INPUTS',obj.Pin,obj.Tin,'CVMASS',fld.handle);
                    end
                    
                    Pout = obj.Pin / obj.pr0 ;
                    Tout = obj.Tin / (obj.pr0 ^ (obj.eta0 * (Gama - 1)/Gama)) ;
                    
                    RHOout = RP1('PT_INPUTS',Pout,Tout,'D',fld) ;
                    RHO0  = 1.225 ; % Density of air at standard conditions
                    n     = 1 ; % This allows you to modify how strongly the cost is reduced
                    scale = (RHOout / RHO0)^n ;
                    if isinf(scale)
                        scale = 1 ;
                    end
                    
                    COST = COST / scale ; 
                    
                case 50 
                    % Cost for an sCO2 expander developed by Carlson et al. 2017
                    % See equation E4 in Q2 solar-PTES report
                    COST = 7790 * (obj.W0/1e3) ^ 0.6842 ;
                    COST = COST * CEind(curr) / CEind(2017) ;
                    
                case 51
                    % sCO2 expander from Morandin 2013 - eq. E5
                    COST = 40e3 + 9e3 * (obj.W0/1e3)^0.69 ;
                    COST = COST * CEind(curr) / CEind(2013) ;
                    
                case 52
                    % sCO2 radial turbine, 8-35 MW
                    % From Weiland et al 2019
                    % Valid up to TIT = 700C, and Pin = 200-260 bar
                    if obj.Tin > 550 + 273.15
                        fact = 1 + 1.137e-5 * (obj.Tin - 550 - 273.15)^2 ;
                    else
                        fact = 1.0 ;
                    end
                    COST = fact * 4.062e6 * (obj.W0/1e6)^0.8 ;
                    COST = COST * CEind(curr) / CEind(2019) ;
                    
                case 53
                    % sCO2 axial turbine, 10-750 MW
                    % From Weiland et al 2019
                    % Valid up to TIT = 730C, and Pin = 240-280 bar
                    if obj.Tin > 550 + 273.15
                        fact = 1 + 1.106e-4 * (obj.Tin - 550 - 273.15)^2 ;
                    else
                        fact = 1.0 ;
                    end
                    COST = fact * 1.826e6 * (obj.W0/1e6)^0.5561 ;
                    COST = COST * CEind(curr) / CEind(2019) ;
                    
                    
                case 60
                    % Pump cost. Centrifugal, carbon steel, includes motor.
                    % - eq. P1
                    COST = 2409.6 + 75.9 * obj.W0 / 1e3 ;
                    COST = COST * CEind(curr) / CEind(1998) ;
                    
                case 61
                    % Pump cost. Centrifugal, includes motor.
                    % - eq. P2
                    COST = 1227.5 + 177.8 * obj.W0 / 1e3 ;
                    COST = COST * CEind(curr) / CEind(1990) ;
                    
                case 62
                    % Pump cost. - eq. P3
                    COST = 50e3 + 1500 * (obj.W0 / 1e3)^0.8 ;
                    COST = COST * CEind(curr) / CEind(2009) ;
                    
                case 70
                    % FANS --> NEED UPDATING! (Just using value from 62)
                    COST = 50e3 + 1500 * (obj.W0 / 1e3)^0.8 ;
                    COST = COST * CEind(curr) / CEind(2009) ;
                    
                case 71
                    % FANS from Couper
                    % Convert mass flow to KSCFM (kilo-standard-cubic feet
                    % per min)
                    V = obj.mdot0 * 60 / 1.225 ; % m3/min based on air at standard conditions
                    V = V * 3.28084^3 / 1e3 ; % kilo-ft3/min
                    
                    nF = 1 ;
                    maxV = 500 ;
                    if V > maxV
                        nF = V / maxV ;
                        V  = maxV ;
                    end
                    
                    mat = 'CS' ;
                    
                    switch mat
                        case 'CS'
                            fm = 2.2 ;
                        case 'Fibreglass'
                            fm = 4.0 ;
                        case 'SS'
                            fm = 5.5 ;
                        case 'Ni'
                            fm = 11 ;
                    end
                    
                    fantype = 'propeller' ;
                    
                    switch fantype
                        case 'radial'
                            a = 0.4682 ;
                            b = 0.1203 ;
                            c = 0.0931 ;
                        case 'backward'
                            a = 0.04 ;
                            b = 0.1821 ;
                            c = 0.0786 ;
                        case 'propeller'
                            a = -0.4456 ;
                            b = 0.2211 ;
                            c = 0.082 ;
                        case 'prop vanes'
                            a = -1.0181 ;
                            b = 0.3332 ;
                            c = 0.0647 ;
                    end
                    
                    fp = 1 ;
                    
                    COST = 1218 * nF * fm * fp * exp(a + b * log(V) + c * (log(V)^2)) ;
                    COST = COST * CEind(curr) / CEind(2010) ;
                    
                case 72
                    % Fans from Benato 2017 - very high
                    COST = exp(6.6547 + 0.79 * log(obj.W0)) ;
                    COST = COST * CEind(curr) / CEind(2017) ;

            end
                                    
            obj.cmpexp_cost.COST = COST ;
            obj.cmpexp_cost.cost = COST / obj.W0 ;
            
        end
        
        
    end
end