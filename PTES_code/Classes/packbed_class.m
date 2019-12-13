classdef packbed_class
    properties
        type    % 
        
        L
        D
        A
        dp
        eps
        Sv
        AR
        MS
        VS
        Vtank
        
        rhoS
        kS
        cS
        
        ceff
        
        rhoF
        kF
        cF
        mu
        Pr
        mdot
        us
        ui
        u
        Pin
        
        tN
        dt
        dx
        DELX
        CFL
        Tprof
        Nprof
        NX
        Nt
        TMAX
        Ncyc
        
        Ltime
        Ltext
        timeC
        timeD
        textC
        textD
        
        keff
        Bi
        Re0
        St0
        Nu0
        h0
        
        Cond
        
        Cf
        
        TC
        TD
        
        TS
        TF
        P
        rho
            
        w       % specific work transfer, J/kg
        q       % specific heat transfer, J/kg
        Dh      % specific enthalpy change, J/kg
        sirr    % specific entropy generation, J/kgK
        
        W        % work transfer, W
        Q        % heat transfer, W
        DH       % enthalpy change, W
        Sirr     % entropy generation, W/K
                
        % Economics
        pckbed_cost = econ_class(0,0,0,0) ;
        
    end
    
    methods
        function obj = packbed_class(type)
            
            obj.type = type ;
            
        end
        
        % Set up packed bed geometry and time steps etc.
        function obj = PB_INITIALISE(obj)
            
            % mass of solid required
            obj.ceff = obj.eps * obj.rhoF * obj.cF + (1-obj.eps) * obj.rhoS * obj.cS ;
            obj.ceff = obj.ceff / (obj.rhoF + obj.rhoS) ;
            obj.MS = obj.tN * obj.mdot * obj.cF / obj.cS ;
            
            % volumes
            obj.VS    = obj.MS / obj.rhoS ; % Solid volume
            obj.Vtank = obj.VS / obj.eps ; % Tank volume
            
            % Geometry
            obj.D  = (4. * obj.Vtank / (pi * obj.AR)) ^ (1./3.) ;
            obj.L  = obj.D * obj.AR ;
            obj.A  = 0.25 * pi * obj.D^2 ;
            obj.Sv = 6. / obj.dp ;
            
            % Velocities
            obj.us = obj.mdot / (obj.rhoF * obj.A) ; % Superficial velocity
            obj.ui = obj.us / obj.eps ; % Interstitial velcoty
            
            % Grid and time-steps
            obj.dx = obj.DELX * obj.L ;
            obj.dt = obj.CFL * obj.dx / obj.ui ;
            obj.NX = obj.L / obj.dx ; % Number of grid steps
            obj.Nt = obj.TMAX * obj.tN / obj.dt ; % Max number of time steps
            obj.TMAX = obj.TMAX * obj.tN ;
            
            % Save temp profiles at the following times
            if obj.Nprof < 3
                obj.Nprof = 3 ;
            end
            obj.Tprof = ones(obj.Nprof,1) ; % Always save profiles at the start and end, so have Nprof-2 profiles to find the times for
            inc       = int64(obj.tN / (obj.dt * (obj.Nprof-1))) ;
            for i = 2 : obj.Nprof-1
                obj.Tprof(i) = i * inc ;
            end
            
            % Effective conductivity
            obj.keff = eff_cond(obj) ;
            
            % Coefficienct of friction
            obj.Cf = friction(obj, obj.rhoF, obj.us) ;
            
            % Stanton number
            obj.Re0 = obj.rhoF * obj.us * obj.dp / obj.mu ;
            [obj.Nu0, obj.h0] = nusselt(obj.Re0, obj.Pr, obj) ;
            obj.St0 = obj.Nu0 / (obj.Re0 * obj.Pr) ;
            
            % Biot number
            obj.Bi = obj.St0 * obj.dp * obj.rhoF * obj.us * obj.cF / obj.kS ;
            
            % Create arrays for properties that are discretized.
            % Some have two columns. One for the current time-step, and one for the previous time-step
            obj.TS  = zeros(obj.NX, 2) ; % Solid temperature
            obj.TF  = zeros(obj.NX, 2) ; % Fluid temperature
            obj.P   = zeros(obj.NX, 1) ; % Pressure
            obj.u   = zeros(obj.NX, 2) ; % Fluid velocity
            obj.rho = zeros(obj.NX, 2) ; % Fluid density
            
            % Set the initial temperature distributions
            obj.TS(:,:)  = obj.TD ;
            obj.TF(:,:)  = obj.TD ;
            obj.P        = obj.Pin ;
            obj.u(:,:)   = obj.us ; 
            obj.rho(:,:) = obj.rhoF ; 
            
            % Error control
            if (obj.Ltime && obj.Ltext)
                error('Packed beds: Ltime and Ltext cannot both be true or false')
            end
            
        end
        
        
        % March forward one timestep
        function obj = PB_TIMESTEP(obj, mode)
            
            switch mode
                case 'chg'
                    Tin = obj.TC ;
                case 'dis'
                    Tin = obj.TD ;
            end
            
            N  = obj.NX ;
            Tf = obj.TF ;
            Ts = obj.TS ;
            
            % Find the heat capacities at each grid step as a function of
            % temperature (to be implemented)
            CF = obj.cF .* ones(N,1) ;
            CS = obj.cS .* ones(N,1) ;
            
            % Calculate Re, Nu, and h for each grid step
            Re = obj.rhoF .* obj.u(:,2) .* obj.dp ./ obj.mu ;
            Nu = zeros(N,1) ;
            h  = zeros(N,1) ;
            for i = 1 : N
                [Nu(i), h(i)] = nusselt(Re(i), obj.Pr, obj) ;
            end
            
            % First calculate some coefficients at each node
            a = obj.u(:,2) .* obj.dt / (obj.eps * obj.dx) ;
            b = obj.keff .* obj.dt ./ (obj.eps .* obj.rhoF .* CF * obj.dx * obj.dx) ;
            c = (1-obj.eps) .* obj.Sv .* h .* obj.dt ./ (2 .* obj.eps .* obj.rhoF .* CF) ;
            
            d = obj.Sv .* h .* obj.dt ./ (2 .* obj.rhoS .* CS) ;
            
            % Calculate the determinant at each node
            DET = 1 ./ ((1 + d) .* (1 + a + c) - c .* d ) ;
            
            % terms for matrix calculation
            t1 = (1 - d) .* Ts(:,2) + d .* Tf(:,2) ;
            t2 = (1 - c) .* Tf(:,2) + c .* Ts(:,2) - 2 .* b .* Tf(:,2) ;
            
            % First node
            X = [(1+a(1)+c(1)), d(1) ; c(1), (1+d(1))] ; % Could probably make a matrix of matrices, urgh
            Y = [ t1(1) ; t2(1) + a(1) * Tin + b(1) * (Tf(2,2) + Tin)] ;
            T(1,:) = [X * Y * DET(1)]' ;
            
            % Remaining nodes
            for i = 2 : N-1
                X = [(1 + a(i) + c(i)), d(i) ; c(i), (1 + d(i))] ; % Could probably make a matrix of matrices, urgh
                Y = [ t1(i) ; t2(i) + a(i) * Tf(i-1,1) + b(i) * (Tf(i+1,2) + Tf(i-1,2))] ;
                T(i,:) = [X * Y * DET(i)]' ;
            end
            
            % Final node
            X = [(1+a(N)+c(N)), d(N) ; c(N), (1+d(N))] ; % Could probably make a matrix of matrices, urgh
            Y = [ t1(N) ; t2(N) + a(N) * Tf(N-1,1) + b(1) * (3*Tf(N,2) - 2.*Tf(N-1,2) + Tf(N-2,2))] ;
            T(N,:) = [X * Y * DET(N)]' ;
            
            % Assign temperatures
            obj.TS(:,1) = T(:,1) ;
            obj.TF(:,1) = T(:,2) ;
            
            % Calculate densities - need to call Coolprop
            obj.rho(:,1) = obj.rho(:,2) ; % temporary call
            
            % Calculate velocities and pressures 
            obj.u(1,1) = obj.us ;
            obj.P(1)   = obj.Pin ;
            fact1 = obj.Sv * (1-obj.eps) * obj.dx / (2. * obj.eps^3) ;
            for i = 2 : N
               obj.u(i,1) = (obj.rho(i-1,1) * obj.u(i-1,1) + obj.eps * obj.dx * (obj.rho(i,2) - obj.rho(i,1))) / obj.rho(i,1) ;
               fact2      = fact1 * obj.rho(i,1) * obj.u(i,1) * obj.u(i,1) ;
               obj.P(i)   = obj.P(i-1) - fact2 * friction(obj,obj.rho(i,1),obj.u(i,1)) ; 
            end
            
            
            % Should now probably calculate lots of losses
                        
            % Move new results into previous results, for use in future steps
            obj.TS(:,2)  = obj.TS(:,1) ;
            obj.TF(:,2)  = obj.TF(:,1) ;
            obj.u(:,2)   = obj.u(:,1) ;
            obj.rho(:,2) = obj.rho(:,1) ;
            
            
        end
        
        % Calculate energy in storage at end of charge and end of discharge
        function [obj, en] = PB_ENERGY(obj)
            
            en = obj.cS * obj.rhoS * (1-obj.eps) * (trapz(obj.dx, obj.TS(:,1)) - obj.TD) ;
            en = en + obj.cF * obj.rhoF * obj.eps * (trapz(obj.dx, obj.TF(:,1)) - obj.TD) ;
            en = en * obj.A / 3600e6 ; % Energy in MWh
            
        end
        
        
        
        
        function keff = eff_cond(obj)
           
            e = obj.eps ;
            kf = obj.kF ;
            ks = obj.kS ;
            % Very simply assume that the effective conductivity is the
            % average of the minimum and maximum estimates
            k1 = kf / (e + (1-e)*kf/ks) ;
            k2 = kf * (e + (1-e)*ks/kf) ;
            
            keff = 0.5 * (k1 + k2) * obj.Cond ;
            
        end
        
        % Coefficient of friction - Carman correlation
        function Cf = friction(obj, rhoF, us)
            
            Rem = obj.rhoF * obj.us / ((1-obj.eps) * obj.Sv * obj.mu) ;
            Cf  = 10 / Rem + 8 / (10 * Rem^(1/10)) ;
            
        end
        
        % Stanton number
        function [Nu, h] = nusselt(Re, Pr, obj)
           
            Nu = 2. + 1.1 * Re^(3/5) * Pr^(1/3) ;
            h  = Nu * obj.kF / obj.dp ;
            
        end
        
        
        % Calculate the economic cost of the packed bed
        function obj = pckbed_econ(obj, CEind)
            mode = obj.pckbed_cost.cost_mode ;
            curr = 2019 ; % Current year
            
            switch mode
                case 0
                    % Random correlation
                    COST = 200. * obj.vol ;
                    COST = COST * CEind(curr) / CEind(2018) ;
                
                case 1
                    error('Not implemented')

            end
              
            obj.pckbed_cost.COST = COST ;
            obj.pckbed_cost.cost = COST / obj.vol ;
            
        end
        
     
        
        
    end
end