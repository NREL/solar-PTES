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
        V       % Tank volume
        
        Sname % Solid name
        sld   % Table that contains solid data
        rhoS
        kS
        cS
                
        rhoF
        kF
        cF
        R
        mu
        Pr
        mdot
        us
        ui
        u
        Pin
        
        CSfacs
        CFfacs  % This is a structure of coefficients
        Hfacs   % coefficients for enthalpy
        Dfacs   % coefficients for density
        Sfacs   % coefficients for entropy
        KFfacs  % coefficients for conductivity
        MUfacs  % coefficients for viscosity
        
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
        
        Ltime % End cycle after a certain time?
        Ltext % End cycle after a certain temperature fraction has been exceeded?
        timeC % Time after which to end charging cycle 
        timeD % Time after which to end discharging cycle
        textC % Exit temperature fraction during charge
        textD % Exit temperature fraction during discharge
        TxC   % Exit temperature during charge
        TxD   % Exit temperature during discharge
        time
        
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
        
        % These are profiles along the axis that are used for calculations
        TS
        TF
        P
        rho
        
        % These are profiles that are saved. Array with 3 dimensions.
        % 1: gridsteps. 2: i-th profile in that cycle. 3: cycle number
        TSprof % Solid temp
        TFprof % Fluid temp
        
        H       % Enthalpy stored in the packed bed (at the end of each cycle)
        S       % Entropy of the packed bed (at the end of each cycle)
        Mflux   % Flux of mass into/out of the packed bed
        Hflux   % Flux of enthalpy into/out of the packed bed
        Sflux   % Flux of entropy into/out of the packed bed
        Hnet    % Net change of enthalpy (should be zero)
        Snet    % Net change of entropy (is sirr)
            
        w       % specific work transfer, J/kg
        q       % specific heat transfer, J/kg
        Dh      % specific enthalpy change, J/kg
        sirr    % specific entropy generation, J/kgK
        
        W        % work transfer, W
        Q        % heat transfer, W
        DH       % enthalpy change, W
        Sirr     % entropy generation, W/K
        
        DHprev   % Save the previous cycles DH here (for assessing whether steady state has been reached)
        Sirrprev % Save the previous cycles Sirr here (for assessing whether steady state has been reached)
        
        Lconst
        Lideal
                
        % Economics
        pckbed_cost = econ_class(0,0,0,0) ;
        
    end
    
    methods
        function obj = packbed_class(type)
            
            obj.type = type ;
            
        end
        
        % Set up packed bed geometry and time steps etc.
        function obj = PB_INITIALISE(obj,fld,p,load)
            
            obj.sld   = create_table(obj.Sname) ;
            obj.kS    = obj.sld(1,6) ;    % Thermal conductivity, W/mK
            obj.rhoS  = obj.sld(1,3) ;   % Density, kg/m3
            
            x  = obj.sld(:,1); %temperatures
            h1 = interp1(x,obj.sld(:,2),obj.TC);
            h2 = interp1(x,obj.sld(:,2),obj.TD);
            
            obj.cS   = (h1 - h2) / (obj.TC - obj.TD) ;    % Specific heat capacity, J/kgK
            
            % Fluid properties
            obj.mdot = load.mdot(1) ; % Fluid mass flow rate, kg/s
            obj.Pin  = p ;
            Tave     = 0.5 * (obj.TC + obj.TD) ;
            obj.kF   = RP1('PT_INPUTS',obj.Pin,Tave,'L',fld); %0.035 ; % Thermal conductivity, W/mK
            obj.rhoF = RP1('PT_INPUTS',obj.Pin,Tave,'D',fld);%857;%12;% % Density, kg/m3 - Need to calculate these properly!
            obj.cF   = RP1('PT_INPUTS',obj.Pin,Tave,'CPMASS',fld); %2300;%1000; % Specific heat capacity, J/kgK
            obj.Pr   = RP1('PT_INPUTS',obj.Pin,Tave,'PRANDTL',fld);%0.7 ;  % Prandtl number
            obj.mu   = RP1('PT_INPUTS',obj.Pin,Tave,'V',fld);%5e-4 ;%1e-5;% Viscosity, Pa.s
            
            % volume of solid required
            obj.V = obj.tN * obj.mdot * obj.cF / ((1-obj.eps)*obj.rhoS * obj.cS) ;
                         
            % Geometry
            obj.D  = (4. * obj.V / (pi * obj.AR)) ^ (1./3.) ;
            obj.L  = obj.D * obj.AR ;
            obj.A  = 0.25 * pi * obj.D^2 ;
            obj.Sv = 6. / obj.dp ;
            
            % Velocities
            obj.us = obj.mdot / (obj.rhoF * obj.A) ; % Superficial velocity
            obj.ui = obj.us / obj.eps ; % Interstitial velcoty
            
            % Grid and time-steps
            obj.dx = obj.DELX * obj.L ;
            obj.dt = obj.CFL * obj.dx / obj.ui ;
            obj.NX = int64(obj.L / obj.dx) ; % Number of grid steps
            obj.Nt = obj.TMAX * obj.tN / obj.dt ; % Max number of time steps
            obj.TMAX = obj.TMAX * obj.tN ; % Maximum time to allow cycle to run for 
            
            % Exit temperature fractions
            obj.TxC = obj.TD + obj.textC * (obj.TC - obj.TD) ;
            obj.TxD = obj.TC + obj.textD * (obj.TD - obj.TC) ;
            
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
            obj.keff = 0.0;%eff_cond(obj) ;
            
            % Coefficienct of friction
            obj.Cf = friction(obj, obj.rhoF, obj.us) ;
            
            % Dimensionless numbers
            obj.Re0 = obj.rhoF * obj.us * obj.dp / obj.mu ;
            [obj.Nu0, obj.h0] = nusselt(obj.Re0, obj.Pr, obj) ;
            obj.St0 = obj.Nu0 / (obj.Re0 * obj.Pr) ;
            
            % Biot number
            obj.Bi = obj.St0 * obj.dp * obj.rhoF * obj.us * obj.cF / obj.kS ;
            
            % Create arrays for properties that are discretized.
            % Two columns. One for the current time-step, and one for the previous time-step
            obj.TS  = zeros(obj.NX, 2) ; % Solid temperature
            obj.TF  = zeros(obj.NX, 2) ; % Fluid temperature
            obj.P   = zeros(obj.NX, 2) ; % Pressure
            obj.u   = zeros(obj.NX, 2) ; % Fluid velocity
            obj.rho = zeros(obj.NX, 2) ; % Fluid density
            
            % Set the initial temperature distributions
            obj.TS(:,:)  = obj.TD ;
            obj.TF(:,:)  = obj.TD ;
            obj.P(:,:)   = obj.Pin ;
            obj.u(:,:)   = obj.us ; 
            obj.rho(:,:) = obj.rhoF ; 
            
            % Use a simple fit to find correlation for temperature
            % dependent properties - FLUID
            if strcmp(fld.read,'TAB')
                lo = min(obj.TD,obj.TC); % Lowest temp
                hi = max(obj.TD,obj.TC); % Highest temp
                
                % Find indices closest to these hi and lo temps
                [~, lo] = min( abs( fld.TAB(:,1) - lo ) ) ;
                [~, hi] = min( abs( fld.TAB(:,1) - hi ) ) ;
                
                % Only find the heat cap function between these temperatures
                X = [ones(size(fld.TAB(lo:hi,1))) fld.TAB(lo:hi,1) fld.TAB(lo:hi,1).*fld.TAB(lo:hi,1) ];
                obj.Hfacs  = X\fld.TAB(lo:hi,2); % coefficients for enthalpy
                obj.Dfacs  = X\fld.TAB(lo:hi,3); % coefficients for density
                obj.Sfacs  = X\fld.TAB(lo:hi,4); % coefficients for entropy
                obj.CFfacs = X\fld.TAB(lo:hi,5); % coefficients for heat capacity
                obj.KFfacs = X\fld.TAB(lo:hi,6); % coefficients for conductivity
                obj.MUfacs = X\fld.TAB(lo:hi,7); % coefficients for viscosity
            elseif strcmp(fld.read,'CP')
                if obj.Lideal
                    obj.CFfacs = [obj.cF; 0; 0] ;
                    obj.R      = RP1('PT_INPUTS',obj.Pin,Tave,'CPMASS',fld) - RP1('PT_INPUTS',obj.Pin,Tave,'CVMASS',fld) ;
                else
                    error('Not implemented')                    
                end
            else
                error('Not implemented')
            end
            
            % Use a simple fit to find correlation for temperature
            % dependent properties - SOLID
            lo = min(obj.TD,obj.TC); % Lowest temp
            hi = max(obj.TD,obj.TC); % Highest temp
                
            % Find indices closest to these hi and lo temps
            [~, lo] = min( abs( obj.sld(:,1) - lo ) ) ;
            [~, hi] = min( abs( obj.sld(:,1) - hi ) ) ;
                
            % Only find the heat cap function between these temperatures
            X = [ones(size(obj.sld(lo:hi,1))) obj.sld(lo:hi,1) obj.sld(lo:hi,1).*obj.sld(lo:hi,1) ];
            obj.CSfacs = X\obj.sld(lo:hi,5); % coefficients for heat capacity
                
            obj.Mflux = zeros(30,2) ;
            obj.Hflux = zeros(30,2) ;
            obj.Sflux = zeros(30,2) ;
            obj.time  = zeros(30,2) ;
            
            
            % Error control
            if (obj.Ltime && obj.Ltext)
                error('Packed beds: Ltime and Ltext cannot both be true or false')
            end
            
        end
        
        % This function obtains the timesteps of each packed bed, picks the
        % shortest, and resets values for the other beds
        function [hot, cold] = PB_TIMINGS(hot, cold)
            
            % Set minimum time step
            mindt = min(hot.dt,cold.dt) ;
            
            hot.dt = mindt;
            cold.dt = mindt ;
            
            % Hot
            hot.Nt = hot.TMAX / hot.dt ; % Max number of time steps
                                
            % Save temp profiles at the following times
            inc       = int64(hot.tN / (hot.dt * (hot.Nprof-1))) ;
            for i = 2 : hot.Nprof-1
                hot.Tprof(i) = i * inc ;
            end
            
            % Cold
            cold.Nt = cold.TMAX / cold.dt ; % Max number of time steps
                                
            % Save temp profiles at the following times
            inc       = int64(cold.tN / (cold.dt * (cold.Nprof-1))) ;
            for i = 2 : cold.Nprof-1
                cold.Tprof(i) = i * inc ;
            end
            
        end
        
        function Lend = PB_STOP_PHASE(hot, cold, time, mode)
           
            Lend = false ;
            
            c = zeros(4,1) ;
            c(1) = time > hot.TMAX ;
            c(2) = time > cold.TMAX ;
                        
            switch mode
                case 'chgPB'
                    
                    if hot.Ltime
                        % Has time elapsed?
                        c(3) = time > hot.timeC ;
                        c(4) = time > cold.timeC ;
                    elseif hot.Ltext
                        % Have exit temperatures been exceeded?
                        c(3) = hot.TS(end,1) > hot.TxC ;
                        c(4) = cold.TS(end,1) < cold.TxC ;
                    end
                    
                case 'disPB'
                
                    if hot.Ltime
                        % Has time elapsed?
                        c(3) = time > hot.timeD ;
                        c(4) = time > cold.timeD ;                
                    elseif hot.Ltext
                        % Have exit temperatures been exceeded?
                        c(3) = hot.TS(end,1) < hot.TxD ;
                        c(4) = cold.TS(end,1) > cold.TxD ;
                    end
                    
            end
            
            % Check if any conditions are true
            if any(c)
                Lend = true ;
            end
            
            
        end
        
        function [hot, cold] = PB_REVERSE(hot, cold, Load, i)
            % If time to end the phase then switch around arrays if the
            % follow load.type is different to the current load type
            % Watch out for being at the final load
            if i < Load.num
                if Load.type(i) ~= Load.type(i+1)
                    
                    hot.TS  = flip(hot.TS,1) ;
                    hot.TF  = flip(hot.TF,1) ;
                    hot.P   = flip(hot.P,1) ;
                    hot.u   = flip(hot.u,1) ;
                    hot.rho = flip(hot.rho,1) ;
                    
                    cold.TS  = flip(cold.TS,1) ;
                    cold.TF  = flip(cold.TF,1) ;
                    cold.P   = flip(cold.P,1) ;
                    cold.u   = flip(cold.u,1) ;
                    cold.rho = flip(cold.rho,1) ;
                    
                    % Set the flux counters to zero after each cycle
                    hot.Mflux = zeros(30,2) ;
                    hot.Hflux = zeros(30,2) ;
                    hot.Sflux = zeros(30,2) ;
                    
                    cold.Mflux = zeros(30,2) ;
                    cold.Hflux = zeros(30,2) ;
                    cold.Sflux = zeros(30,2) ;
                    
                end
            end
            
        end
        
        % Check to see whether the hot and cold packed beds have changed
        % significantly compared to the previous cycle. If not, then cyclic
        % operation has probably been reached
        function [hot, cold, iL, Icyc] = PB_STOP_CYC(hot, cold, Load, iL, Icyc, Ncyc, Lcyclic)
            
            if Lcyclic
                error = 1e-2 ;
                if isempty(hot.DHprev)
                    hot.DHprev = hot.DH*100;
                    hot.Sirrprev = hot.Sirr*100;
                    cold.DHprev = cold.DH*100;
                    cold.Sirrprev = cold.Sirr*100;
                end
                
                if Icyc < Ncyc
                    cyclic = zeros(1,4) ;
                   
                    cyclic(1) = 100 * abs((hot.DH(1) - hot.DHprev(1)) / hot.DHprev(1)) < error ;
                    cyclic(2) = 100 * abs((cold.DH(1) - cold.DHprev(1)) / cold.DHprev(1)) < error ;
                    
                    StotH  = hot.Sirr(1) + hot.Sirr(2) ;
                    StotHp = hot.Sirrprev(1) + hot.Sirrprev(2) ;
                    
                    StotC  = cold.Sirr(1) + cold.Sirr(2) ;
                    StotCp = cold.Sirrprev(1) + cold.Sirrprev(2) ;
                                        
                    cyclic(3) = 100 * abs((StotH - StotHp) / StotHp) < error ;
                    cyclic(4) = 100 * abs((StotC - StotCp) / StotCp) < error ;
                    
                    if all(cyclic)
                        fprintf("\n\nSTEADY STATE OPERATION REACHED!\n\n");
                        iL=iL+1;
                    else
                        Icyc = Icyc + 1 ;
                        iL   = 1 ;
                        
                        % Flip arrays over
                        [hot, cold] = PB_REVERSE(hot, cold, Load, iL);
                    end
                else
                    warning('Maximum number of cycles (%i) has elapsed without reaching steady-state operation',Ncyc);
                    iL=iL+1;
                end
                
                hot.DHprev   = hot.DH ;
                hot.Sirrprev = hot.Sirr ;
                
                cold.DHprev   = cold.DH ;
                cold.Sirrprev = cold.Sirr ;
                
                
            end
            
        end
        
        function obj = PB_FLUX(obj, T0, P0, fld, iCYC)
            
            % Reference points
            HF0 = RP1('PT_INPUTS',P0,T0,'H',fld) ;
            SF0 = RP1('PT_INPUTS',P0,T0,'S',fld) ;

            % Calculate the mass, enthalpy, and entropy flux INTO the storage in that timestep
            Min = obj.u(1,1) * obj.rho(1,1) * obj.A * obj.dt ;
            obj.Mflux(iCYC, 1) = obj.Mflux(iCYC, 1) + Min ;
            obj.Hflux(iCYC, 1) = obj.Hflux(iCYC, 1) + Min * (RP1('PT_INPUTS',obj.P(1),obj.TF(1,1),'H',fld) - HF0);
            %obj.Hflux(iCYC, 1) = obj.Hflux(iCYC, 1) + Min * (obj.cF *(obj.TF(1,1)-T0));
            obj.Sflux(iCYC, 1) = obj.Sflux(iCYC, 1) + Min * (RP1('PT_INPUTS',obj.P(1),obj.TF(1,1),'S',fld) - SF0);
            
            % Calculate the mass, enthalpy, and entropy flux OUT OF the storage in that timestep
            Mout = obj.u(end,1) * obj.rho(end,1) * obj.A * obj.dt ;
            obj.Mflux(iCYC, 2) = obj.Mflux(iCYC, 2) + Mout ;
            obj.Hflux(iCYC, 2) = obj.Hflux(iCYC, 2) + Mout * (RP1('PT_INPUTS',obj.P(end),obj.TF(end,1),'H',fld) - HF0) ;
            %obj.Hflux(iCYC, 2) = obj.Hflux(iCYC, 2) + Mout * (obj.cF *(obj.TF(end,1)-T0));
            obj.Sflux(iCYC, 2) = obj.Sflux(iCYC, 2) + Mout * (RP1('PT_INPUTS',obj.P(end),obj.TF(end,1),'S',fld) - SF0) ;
             
            
        end
        
        % Energy and exergy balances for phase iCYC
        function [obj] = PB_LOSS(obj, iCYC)
            
            obj.W(iCYC)    = 0.0 ;
            obj.Q(iCYC)    = obj.Hflux(iCYC,1) - obj.Hflux(iCYC,2) ; % Heat transferred into/out of storage
            obj.DH(iCYC)   = obj.H(iCYC,1) - obj.H(iCYC,2) ; % Change in enthalpy of packed bed. Should equal Q (small errors occur though)
            obj.Sirr(iCYC) = obj.S(iCYC,1) - obj.S(iCYC,2) ; % Change in entropy of packed bed
            obj.Sirr(iCYC) = obj.Sirr(iCYC) - (obj.Sflux(iCYC,1) - obj.Sflux(iCYC,2)) ; % Net entropy increase of packed bed
            obj.Sirr(iCYC) = abs(obj.Sirr(iCYC)) ; % Has to be positive
            
        end
             
        % Calculate energy in storage at end of charge and end of discharge
        function [obj] = PB_ENERGY(obj, fld, iCYC, mode)
            
            N  = obj.NX ;
            Tf = obj.TF ;
            Ts = obj.TS ;
            
            HF  = zeros(N,1) ;
            SF  = zeros(N,1) ;
            
            HS  = zeros(N,1) ;
            SS  = zeros(N,1) ;
            
            if obj.TC > obj.TD
                T0 = obj.TD ;
            else
                T0 = obj.TC ;
            end
            
            if obj.Lconst
                CF  = obj.cF .* ones(N,1) ;
                CS  = obj.cS .* ones(N,1) ;
                
                en = obj.rhoS * (1-obj.eps) * (trapz(obj.dx, CS.*(obj.TS(:,1)- T0)) ) ;
                en = en + obj.eps * (trapz(obj.dx, CF.*obj.rho(:,1).*(obj.TF(:,1)- T0)) ) ;
                en = en * obj.A ; % Energy in MWh
                
                ent = 0.0;
            
            else
            
                % Calculate solid and fluid entropies and entropies by
                % interpolating original data
                
                % Reference points
                HF0 = RP1('PT_INPUTS',1e5,T0,'H',fld) ;
                SF0 = RP1('PT_INPUTS',1e5,T0,'S',fld) ;
                
                x = obj.sld(:,1) ;
                HS0 = interp1(x,obj.sld(:,2),T0) ;
                SS0 = interp1(x,obj.sld(:,4),T0) ;
                
                for i = 1 : N
                    
                    HF(i) = RP1('PT_INPUTS',1e5,obj.TF(i,1),'H',fld) ; % Enthalpy
                    SF(i) = RP1('PT_INPUTS',1e5,obj.TF(i,1),'S',fld) ; % Entropy
                    
                    HS(i) = interp1(x,obj.sld(:,2),obj.TS(i,1)) ; % Enthalpy
                    SS(i) = interp1(x,obj.sld(:,4),obj.TS(i,1)) ; % Entropy
                    
                end
                
                en = obj.rhoS * (1-obj.eps) * (trapz(obj.dx, (HS(:) - HS0))) ;
                en = en + obj.eps * (trapz(obj.dx, obj.rho(:,1).*(HF(:) - HF0))) ;
                en = en * obj.A ; % Energy 
                
                ent = obj.rhoS * (1-obj.eps) * (trapz(obj.dx, (SS(:) - SS0))) ;
                ent = ent + obj.eps * (trapz(obj.dx, obj.rho(:,1).*(SF(:) - SF0))) ;
                ent = ent * obj.A  ; % Entropy 
                                                
            end
            
            % Assign energy and entropy to correct term
            switch mode
                case 'start'
                    obj.H(iCYC,1) = en ;
                    obj.S(iCYC,1) = ent ;
                case 'end'
                    obj.H(iCYC,2) = en ;
                    obj.S(iCYC,2) = ent ;
            end
            
        end
        
        
        % March forward one timestep
        function obj = PB_TIMESTEP(obj, fld, mode)
            
            N  = obj.NX ;
            Tf = obj.TF ;
            Ts = obj.TS ;
                        
            % Some arrays
            Nu  = zeros(N,1) ;
            h   = zeros(N,1) ;
            CF  = zeros(N,1) ;
            CS  = zeros(N,1) ;
            %muf = zeros(N,1) ;
            
            switch mode
                case 'chg'
                    Tin = obj.TC ;
                case 'dis'
                    Tin = obj.TD ;
            end
            
            muf = obj.mu .* ones(N,1) ;
            % Find properties that vary with temperature
            if obj.Lconst
                CF  = obj.cF .* ones(N,1) ;
                CS  = obj.cS .* ones(N,1) ;
                muf = obj.mu .* ones(N,1) ;
            elseif strcmp(fld.read,'CP') && ~obj.Lideal
                for i = 1 : N
                    % Temperature array from previous timestep
                    Xs = [1 Ts(i,2) Ts(i,2)*Ts(i,2)] ;
                    CS(i) = Xs * obj.CSfacs ;
                    CF(i) = CP1('PT_INPUTS',obj.P(i,2),obj.TF(:,1),'CPMASS',fld.handle) ;
                    %muf(i) = Xf * obj.MUfacs ;                    
                end
            else % This accounts for fld.read='TAB' and also ideal gases
                for i = 1 : N
                    % Temperature array from previous timestep
                    Xf = [1 Tf(i,2) Tf(i,2)*Tf(i,2)] ;
                    Xs = [1 Ts(i,2) Ts(i,2)*Ts(i,2)] ;
                    CF(i) = Xf * obj.CFfacs ;
                    CS(i) = Xs * obj.CSfacs ;
                    %muf(i) = Xf * obj.MUfacs ;                    
                end
            end
            
            % Reynolds number, Nusselt number and heat transfer coef
            for i = 1 : N
                Re  = obj.rho(i,2) .* obj.u(i,2) .* obj.dp ./ muf(i) ;
                [Nu(i), h(i)] = nusselt(Re, obj.Pr, obj) ;
            end
            
            % First calculate some coefficients at each node
            a = obj.u(:,2) .* obj.dt / (obj.eps * obj.dx) ;
            b = obj.keff .* obj.dt ./ (obj.eps .* obj.rho(:,2) .* CF * obj.dx * obj.dx) ;
            c = (1-obj.eps) .* obj.Sv .* h .* obj.dt ./ (2 .* obj.eps .* obj.rho(:,2) .* CF) ;
            
            d = obj.Sv .* h .* obj.dt ./ (2 .* obj.rhoS .* CS) ;
            
            % Calculate the determinant at each node
            DET = ((1 + d) .* (1 + a + c) - c .* d ) ;
            
            % terms for matrix calculation
            t1 = (1 - d) .* Ts(:,2) + d .* Tf(:,2) ;
            t2 = (1 - c) .* Tf(:,2) + c .* Ts(:,2) - 2 .* b .* Tf(:,2) ;
            
            % Bit annoying Matlab doesn't do zero indices - otherwise the
            % first node could be easily wrapped into the for loop below
            % First node
            X = [(1+a(1)+c(1)), d(1) ; c(1), (1+d(1))] ; % Could probably make a matrix of matrices, urgh
            Y = [ t1(1) ; t2(1) + a(1) * Tin + b(1) * (Tf(2,2) + Tin)] ;
            T(1,:) = X * Y / DET(1) ;
            
            % For final node extrapolate Tf at N+1
            Tf(N+1,2) = 2.*Tf(N,2) - Tf(N-1,2) ;
            
            % Remaining nodes
            for i = 2 : N
                X = [(1 + a(i) + c(i)), d(i) ; c(i), (1 + d(i))] ; % Could probably make a matrix of matrices, urgh
                Y = [ t1(i) ; t2(i) + a(i) * Tf(i-1,1) + b(i) * (Tf(i+1,2) + Tf(i-1,2))] ;
                T(i,:) = X * Y / DET(i) ;
            end
            
            % Final node
            %X = [(1+a(N)+c(N)), d(N) ; c(N), (1+d(N))] ; % Could probably make a matrix of matrices, urgh
            %Y = [ t1(N) ; t2(N) + a(N) * Tf(N-1,1) + b(1) * (3*Tf(N,2) - 2.*Tf(N-1,2) + Tf(N-2,2))] ;
            %T(N,:) = X * Y / DET(N) ;
            
            % Assign temperatures
            obj.TS(:,1) = T(:,1) ;
            obj.TF(:,1) = T(:,2) ;
            
            % Calculate densities 
            if obj.Lconst
                obj.rho(:,1) = obj.rho(:,2)  ;
            else
                if strcmp(fld.read,'TAB')
                    for i = 1 : N
                        Xf = [1 T(i,2) T(i,2)*T(i,2) ] ;
                        obj.rho(i,1) = Xf * obj.Dfacs ;
                    end
                elseif obj.Lideal
                    for i = 1 : N
                        obj.rho(i,1) = obj.P(i,2) / (T(i,2) * obj.R) ;
                    end
                elseif strcmp(fld.read,'CP')
                    for i = 1 : N
                        obj.rho(i,1) = CP1('PT_INPUTS',obj.P(i,2),obj.TF(:,1),'D',fld.handle) ;
                    end
                end
            end
            
            
            % Calculate velocities and pressures 
            obj.u(1,1) = obj.us ;
            obj.P(1,1) = obj.Pin ;
            fact1 = obj.Sv * (1-obj.eps) * obj.dx / (2. * obj.eps^3) ;
            for i = 2 : N
               obj.u(i,1) = (obj.rho(i-1,1) * obj.u(i-1,1) + obj.eps * obj.dx * (obj.rho(i,2) - obj.rho(i,1))/obj.dt) / obj.rho(i,1) ;
               fact2      = fact1 * obj.rho(i,1) * obj.u(i,1) * obj.u(i,1) ;
               obj.P(i,1) = obj.P(i-1,1) - fact2 * friction(obj,obj.rho(i,1),obj.u(i,1)) ; 
            end
            
            % Move new results into previous results, for use in future steps
            obj.TS(:,2)  = obj.TS(:,1) ;
            obj.TF(:,2)  = obj.TF(:,1) ;
            obj.u(:,2)   = obj.u(:,1) ;
            obj.P(:,2)   = obj.P(:,1) ;
            obj.rho(:,2) = obj.rho(:,1) ;
            
            
        end
        
        
        % March forward one timestep
        % Modified to include pressure losses
        function obj = PB_TIMESTEP2(obj, fld, mode)
            
            N  = obj.NX ;
            Tf = obj.TF ;
            Ts = obj.TS ;
                        
            % Some arrays
            Nu  = zeros(N,1) ;
            h   = zeros(N,1) ;
            CF  = zeros(N,1) ;
            CS  = zeros(N,1) ;
            %muf = zeros(N,1) ;
            
            switch mode
                case 'chg'
                    Tin = obj.TC ;
                case 'dis'
                    Tin = obj.TD ;
            end
            
            muf = obj.mu .* ones(N,1) ;
            % Find properties that vary with temperature
            if obj.Lconst
                CF  = obj.cF .* ones(N,1) ;
                CS  = obj.cS .* ones(N,1) ;
                muf = obj.mu .* ones(N,1) ;
            elseif strcmp(fld.read,'CP') && ~obj.Lideal
                for i = 1 : N
                    % Temperature array from previous timestep
                    Xs = [1 Ts(i,2) Ts(i,2)*Ts(i,2)] ;
                    CS(i) = Xs * obj.CSfacs ;
                    CF(i) = CP1('PT_INPUTS',obj.P(i,2),obj.TF(:,1),'CPMASS',fld.handle) ;
                    %muf(i) = Xf * obj.MUfacs ;                    
                end
            else % This accounts for fld.read='TAB' and also ideal gases
                for i = 1 : N
                    % Temperature array from previous timestep
                    Xf = [1 Tf(i,2) Tf(i,2)*Tf(i,2)] ;
                    Xs = [1 Ts(i,2) Ts(i,2)*Ts(i,2)] ;
                    CF(i) = Xf * obj.CFfacs ;
                    CS(i) = Xs * obj.CSfacs ;
                    %muf(i) = Xf * obj.MUfacs ;                    
                end
            end
            
            % Reynolds number, Nusselt number and heat transfer coef
            for i = 1 : N
                Re  = obj.rho(i,2) .* obj.u(i,2) .* obj.dp ./ muf(i) ;
                [Nu(i), h(i)] = nusselt(Re, obj.Pr, obj) ;
            end
            
            % First calculate some coefficients at each node
            a = obj.u(:,2) .* obj.dt / (obj.eps * obj.dx) ;
            b = obj.keff .* obj.dt ./ (obj.eps .* obj.rho(:,2) .* CF * obj.dx * obj.dx) ;
            c = (1-obj.eps) .* obj.Sv .* h .* obj.dt ./ (2 .* obj.eps .* obj.rho(:,2) .* CF) ;
            d = obj.Sv .* h .* obj.dt ./ (2 .* obj.rhoS .* CS) ;
            e = obj.dx ./ (obj.rho(:,2) .* CF(:)) ;
            
            fact1 = obj.Sv * (1-obj.eps) * obj.dx / (2. * obj.eps^3) ;
            
            % Calculate the determinant at each node
            DET = ((1 + d) .* (1 + a + c) - c .* d ) ;
            
            % terms for matrix calculation
            t1 = (1 - d) .* Ts(:,2) + d .* Tf(:,2) ;
            t2 = (1 - c) .* Tf(:,2) + c .* Ts(:,2) - 2 .* b .* Tf(:,2) ;
            
            % Bit annoying Matlab doesn't do zero indices - otherwise the
            % first node could be easily wrapped into the for loop below
            % First node
            X = [(1+a(1)+c(1)), d(1) ; c(1), (1+d(1))] ; % Could probably make a matrix of matrices, urgh
            Y = [ t1(1) ; t2(1) + a(1) * Tin + b(1) * (Tf(2,2) + Tin)] ;
            T(1,:) = X * Y / DET(1) ;
            
            % Calculate density at first node
            if obj.Lconst
                obj.rho(1,1) = obj.rho(1,2)  ;
            else
                if strcmp(fld.read,'TAB')
                    Xf = [1 T(1,1) T(1,1)*T(1,1) ] ;
                    obj.rho(1,1) = Xf * obj.Dfacs ;
                elseif obj.Lideal
                    obj.rho(1,1) = obj.P(1,2) / (T(1,1) * obj.R) ;
                elseif strcmp(fld.read,'CP')
                    obj.rho(1,1) = CP1('PT_INPUTS',obj.P(1,2),obj.TF(1,1),'D',fld.handle) ;
                end
            end
                
            % Calculate velocities and pressures
            obj.u(1,1) = obj.us ;
            obj.P(1,1) = obj.Pin ;
            
            % For final node extrapolate Tf at N+1
            Tf(N+1,2) = 2.*Tf(N,2) - Tf(N-1,2) ;
            
            % Remaining nodes
            for i = 2 : N
                f = e(i) * (obj.P(i-1,1) - obj.P(i-1,2)) ;
                X = [(1 + a(i) + c(i)), d(i) ; c(i), (1 + d(i))] ; % Could probably make a matrix of matrices, urgh
                Y = [ t1(i) ; t2(i) + a(i) * Tf(i-1,1) + b(i) * (Tf(i+1,2) + Tf(i-1,2)) + f] ;
                T(i,:) = X * Y / DET(i) ;
                
                % Calculate densities
                if obj.Lconst
                    obj.rho(i,1) = obj.rho(i,2)  ;
                else
                    if strcmp(fld.read,'TAB')
                        Xf = [1 T(i,1) T(i,1)*T(i,1) ] ;
                        obj.rho(i,1) = Xf * obj.Dfacs ;
                    elseif obj.Lideal
                        obj.rho(i,1) = obj.P(i,2) / (T(i,1) * obj.R) ;
                    elseif strcmp(fld.read,'CP')
                        obj.rho(i,1) = CP1('PT_INPUTS',obj.P(i,2),obj.TF(:,1),'D',fld.handle) ;
                    end
                end
                
                % Calculate velocities and pressures
                obj.u(i,1) = (obj.rho(i-1,1) * obj.u(i-1,1) + obj.eps * obj.dx * (obj.rho(i,2) - obj.rho(i,1))/obj.dt) / obj.rho(i,1) ;
                fact2      = fact1 * obj.rho(i,1) * obj.u(i,1) * obj.u(i,1) ;
                obj.P(i,1) = obj.P(i-1,1) - fact2 * friction(obj,obj.rho(i,1),obj.u(i,1)) ;
                
            end
                       
            % Assign temperatures
            obj.TS(:,1) = T(:,1) ;
            obj.TF(:,1) = T(:,2) ;
            
            
            
            % Move new results into previous results, for use in future steps
            obj.TS(:,2)  = obj.TS(:,1) ;
            obj.TF(:,2)  = obj.TF(:,1) ;
            obj.u(:,2)   = obj.u(:,1) ;
            obj.P(:,2)   = obj.P(:,1) ;
            obj.rho(:,2) = obj.rho(:,1) ;
            
            
        end
        
        %***** ALTERNATIVE !! ******
        % March forward one timestep. Routine is intended for ideal gases
        function obj = PB_TIMESTEP_IDEAL(obj, fld, iL, iG, mode)
            
            N  = obj.NX ;
            Tf = obj.TF ;
            Ts = obj.TS ;
                        
            % Some arrays
            St  = zeros(N,1) ;
            CF  = zeros(N,1) ;
            CS  = zeros(N,1) ;
            %muf = zeros(N,1) ;
                        
            Tin = fld.state(iL,iG).T;
            pin = fld.state(iL,iG).p;
            Rin = fld.state(iL,iG).rho;
            
            muf = obj.mu .* ones(N,1) ;
            % Find properties that vary with temperature
            if obj.Lconst
                CF  = obj.cF .* ones(N,1) ;
                CS  = obj.cS .* ones(N,1) ;
                muf = obj.mu .* ones(N,1) ;
            elseif strcmp(fld.read,'CP') && ~obj.Lideal
                for i = 1 : N
                    % Temperature array from previous timestep
                    Xs = [1 Ts(i,2) Ts(i,2)*Ts(i,2)] ;
                    CS(i) = Xs * obj.CSfacs ;
                    CF(i) = CP1('PT_INPUTS',obj.P(i,2),obj.TF(:,1),'CPMASS',fld.handle) ;
                    %muf(i) = Xf * obj.MUfacs ;                    
                end
            else % This accounts for fld.read='TAB' and also ideal gases
                for i = 1 : N
                    % Temperature array from previous timestep
                    Xf = [1 Tf(i,2) Tf(i,2)*Tf(i,2)] ;
                    Xs = [1 Ts(i,2) Ts(i,2)*Ts(i,2)] ;
                    CF(i) = Xf * obj.CFfacs ;
                    CS(i) = Xs * obj.CSfacs ;
                    %muf(i) = Xf * obj.MUfacs ;                    
                end
            end
            
            % Reynolds number, Nusselt number and heat transfer coef
            for i = 1 : N
                Re      = obj.rho(i,2) .* obj.u(i,2) .* obj.dp ./ muf(i) ;
                [Nu, ~] = nusselt(Re, obj.Pr, obj) ;
                St(i)   = Nu / (Re * obj.Pr) ;
            end
            
            % First calculate some coefficients at each node
            ell = 1 ./ ((1-obj.eps) * obj.Sv * St) ;
            tau = obj.rhoS .* CS ./ (obj.rho(:,2) .* obj.u(:,2) .* CF .* St .* obj.Sv) ;
            alp = obj.keff / ((1-obj.eps) .* CS .* obj.rhoS) ;
            C   = (obj.eps .* obj.dx ./ (obj.u(:,2) * obj.dt)) ; % Unsteady term. Gets updated as go along.
            
            f1 = obj.dx ./ (2 .* ell) ;
            f2 = obj.dt ./ (2 .* tau) ;
            f3 = alp .* obj.dt ./ (obj.dx * obj.dx) ;
            f4 = obj.Sv * (1-obj.eps) * obj.dx / (2. * obj.eps^3) ;
            
            % Calculate the determinant at each node
            DET = (1. + f1) .* (1. + f2) - (f1 .* f2) ;
            
            % terms for matrix calculation
            t1 = (1 - f2 - 2*f3) .* Ts(:,2) + f2 .* Tf(:,2) ;
                       
            % Inlet temperatures of gas and solid
            TinG = Tin ;
            TinS = exp(-2. * f2(1)) ;
            TinS = TinG * (1. - TinS) + obj.TS(1,2) * TinS ; % Based on analytical solution to a simple differential equation
            
            % Bit annoying Matlab doesn't do zero indices - otherwise the
            % first node could be easily wrapped into the for loop below
            % First node - assume C(1) = 0
            X = [(1 + f2(1)), f1(1) ; f2(1), (1 + f1(1))] ; 
            Y = [(1.- f1(1))*TinG + f1(1)*TinS ; t1(1) + f3(1) * (TinS + Ts(2,2)) ] ;   
            
            T(1,:) = X * Y / DET(1) ;
            % Assign temperatures
            Ts(1,1) = T(1,2) ;
            Tf(1,1) = T(1,1) ;
            
            % Calculate pressures
            f5         = f4 * obj.rho(1,2) * obj.u(1,2) * obj.u(1,2) ;
            obj.P(1,1) = pin - f5 * friction(obj,obj.rho(1,2),obj.u(1,2)) ;
            
            % Calculate density at first node
            if obj.Lconst
                obj.rho(1,1) = obj.rho(1,2)  ;
            else
                if strcmp(fld.read,'TAB')
                    Xf = [1 T(1,1) T(1,1)*T(1,1) ] ;
                    obj.rho(1,1) = Xf * obj.Dfacs ;
                elseif obj.Lideal
                    obj.rho(1,1) = obj.P(1,1) / (T(1,1) * obj.R) ;
                elseif strcmp(fld.read,'CP')
                    obj.rho(1,1) = CP1('PT_INPUTS',obj.P(1,1),obj.TF(1,1),'D',fld.handle) ;
                end
            end
            
            % Calculate mass flux and velocity
            obj.u(1,1) = obj.mdot / (obj.A * Rin) ;
            %obj.u(1,1) = (Rin * obj.us + obj.eps * obj.dx * (obj.rho(1,2) - obj.rho(1,1))/obj.dt) / obj.rho(1,1) ;
            
            % For final node extrapolate Tf at N+1
            Ts(N+1,2) = 2.*Ts(N,2) - Ts(N-1,2) ;
            
            % Remaining nodes
            for i = 2 : N
                cc     = (obj.P(i-1,1) - obj.P(i-1,2)) / (obj.rho(i-1,1) * CF(i)) ;
                cc     = C(i) * (cc - (Tf(i-1,1) - Tf(i-1,2))) ;
                X      = [(1 + f2(i)), f1(i) ; f2(i), (1 + f1(i))] ;
                Y      = [(1-f1(i))*Tf(i-1,1) + f1(i) * Ts(i-1,1) + cc ; t1(i) + f3(i) * (Ts(i-1,2) + Ts(i+1,2)) ] ;
                T(i,:) = X * Y / DET(i) ;
                                
                % Assign temperatures
                Ts(i,1) = T(i,2) ;
                Tf(i,1) = T(i,1) ;
                
                % Calculate pressures
                f5         = f4 * obj.rho(i,2) * obj.u(i,2) * obj.u(i,2) ;
                obj.P(i,1) = obj.P(i-1,1) - f5 * friction(obj,obj.rho(i,2),obj.u(i,2)) ;
                if imag(obj.P(i,1))>0
                    keyboard
                end
                
                % Calculate densities
                if obj.Lconst
                    obj.rho(i,1) = obj.rho(i,2)  ;
                else
                    if strcmp(fld.read,'TAB')
                        Xf = [1 T(i,1) T(i,1)*T(i,1) ] ;
                        obj.rho(i,1) = Xf * obj.Dfacs ;
                    elseif obj.Lideal
                        obj.rho(i,1) = obj.P(i,1) / (T(i,1) * obj.R) ;
                    elseif strcmp(fld.read,'CP')
                        obj.rho(i,1) = CP1('PT_INPUTS',obj.P(i,1),obj.TF(:,1),'D',fld.handle) ;
                    end
                end
                
                % Calculate velocities
                obj.u(i,1) = (obj.rho(i-1,1) * obj.u(i-1,1) + obj.eps * obj.dx * (obj.rho(i,2) - obj.rho(i,1))/obj.dt) / obj.rho(i,1) ;
    
            end
            
            % Assign temperatures
            obj.TS(:,1) = T(:,2) ;
            obj.TF(:,1) = T(:,1) ;
            
            % Move new results into previous results, for use in future steps
            obj.TS(:,2)  = obj.TS(:,1) ;
            obj.TF(:,2)  = obj.TF(:,1) ;
            obj.u(:,2)   = obj.u(:,1) ;
            obj.P(:,2)   = obj.P(:,1) ;
            obj.rho(:,2) = obj.rho(:,1) ;
            
        end
        
        
        
%         %***** ALTERNATIVE V2 !! ******
%         % This follows the same methodology as PB_TIMESTEP but the routine
%         % has been updated to (hopefully) capture the unsteady energy term
%         % of liquids more accurately.
%         % Doesn't really work :(
%         % March forward one timestep. Routine is intended for ideal gases
%         function obj = PB_TIMESTEP_IDEAL2(obj, fld, mode)
%             
%             N  = obj.NX ;
%             Tf = obj.TF ;
%             Ts = obj.TS ;
%                         
%             % Some arrays
%             St  = zeros(N,1) ;
%             CF  = zeros(N,1) ;
%             CS  = zeros(N,1) ;
%             %muf = zeros(N,1) ;
%             
%             switch mode
%                 case 'chg'
%                     Tin = obj.TC ;
%                 case 'dis'
%                     Tin = obj.TD ;
%             end
%             
%             muf = obj.mu .* ones(N,1) ;
%             % Find properties that vary with temperature
%             if obj.Lconst
%                 CF  = obj.cF .* ones(N,1) ;
%                 CS  = obj.cS .* ones(N,1) ;
%                 muf = obj.mu .* ones(N,1) ;
%             elseif strcmp(fld.read,'CP') && ~obj.Lideal
%                 for i = 1 : N
%                     % Temperature array from previous timestep
%                     Xs = [1 Ts(i,2) Ts(i,2)*Ts(i,2)] ;
%                     CS(i) = Xs * obj.CSfacs ;
%                     CF(i) = CP1('PT_INPUTS',obj.P(i,2),obj.TF(:,1),'CPMASS',fld.handle) ;
%                     %muf(i) = Xf * obj.MUfacs ;                    
%                 end
%             else % This accounts for fld.read='TAB' and also ideal gases
%                 for i = 1 : N
%                     % Temperature array from previous timestep
%                     Xf = [1 Tf(i,2) Tf(i,2)*Tf(i,2)] ;
%                     Xs = [1 Ts(i,2) Ts(i,2)*Ts(i,2)] ;
%                     CF(i) = Xf * obj.CFfacs ;
%                     CS(i) = Xs * obj.CSfacs ;
%                     %muf(i) = Xf * obj.MUfacs ;                    
%                 end
%             end
%             
%             % Reynolds number, Nusselt number and heat transfer coef
%             for i = 1 : N
%                 Re      = obj.rho(i,2) .* obj.u(i,2) .* obj.dp ./ muf(i) ;
%                 [Nu, ~] = nusselt(Re, obj.Pr, obj) ;
%                 St(i)   = Nu / (Re * obj.Pr) ;
%             end
%             
%             % First calculate some coefficients at each node
%             ell = 1 ./ ((1-obj.eps) * obj.Sv * St) ;
%             tau = obj.rhoS .* CS ./ (obj.rho(:,2) .* obj.u(:,2) .* CF .* St .* obj.Sv) ;
%             alp = obj.keff / ((1-obj.eps) .* CS .* obj.rhoS) ;
%                         
%             a = obj.dx ./ (2 .* ell) ;
%             b = obj.eps * obj.dx / (2 .* obj.u(:,2) * obj.dt) ;
%             C = (obj.eps .* obj.dx ./ (obj.u(:,2) * obj.dt)) ; % Unsteady term. Gets updated as go along.
%             d = obj.dt ./ (2 .* tau) ;
%             e = alp .* obj.dt ./ (obj.dx * obj.dx) ;
%             f = obj.Sv * (1-obj.eps) * obj.dx / (2. * obj.eps^3) ;
%             
%             % Calculate the determinant at each node
%             DET = (1 + a + b).*(1 + d) - a .* d ;
%             
%             % First node - assume C(1) = 0
%             X = [(1 + d(1)), a(1) ; d(1), (1 + a(1) + b(1))] ; 
%             y1 = Tin + b(1) * Tf(1,2) ;
%             y2 = (1-d(1))*Ts(1,2) + d(1)*Tf(1,2) + e(1)*(Ts(2,2) - 2*Ts(1,2) + Tin) ; 
%             Y = [y1 ; y2 ] ;   
%             
%             T(1,:) = X * Y / DET(1) ;
%             % Assign temperatures
%             Ts(1,1) = T(1,2) ;
%             Tf(1,1) = T(1,1) ;
%             
%             % Calculate density at first node
%             if obj.Lconst
%                 obj.rho(1,1) = obj.rho(1,2)  ;
%             else
%                 if strcmp(fld.read,'TAB')
%                     Xf = [1 T(1,1) T(1,1)*T(1,1) ] ;
%                     obj.rho(1,1) = Xf * obj.Dfacs ;
%                 elseif obj.Lideal
%                     obj.rho(1,1) = obj.P(1,2) / (T(1,1) * obj.R) ;
%                 elseif strcmp(fld.read,'CP')
%                     obj.rho(1,1) = CP1('PT_INPUTS',obj.P(1,2),obj.TF(1,1),'D',fld.handle) ;
%                 end
%             end
%                 
%             % Calculate velocities and pressures
%             obj.u(1,1) = obj.us ;
%             obj.P(1,1) = obj.Pin ;
%             
%             % For final node extrapolate Tf at N+1
%             Ts(N+1,2) = 2.*Ts(N,2) - Ts(N-1,2) ;
%             
%             % Remaining nodes
%             for i = 2 : N
%                 
%                 c  = C(i) * (obj.P(i-1,1) - obj.P(i-1,2)) / (obj.rho(i-1,1) * CF(i)) ;
%                 y1 = (1-a(i)-b(i))*Tf(i-1,1) + a(i)*Ts(i-1,1) + b(i)*(Tf(i-1,2) + Tf(i,2)) +  c ;
%                 y2 = (1-d(i))*Ts(i,2) + d(i)*Tf(i,2) + e(i)*(Ts(i+1,2) - 2*Ts(i,2) + Ts(i-1,2)) ;
%                 
%                 X = [(1 + d(i)), a(i) ; d(i), (1 + a(i) + b(i))] ;
%                 Y = [y1 ; y2 ] ;
%                 T(i,:) = X * Y / DET(i) ;
%                 
%                 % Assign temperatures
%                 Ts(i,1) = T(i,2) ;
%                 Tf(i,1) = T(i,1) ;
%                 
%                 % Calculate densities
%                 if obj.Lconst
%                     obj.rho(i,1) = obj.rho(i,2)  ;
%                 else
%                     if strcmp(fld.read,'TAB')
%                         Xf = [1 T(i,1) T(i,1)*T(i,1) ] ;
%                         obj.rho(i,1) = Xf * obj.Dfacs ;
%                     elseif obj.Lideal
%                         obj.rho(i,1) = obj.P(i,2) / (T(i,1) * obj.R) ;
%                     elseif strcmp(fld.read,'CP')
%                         obj.rho(i,1) = CP1('PT_INPUTS',obj.P(i,2),obj.TF(:,1),'D',fld.handle) ;
%                     end
%                 end
%                 
%                 % Calculate velocities and pressures
%                 obj.u(i,1) = (obj.rho(i-1,1) * obj.u(i-1,1) + obj.eps * obj.dx * (obj.rho(i,2) - obj.rho(i,1))/obj.dt) / obj.rho(i,1) ;
%                 g          = f * obj.rho(i,1) * obj.u(i,1) * obj.u(i,1) ;
%                 obj.P(i,1) = obj.P(i-1,1) - g * friction(obj,obj.rho(i,1),obj.u(i,1)) ;
%                 
%             end
%             
%             % Assign temperatures
%             obj.TS(:,1) = T(:,2) ;
%             obj.TF(:,1) = T(:,1) ;
%             
%             % Move new results into previous results, for use in future steps
%             obj.TS(:,2)  = obj.TS(:,1) ;
%             obj.TF(:,2)  = obj.TF(:,1) ;
%             obj.u(:,2)   = obj.u(:,1) ;
%             obj.P(:,2)   = obj.P(:,1) ;
%             obj.rho(:,2) = obj.rho(:,1) ;
%             
%         end
        
           
        
        
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
            
            Rem = rhoF * us / ((1-obj.eps) * obj.Sv * obj.mu) ;
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
    
    
    % Other methods
    methods (Static)
        
        %%% SEE IF THIS CAN BE REMOVED (METHOD HAS BEEN CONSOLIDATED WITH
        %%% CREATE_TABLE FUNCTION)
        %%% --------------------->
        %{
        function [A] = create_solid_table(Sname)
            % THis function is based on CREATE_TABLE
            % Except it has been modified to provide data for a solid
            % material
            % The data columns are the same as a fluid (even though it
            % includes Prandtl number etc.) so that the same routines can
            % call this table and fluid tables
            
            %   Create a table with thermophysical properties as a
            %   function of temperature. Tbot and Ttop define the bottom and top
            %   temperatures within which data is correct. However, a broader range
            %   between Tmin and Tmax (with constant properties) is allowed for
            %   screening purposes (i.e. program does not crash when reading outside
            %   Tbot and Ttop).
            
            fileName  = strcat('./Data/',Sname,'.dat');
            
            % Check if file already exists. If it does, return control to invoking
            % function. Otherwise, create file and proceed.
            if isfile(fileName)
                A = load(fileName);
                return
            else
                fileID = fopen(fileName,'w');
            end
            
            % Set temperature limits
            Tmin = 50;    % Not correct. Used only for screening purposes
            Tmax = 1300;  % Not correct. Used only for screening purposes
            if strcmp(Sname,'Magnetite') % Magnetite
                Tbot  = 100;
                Ttop  = 1200;
                
            elseif strcmp(Sname,'SiO2') % Silicon Oxide
                Tbot  = 100;
                Ttop  = 1200;
                
            elseif strcmp(Sname,'Generic') % Imaginary (constant prop.) solid
                Tbot = Tmin;
                Ttop = Tmax;
                
            else
                fprintf(1,'substance == %s',Sname);
                error('Unknown substance')
            end
            
            % Set arrays
            T  = ((Tmin):1:(Tmax))'; % temperature array, K
            n  = length(T);          % get dimensions of T array
            Cp = zeros(n,1);         % specific heat capacity, J/kg/K
            h  = zeros(n,1);         % specific enthalpy, J/kg
            rho= zeros(n,1);         % density, kg/m3
            s  = zeros(n,1);         % specific entropy, J/kg/K
            k  = zeros(n,1);         % conductivity, W/(m.K)
            mu = zeros(n,1);         % viscosity, Pa.s
            
            % Obtain thermophysical properties
            bot = T<Tbot;
            top = T>Ttop;
            mid = ~bot & ~top;
            
            % Would really prefer to read this data in from a file to make
            % more user friendly
            if strcmp(Sname,'Magnetite')
                % Data from JANEF Tables and McTigue thesis
                T_dat   = [ 10,15,20,25,30,35,40,45,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,350,400,500,600,700,800,900,1000,1100,1200,1300,1400];
                Cp_dat  = [ 0.487877348,1.825022673,4.46317426,9.341947743,16.76852516,26.38151587,37.62076441,50.39592313,64.29139279,94.52171885,127.2998488,161.9391406,198.4757331,236.0422889,280.384919,312.0337033,343.6824876,370.6060894,395.5420428,419.3938242,442.1614338,463.6641762,484.263442,503.5978406,522.209458,539.5562082,556.1801771,572.0813647,587.2597711,601.8960916,615.8096307,629.0003887,641.4683654,653.3942561,706.6993738,739.106,830.904,918.007,1005.109,1092.212,1179.324,1179.324,1179.324,1179.324,1179.324,1179.324];
                rho_dat = 5175 .* ones(size(T_dat)) ;
                k_dat   = 5.0 .* ones(size(T_dat)) ;
                Pr_dat  = zeros(size(T_dat)) ;
                mu_dat  = zeros(size(T_dat)) ;
                
                Cp(mid)  = interp1(T_dat,Cp_dat ,T(mid));
                rho(mid) = interp1(T_dat,rho_dat,T(mid));
                k(mid)   = interp1(T_dat,k_dat  ,T(mid));
                mu(mid)  = 0.0 ;
                
            elseif strcmp(Sname,'SiO2')
                % Data from JANEF Tables and McTigue thesis
                T_dat   = [0, 100, 200, 298.15, 300, 400, 500, 600, 700, 800, 847];
                Cp_dat  = [0, 261.072, 543.232, 742.123, 745.119, 889.270, 992.677, 1072.134, 1144.55, 1226.653, 1273.388];
                rho_dat = 2660 .* ones(size(T_dat)) ;
                k_dat   = 5.0 .* ones(size(T_dat)) ;
                Pr_dat  = zeros(size(T_dat)) ;
                mu_dat  = zeros(size(T_dat)) ;
                
                Cp(mid)  = interp1(T_dat,Cp_dat ,T(mid));
                rho(mid) = interp1(T_dat,rho_dat,T(mid));
                k(mid)   = interp1(T_dat,k_dat  ,T(mid));
                mu(mid)  = 0.0 ;
                
            elseif strcmp(Sname,'Generic')
                Cp  = 1500.*ones(size(T));
                rho = 2000.*ones(size(T));
                k   =  0.5.*ones(size(T));
                mu  = 2e-3.*ones(size(T));
            end
            
            % Set constant properties outside of normal range
            Cp(bot)  = Cp( find(mid,1,'first')).*ones(size(T(bot)));
            rho(bot) = rho(find(mid,1,'first')).*ones(size(T(bot)));
            k(bot)   = k(  find(mid,1,'first')).*ones(size(T(bot)));
            mu(bot)  = mu( find(mid,1,'first')).*ones(size(T(bot)));
            
            Cp(top)  = Cp( find(mid,1,'last')).*ones(size(T(top)));
            rho(top) = rho(find(mid,1,'last')).*ones(size(T(top)));
            k(top)   = k(  find(mid,1,'last')).*ones(size(T(top)));
            mu(top)  = mu( find(mid,1,'last')).*ones(size(T(top)));
            
            % Compute derived parameters
            v  = 1./rho;    % specific volume, m3/kg
            Pr = Cp.*mu./k; % Prandtl number
            h(1) = 0;
            s(1) = 0;
            for i=1:(n-1)
                h(i+1)  = h(i) + 0.5*(Cp(i) + Cp(i+1))*(T(i+1)-T(i));    % enthalpy increase in isobaric process
                s(i+1)  = s(i) + (h(i+1) - h(i))/(0.5*(T(i+1) + T(i)));  % entropy  increase in isobaric process
            end
            
            %Save in a single matrix for printing to file
            A = [ T, h, v, s, Cp, k, mu, Pr];
            
            %Print to file
            fprintf(fileID,'%%%7s %20s %20s %20s %20s %20s %20s %20s\n',...
                'T[K]','h[J/kg]','v[m3/kg]','s[J/(kg.K)]','Cp[J/kg/K]','k[W/m.K]','mu[Pa.s]','Pr[-]');
            fprintf(fileID,' %7.2f %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n',A');
            fclose(fileID);
            
        end
        %}
        %%% <--------------------
        
     
    end
end