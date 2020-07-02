%**** STEP FORWARDS IN TIME ****%
n = 1; % Timestep number
time = 0.;
nout = 1 ;
while time <= Tload(Iload) && Lend == false
    
    % Constants for this time step
    ap = Ux(:,1) * Apcm * dt ./ (rhop(:,2)) ;
    as = Ar * dx;
    
    %** CALCULATE PCM PROPERTIES **%
    % Guess values of Tp and Ts for this timestep
    Ts(1,1) = Ts(1,2) ;
    Tp(1,1) = Tp(1,2) ;
    for i = 2 : Nx
        Ts(i,1) = 0.5 * (Ts(i,2) + Ts(i-1,2));
        Tp(i,1) = 0.5 * (Tp(i,2) + Tp(i-1,2));
    end
    Tpguess = Tp(:,1) ;
    
    % Calculate PCM properties along pipe
    for i = 1 : Nx
        
        % Can iterate a couple of times to improve convergence
        for jj = 1 : PCMiterations
            hp(i,1) = hp(i,2) + ap(i) *(Ts(i,1) - Tp(i,1)) ;
            
            % Calculate PCM properties
            if hp(i,1) > hpl % melted
                xp(i,1)   = 1. ;
                Tp(i,1)   = Tpm + (hp(i,1) - hpl) / cpl ;
                rhop(i,1) = rhopl;
            elseif hp(i,1) < hps % solid
                xp(i,1)   = 0. ;
                Tp(i,1)   = hp(i,1) / cps ;
                rhop(i,1) = rhops ;
            else % partially melted
                xp(i,1) = (hp(i,1) - hps) / (hpl - hps) ; % Estimate dryness fraction from lever rule
                Tp(i,1) = Tpm ;
                vp = xp(i,1) * (1./rhopl - 1./rhops) + 1./rhops ; % Specific volume from lever rule
                rhop(i,1) = 1./vp ;
            end
            
            TPerr(i) = 100 * (Tp(i,1) - Tpguess(i)) / Tpguess(i) ;
            Tpguess(i) = Tp(i,1) ;
        end
    end
    
    %%% STEAM CALCULATIONS %%%
    % Heat transfer conditions at pipe inlet
    if xs(1,1) >= 0.99
        Re(1,1) = Gs(1,1) * dp / musv ;
        if Re(1,1) < 3000
            hc(1,1) = 4.36 * ktv / dp ;
        else
            hc(1,1) = 0.023 * Re(1,1)^0.8 * PrV^0.4 * ktv / dp ;
        end
        fac(1)  = 1 ;
        hc(1,1) = hc(1,1) * fac(1) ;
        Nu(1,1) = hc(1,1) * dp / ktv ;
        St(1,1) = Nu(1,1) / (Re(1,1) * PrV) ;
    else
        Re(1,1) = Gs(1,1) * dp / musl ; % Reynolds number based on all liquid flow
        if Re(1,1) < 3000
            hc(1,1) = 4.36 * ktl / dp ;
        else
            hc(1,1) = 0.023 * Re(1,1)^0.8 * PrL^0.4 * ktl / dp ;
        end
        
        %fac(i)  = hx_mult(xs(i,1),Ps(i,1),Psc,mode) ;
        if strcmp(Load(Iload),"c") % This is temporary move to a function
            fac(1)  = (1-xs(1,1))^0.8 + 3.8 * xs(1,1)^0.76 * (1-xs(1,1))^0.04 / (Ps(1,1) / Psc)^0.38 ;
        else
            Fr = (Gs(1,1) * vl )^2 / (9.81 * dp) ; % Froud number if all flow is liquid
            Fr = max(Fr, (25.*Fr)^0.3) ;
            q  = hc(1,2) * abs(Ts(1,2) - Tp(1,2)) ; % Heat flux
            Bo = q / (Gs(1,1) * (hv - hl)); % Boiling number
            Co = ((1-xs(1,1))/xs(1,1))^0.8 * (vl / vv)^0.5 ; % Convection number
            
            fconv = 1.136 * Fr * Co^-0.9 + 667.2 * Bo^0.7 ;
            fnucl = 0.6683 * Fr * Co^-0.2 + 1058 * Bo^0.7 ;
            fac(1) = max(fconv,fnucl) ;
            %fac(i) = 50.;%
        end
        
        hc(1,1) = hc(1,1) * fac(1) ;
        Nu(1,1) = hc(1,1) * dp / ktl ;
        St(1,1) = Nu(1,1) / (Re(1,1) * PrL) ;
    end
    
    Ux(1,1)   = 1. / (1./hc(1,1) + Rpipe) ;
    
    % Step through remaining nodes to calculate steam properties
    for i = 2 : Nx
        
        Gs(i,1) = Gs(i,2) ; % Steady-state
        Ps(i,1) = Ps(i-1,1) ; % Pressure
        % Guess current value of Ts
        Ts(i,1) = 0.5 * (Ts(i,2) + Ts(i-1,2)) ;
        xs(i,1) = xs(i,2) ; % Guess
        Tguess  = Ts(i,1) ;
        Xguess  = xs(i,1) ;
        
        % Now iterate on the equation until the new Ts value remains constant
        err = 1e6 ;
        tol = 0.01 ;
        cnt = 1 ;
        cntM = 100 ;
        
        while err > tol && cnt < cntM
            
            
            if xs(i,1) >=0.99
                Re(i,1) = Gs(i,1) * dp / musv ;
                if Re(i,1) < 3000
                    hc(i,1) = 4.36 * ktv / dp ;
                else
                    hc(i,1) = 0.023 * Re(i,1)^0.8 * PrV^0.4 * ktv / dp ;
                end
                fac(i)  = 1 ;
                hc(i,1) = hc(i,1) * fac(i) ;
                Nu(i,1) = hc(i,1) * dp / ktv ;
                St(i,1) = Nu(i,1) / (Re(i,1) * PrV) ;
            else
                Re(i,1) = Gs(i,1) * dp / musl ; % Reynolds number based on all liquid flow
                if Re(i,1) < 3000
                    hc(i,1) = 4.36 * ktl / dp ;
                else
                    hc(i,1) = 0.023 * Re(i,1)^0.8 * PrL^0.4 * ktl / dp ;
                end
                %fac(i)  = hx_mult(xs(i,1),Ps(i,1),Psc,mode) ;
                if strcmp(Load(Iload),"c") % This is temporary move to a function
                    fac(i)  = (1-xs(i,1))^0.8 + 3.8 * xs(i,1)^0.76 * (1-xs(i,1))^0.04 / (Ps(i,1) / Psc)^0.38 ;
                else
                    Fr = (Gs(i,1) * vl )^2 / (9.81 * dp) ; % Froud number if all flow is liquid
                    Fr = max(Fr, (25.*Fr)^0.3) ;
                    q  = hc(i,2) * abs(Ts(i,2) - Tp(i,2)) ; % Heat flux
                    Bo = q / (Gs(i,1) * (hv - hl)); % Boiling number
                    Co = ((1-xs(i,1))/xs(i,1))^0.8 * (vl / vv)^0.5 ; % Convection number
                    
                    fconv = 1.136 * Fr * Co^-0.9 + 667.2 * Bo^0.7 ;
                    fnucl = 0.6683 * Fr * Co^-0.2 + 1058 * Bo^0.7 ;
                    fac(i) = max(fconv,fnucl) ;
                    %fac(i) = 50.;%
                end
                hc(i,1) = hc(i,1) * fac(i) ;
                Nu(i,1) = hc(i,1) * dp / ktl ;
                St(i,1) = Nu(i,1) / (Re(i,1) * PrL) ;
            end
            
            Ux(i,1)   = 1. / (1./hc(i,1) + Rpipe) ;
            
            % New enthalpy
            hs(i,1) = hs(i-1,1) + (Tp(i,1) - Ts(i,1)) * as * Ux(i,1) / Gs(i,1) ;
            
            % Calculate steam properties
            if hs(i,1) > hv % steam
                xs(i,1)   = 1. ;
                %Ts(i,1)   = RP1('HmassP_INPUTS',hs(i,1),Ps(i,1),'T',steam);
                Ts(i,1)   = Tsat + (hs(i,1)-hv) / csv ;
                %rhos(i,1) = RP1('HmassP_INPUTS',hs(i,1),Ps(i,1),'D',steam);% rhosv;
                
                %Ts(i,1) = Tcoef(1) + Tcoef(2) * hs(i,1) ;
                rhos(i,1) = 1. / (Vcoef(1) + Vcoef(2) * hs(i,1)) ;
                
            elseif hs(i,1) < hl % water
                xs(i,1)   = 0. ;
                Ts(i,1)   = hs(i,1) / csl ;
                rhos(i,1) = rhosl ;
            else % steam-water mixture
                xs(i,1) = (hs(i,1) - hl) / (hv - hl) ; % Estimate dryness fraction from lever rule
                Ts(i,1) = Tsat ;
                vs = xs(i,1) * (vv - vl) + vl ; % Specific volume from lever rule
                rhos(i,1) = 1./vs ;
            end
            
            cnt = cnt + 1 ;
            %err = 100 * abs(Ts(i,1) - Tguess) / Tguess ;
            %Tguess = Ts(i,1) ;
            err = 100 * abs(xs(i,1) - Xguess) / Xguess ;
            Xguess = xs(i,1) ;
            
        end
        
        
    end
    
    rhos(:,2) = rhos(:,1) ;
    Gs(:,2)   = Gs(:,1) ;
    hs(:,2)   = hs(:,1) ;
    xs(:,2)   = xs(:,1) ;
    Ts(:,2)   = Ts(:,1) ;
    
    rhop(:,2) = rhop(:,1) ;
    hp(:,2)   = hp(:,1) ;
    xp(:,2)   = xp(:,1) ;
    Tp(:,2)   = Tp(:,1) ;
    cp(:,2)   = cp(:,1) ;
    
    Ux(:,2)   = Ux(:,1) ;
    
    % Save data points every so often
    if mod(n,Nt/(Nprof+1))==0
        fprintf("TIMESTEP #%i COMPLETED.\n",n);
        GSsave(:,nout,Iload)   = Gs(:,1) ;
        HSsave(:,nout,Iload)   = hs(:,1) ;
        TSsave(:,nout,Iload)   = Ts(:,1) ;
        XSsave(:,nout,Iload)   = xs(:,1) ;
        RHOSsave(:,nout,Iload) = rhos(:,1) ;
        SSsave(:,nout,Iload)   = ss(:,1) ;
        
        TPsave(:,nout,Iload)   = Tp(:,1) ;
        HPsave(:,nout,Iload)   = hp(:,1) ;
        XPsave(:,nout,Iload)   = xp(:,1) ;
        RHOPsave(:,nout,Iload) = rhop(:,1) ;
        
        TERRsave(:,nout,Iload) = TPerr(:) ;
        
        UXsave(:,nout,Iload) = Ux(:,1) ;
        
        nout = nout + 1 ;
    end
    
    % Check whether the cycle should be ended
    if strcmp(Load(Iload),'c')
        if xp(end) > XPend_chg
            Lend = true ;
            Tload(Iload) = time ;
        end
    elseif strcmp(Load(Iload),'d')
        if xp(end) < XPend_dis
            Lend = true ;
            Tload(Iload) = time ;
        end
    end
    
    n    = n + 1 ;
    time = time + dt ;
end

% If the next load cycle is opposite to this one, then reverse each of
% the arrays
if ~strcmp(Load(Iload),Load(Iload+1))
    Gs   = flip(Gs,1);
    hs   = flip(hs,1);
    Ts   = flip(Ts,1);
    xs   = flip(xs,1);
    rhos = flip(rhos,1);
    ss   = flip(ss,1);
    
    Tp   = flip(Tp,1);
    hp   = flip(hp,1);
    xp   = flip(xp,1);
    rhop = flip(rhop,1);
    
    Ux = flip(Ux,1);
    
end

