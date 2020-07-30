%**** STEP FORWARDS IN TIME ****%
if Itank(Iload) == -1
    Nn = Ntank ;
else
    Nn  = 1 ; 
end

for i = 1 : Nn
    n = 1; % Timestep number
    time = 0.;
    nout = 1 ;
    if Itank(Iload) == -1
        Itk = i ;
    else
        Itk = Itank(Iload) ;
    end
    
    % Final profile of previous load step
    if Iload > 1
        if ~strcmp(Load(Iload),Load(Iload-1))
            Ux(:,2)   = flip(UXlast(:,Itk),1) ;
            rhop(:,2) = flip(RHOPlast(:,Itk),1) ;
            hp(:,2)   = flip(HPlast(:,Itk),1) ;
            Tp(:,2)   = flip(TPlast(:,Itk),1) ;
            xp(:,2)   = flip(XPlast(:,Itk),1) ;
        else
            Ux(:,2)   = UXlast(:,Itk) ;
            rhop(:,2) = RHOPlast(:,Itk) ;
            hp(:,2)   = HPlast(:,Itk) ;
            Tp(:,2)   = TPlast(:,Itk) ;
            xp(:,2)   = XPlast(:,Itk) ;
        end
    end
    
    while time <= Tload(Iload) && Lend == false
        
        % Constants for this time step
        ap = Ux(:,2) * Apcm * dt ./ (rhop(:,2)) ;
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
            hc(1,1) = hx_coef(Re(1,1), ktv, dp, PrV) ;
            
            fac(1)  = 1 ;
            hc(1,1) = hc(1,1) * fac(1) ;
            Nu(1,1) = hc(1,1) * dp / ktv ;
            St(1,1) = Nu(1,1) / (Re(1,1) * PrV) ;
        else
            Re(1,1) = Gs(1,1) * dp / musl ; % Reynolds number based on all liquid flow
            hc(1,1) = hx_coef(Re(1,1), ktl, dp, PrL) ;
            
            %fac(i)  = hx_mult(xs(i,1),Ps(i,1),Psc,mode) ;
            if strcmp(Load(Iload),"c") % This is temporary move to a function
                fac(1)  = hx_cond(xs(1,1),Ps(1,1), Psc) ;
            else
                fac(1)  = hx_boil(Gs(1,1),xs(1,1),hc(1,1),Ts(1,1),Tp(1,1),dp,hv,hl,vv,vl) ;
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
            tol = 1e-4 ;
            cnt = 1 ;
            cntM = 100 ;
            
            while err > tol && cnt < cntM
                
                
                if xs(i,1) >=0.99
                    Re(i,1) = Gs(i,1) * dp / musv ;
                    hc(i,1) = hx_coef(Re(i,1), ktv, dp, PrV) ;
                    
                    fac(i)  = 1 ;
                    hc(i,1) = hc(i,1) * fac(i) ;
                    Nu(i,1) = hc(i,1) * dp / ktv ;
                    St(i,1) = Nu(i,1) / (Re(i,1) * PrV) ;
                    
                else
                    Re(i,1) = Gs(i,1) * dp / musl ; % Reynolds number based on all liquid flow
                    hc(i,1) = hx_coef(Re(i,1), ktl, dp, PrL) ;
                    
                    %fac(i)  = hx_mult(xs(i,1),Ps(i,1),Psc,mode) ;
                    if strcmp(Load(Iload),"c") % This is temporary move to a function
                        fac(i)  = hx_cond(xs(i,1),Ps(i,1), Psc) ;
                    else
                        fac(i)  = hx_boil(Gs(i,1),xs(i,1),hc(i,1),Ts(i,1),Tp(i,1),dp,hv,hl,vv,vl) ;
                    end
                    hc(i,1) = hc(i,1) * fac(i) ;
                    Nu(i,1) = hc(i,1) * dp / ktl ;
                    St(i,1) = Nu(i,1) / (Re(i,1) * PrL) ;
                end
                
                Ux(i,1)   = 1. / (1./hc(i,1) + Rpipe) ;
                
                % New enthalpy
                %hs(i,1) = hs(i-1,1) + (Tp(i,1) - Ts(i,1)) * as * Ux(i,1) / Gs(i,1) ;
                hs(i,1) = hs(i-1,1) + 0.5*(Tp(i,1)+Tp(i-1,1) - Ts(i,1)- Ts(i-1,1)) * as * 0.5*(Ux(i,1)+Ux(i-1,1)) / Gs(i,1) ;
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
                    Ts(i,1)   = Tsat - (hl - hs(i,1))/csl ; %hs(i,1) / csl ;
                    rhos(i,1) = rhosl ;
                else % steam-water mixture
                    xs(i,1) = (hs(i,1) - hl) / (hv - hl) ; % Estimate dryness fraction from lever rule
                    Ts(i,1) = Tsat ;
                    vs = xs(i,1) * (vv - vl) + vl ; % Specific volume from lever rule
                    rhos(i,1) = 1./vs ;
                end
                
                %err = 100 * abs(Ts(i,1) - Tguess) / Tguess ;
                %Tguess = Ts(i,1) ;
                err = abs(xs(i,1) - Xguess) ;
                Xguess  = (1.-cnt/cntM) * xs(i,1) + (cnt/cntM) * Xguess;
                xs(i,1) = Xguess ;
                cnt = cnt + 1 ;
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
        
        % Check whether the cycle should be ended. Also calculate the mass
        % flow rate to give the required power
        if strcmp(Load(Iload),'c')
            if mean(xp(:,1)) > XPend_chg
                Lend = true ;
                Tload(Iload) = time ;
            end
            mdotI = abs(tes_pow(Iload) / (hs(end,1) - hs(1,1))) / Npipe ;
        elseif strcmp(Load(Iload),'d')
            if mean(xp(:,1)) < XPend_dis
                Lend = true ;
                Tload(Iload) = time ;
            end
            mdotI = abs(tes_pow(Iload) / (hs(end,1) - hs(1,1))) * Ntank / Npipe ;
        end
        Gs    = mdotI * ones(Nx,2) / Ap ; % Mass flux
        us    = Gs ./ rhos ;
        
        n    = n + 1 ;
        time = time + dt ;
        
        % Save data points every so often
        if mod(n,Nt/(Nprof+1))==0 || time > Tload(Iload) || Lend
            
            fprintf("TIMESTEP #%i COMPLETED FOR TANK #%i.\n",n,Itk);
            GSsave(:,nout,Iload,Itk)   = Gs(:,1) ;
            HSsave(:,nout,Iload,Itk)   = hs(:,1) ;
            TSsave(:,nout,Iload,Itk)   = Ts(:,1) ;
            XSsave(:,nout,Iload,Itk)   = xs(:,1) ;
            RHOSsave(:,nout,Iload,Itk) = rhos(:,1) ;
            SSsave(:,nout,Iload,Itk)   = ss(:,1) ;
            
            TPsave(:,nout,Iload,Itk)   = Tp(:,1) ;
            HPsave(:,nout,Iload,Itk)   = hp(:,1) ;
            XPsave(:,nout,Iload,Itk)   = xp(:,1) ;
            RHOPsave(:,nout,Iload,Itk) = rhop(:,1) ;
            
            % These are the last recorded values
            TPlast(:,Itk) = Tp(:,1) ;
            HPlast(:,Itk) = hp(:,1) ;
            XPlast(:,Itk) = xp(:,1) ;
            RHOPlast(:,Itk) = rhop(:,1) ;
            UXlast(:,Itk) = Ux(:,1) ;
            
            TERRsave(:,nout,Iload,Itk) = TPerr(:) ;
            
            UXsave(:,nout,Iload,Itk) = Ux(:,1) ;
            
            MDOTsave(nout,Iload,Itk) = mdotI ;
            
            HINsave(nout,Iload,Itk)  = Npipe * mdotI * hs(1,1) / Ntank ;
            HOUTsave(nout,Iload,Itk) = Npipe * mdotI * hs(end,1) / Ntank ;
            
            nout = nout + 1 ;
        end
    end
    
    % If the next load cycle is opposite to this one, then reverse each of
    % the arrays
    if Iload < numel(Load)
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
    end

end

% Calculate the heat transfer coefficient
function h = hx_coef(Re, k, d, Pr)

if Re < 3000
    h = 4.36 * k / d ;
else
    h = 0.023 * Re^0.8 * Pr^0.4 * k / d ;
end

end

% Calculate the heat transfer coefficient multiplying factor for condensing flow
function fac = hx_cond(x, Ps, Psc)

fac  = (1-x)^0.8 + 3.8 * x^0.76 * (1-x)^0.04 / (Ps / Psc)^0.38 ;

end

% Calculate the heat transfer coefficient multiplying factor for boiling flow
function fac = hx_boil(Gs,x,hc,Ts,Tp,dp,hv,hl,vv,vl)

Fr = (Gs * vl )^2 / (9.81 * dp) ; % Froud number if all flow is liquid
Fr = max(1, (25.*Fr)^0.3) ;
q  = hc * abs(Ts - Tp) ; % Heat flux
Bo = q / (Gs * (hv - hl)); % Boiling number
Co = ((1-x)/x)^0.8 * (vl / vv)^0.5 ; % Convection number

fconv = 1.136 * Fr * Co^-0.9 + 667.2 * Bo^0.7 ;
fnucl = 0.6683 * Fr * Co^-0.2 + 1058 * Bo^0.7 ;
fac = max(fconv,fnucl) * (1-x)^0.8 ;

end