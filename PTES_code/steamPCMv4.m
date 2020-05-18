%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCM with steam
% This code models a Phase Change Material (PCM) which surrounds a pipe 
% through which evaporating or condensing steam passes. The code uses 
% routines developed by Pau Farres-Antunez (Cambridge University) and
% Josh McTigue (NREL). The method is based partly on the work of Sharan et
% al. (2018) and partly on McTigue's PhD thesis.
%
% In v2 we no longer assume that both steam and PCM are changing phase at
% the same time. They may have variable temperatures.
%
% In v3, unsteady steam terms seem to be small - try a new steam routine
%
% In v4, we make further simplifications. Assume a steady-state solution
% for each instance in time. Step forward only the PCM equation in time,
% and iterate the steam profile based on this. Also include proper
% calculation of the heat transfer coefficients.
%
% Author: Josh McTigue, Pau Farres-Antunez
% 13 May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clear;

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/','./Other/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\','.\Other\')
end
% Set properties for plots
set_graphics
% Load CoolProp library (low-level interface)
load_coolprop

%**** INPUT PARAMETERS ****%
% PIPE
dp    = 0.02 ; % Pipe diameter
Lp    = 1.0 ; % Pipe length
U     = 110 ; % Overall heat transfer coefficient, W/m2K

% STEAM
steam = fluid_class('Water','WF','CP','HEOS',2,30);
Psat  = 20e5 ; % Saturation pressure
xs0   = 1.0 ;  % Initial dryness fraction of steam
Tsat  = RP1('PQ_INPUTS',Psat,1.0,'T',steam) ; % Saturation temperature
Ts0   = Tsat ;
Ps0   = Psat ;

% PCM
Tpm   = Tsat - 50 ; % Melting point of PCM
Tp0   = Tsat - 51 ; % Initial temperature of PCM
xp0   = 0.0 ;
cpl   = 2.5e3 ;     % Specific heat capacity of liquid PCM
cps   = 1.5e3 ;     % Specific heat capacity of solid PCM
rhopl = 1e3 ;       % Density of liquid PCM
rhops = 0.8e3 ;     % Density of solid PCM
Lpcm  = 2e6 ;       % Latent heat of PCM
kpcm  = 5 ;         % Thermal conductivity of PCM
cp0   = cps;
rhop0 = rhops ;

dTm   = 1.00 ;           % Melting temperature range of PCM
cpm   = Lpcm / (2.*dTm) ; % specific heat capacity of PCM while melting/freezing
rhopm = 0.5 * (rhopl + rhops) ;

% GRID AND TIME STEPS
CFL   = 1 ; % Courant-Freidrich-Lewy number
Nx    = 100 ; % Number of gridsteps
tN    = 4*3600;  % Duration of simultion, s 

%**** SET UP SOME PARAMETERS ****%

% Initial properties of steam 
hv = RP1('PQ_INPUTS',Psat,1.0,'H',steam) ; % Enthalpy of vapour
hl = RP1('PQ_INPUTS',Psat,0.0,'H',steam) ; % Enthalpy of liquid

vv = 1./RP1('PQ_INPUTS',Psat,1.0,'D',steam) ; % Specific volume of vapour
vl = 1./RP1('PQ_INPUTS',Psat,0.0,'D',steam) ; % Specific volume of liquid

hs0   = RP1('PQ_INPUTS',Ps0,xs0,'H',steam) ;
ss0   = RP1('PQ_INPUTS',Ps0,xs0,'S',steam) ; 
rhos0 = RP1('PQ_INPUTS',Ps0,xs0,'D',steam) ; 
mus0  = RP1('PQ_INPUTS',Ps0,xs0,'VISCOSITY',steam) ; 
csl   = RP1('PQ_INPUTS',Ps0,0.0,'CPMASS',steam) ; 
rhosl = RP1('PQ_INPUTS',Ps0,0.0,'D',steam) ; 
musl  = RP1('PQ_INPUTS',Ps0,0.0,'VISCOSITY',steam) ; 
ktl   = RP1('PQ_INPUTS',Ps0,0.0,'L',steam) ;  % Thermal conductivity of liquid
Psc   = RP1(0,0,0,'PCRIT',steam);      % Critical pressure of steam

csv   = RP1('PQ_INPUTS',Ps0,1.0,'CPMASS',steam) ; 
rhosv = RP1('PQ_INPUTS',Ps0,1.0,'D',steam) ; 
musv  = RP1('PQ_INPUTS',Ps0,1.0,'VISCOSITY',steam) ; 

mdot0 = U * pi * dp * Lp * (Tsat - Tpm) / (hv-hl) ;

Apcm = mdot0 * tN * (hv - hl) / (rhops * Lp * Lpcm) ;
do   = sqrt(4.*Apcm/pi + dp^2) ;

% Areas
Ap    = 0.25 * pi * dp^2 ; % Pipe area
Ar    = 4.0 / dp ; % Heat transfer area
Apcm  = 4.0 / (dp * ((do/dp)^2 - 1.)) ;

% Thermal resistance of pipe
Rpipe = 0.5 * dp * log (do/dp) / kpcm ;

Gs0   = mdot0 / Ap ; % Mass flux per unit area
us0   = Gs0 / rhos0 ;

% GRID
dx   = Lp / Nx ; % Grid step size
dt   = CFL * dx / us0 ; % Time step size
Nt   = int64(tN / dt) ; % Number of timesteps

time = 0.;
nout = 1 ;

% Set up some property arrays. 1st column is new value, 2nd column is old.
% steam properties
rhos  = rhos0 * ones(Nx,2) ; % 's' stands for steam
Gs    = Gs0   * ones(Nx,2) ; % Mass flux
Ts    = Ts0   * ones(Nx,2) ; % Temperature
Ps    = Ps0   * ones(Nx,2) ; % Pressure
hs    = hs0   * ones(Nx,2) ; % Enthalpy
ss    = ss0   * ones(Nx,2) ; % Entropy
xs    = xs0   * ones(Nx,2) ; % Dryness fraction 
mus   = mus0  * ones(Nx,2) ; % Viscosity

% Some dimensionless numbers
Re0   = Gs0 * dp / musl ; % Reynolds number based on all liquid flow
Pr0   = RP1('PQ_INPUTS',Ps0,0.0,'PRANDTL',steam) ; % Prandtl number of liquid
Nu0   = 0.023 * Re0^0.8 * Pr0^0.4 ; % Dittus-Boelter
hc0   = Nu0 * ktl / dp ;
St0   = Nu0 / (Re0 * Pr0) ;
U0    = 1. / (1./hc0 + Rpipe) ;

Re    = Re0 * ones(Nx,2) ; % Reynolds number
Pr    = Pr0 * ones(Nx,2) ; % Prandtl number
Nu    = Nu0 * ones(Nx,2) ; % Nusselt number
St    = St0 * ones(Nx,2) ; % Stanton number
hc    = hc0 * ones(Nx,2) ; % Steam heat transfer coefficient

% pcm properties
rhop  = rhop0 * ones(Nx,2) ; % 'p' stands for PCM
Tp    = Tp0   * ones(Nx,2) ; 
%hp    = hp0   * ones(Nx,2) ; 
%sp    = sp0   * ones(Nx,2) ;
xp    = xp0   * ones(Nx,2) ;
cp    = cp0   * ones(Nx,2) ; 

Ux    = U0   * ones(Nx,2) ; 

tic



%**** STEP FORWARDS IN TIME ****%
for n = 1 : Nt
    
    % constants
    ap = U * Apcm * dt ./ (rhop(:,2) .* cp(:,2)) ;
    %ap = Ux(:,1) * Apcm * dt ./ (rhop(:,2) .* cp(:,2)) ;
    as = U * Ar * dx;
    %as = Ar * dx;
    
    % Calculate PCM properties along pipe
    for i = 1 : Nx
       Tp(i,1) = (Tp(i,2) + ap(i) .* Ts(i,2)) ./ (1. + ap(i)) ; % PCM temperature
       
       % New heat capacity, density, and melted fraction
       if Tp(i,1) < Tpm - dTm
           cp(i,1)   = cps ;
           rhop(i,1) = rhops ;
           xp(i,1)   = 0.0 ;
       elseif Tp(i,1) > Tpm + dTm
           cp(i,1)   = cpl ;
           rhop(i,1) = rhopl ;
           xp(i,1)   = 1.0 ;
       else
           cp(i,1)   = cpm ;
           rhop(i,1) = rhops ;
           xp(i,1)   = (Tp(i,1) - Tpm + dTm) / (2. * dTm) ;
       end
        
    end
       
    % Steam conditions at the pipe inlet
    rhos(1,1) = rhos0 ;
    Gs(1,1)   = Gs0 ;
    hs(1,1)   = hs0 ;
    Ts(1,1)   = Ts0 ;
    xs(1,1)   = xs0 ;
    Ps(1,1)   = Ps0 ;
    us(1,1)   = Gs(1,1) / rhos(1,1) ;
    
    % Heat transfer conditions at pipe inlet
    Re(1,1) = Gs(1,1) * dp / musl ; % Reynolds number based on all liquid flow
    hc(1,1)   = 0.023 * Re(1,1)^0.8 * Pr0^0.4 * ktl / dp ;
    fac       = (1-xs(1,1))^0.8 + 3.8 * xs(1,1)^0.76 * (1-xs(1,1))^0.04 / (Ps(1,1) / Psc)^0.38 ;
    if xs(1,1) >=1
        fac = 8 ;
    end
    hc(1,1)   = hc(1,1) * fac ;
    Nu(1,1)   = hc(1,1) * dp / ktl ;
    St(1,1)   = Nu(1,1) / (Re(1,1) * Pr0) ;
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
            
            % Heat transfer conditions at pipe inlet
            Re(i,1)   = Gs(i,1) * dp / musl ; % Reynolds number based on all liquid flow
            hc(i,1)   = 0.023 * Re(i,1)^0.8 * Pr0^0.4 * ktl / dp ;
            fac       = (1-xs(i,1))^0.8 + 3.8 * xs(i,1)^0.76 * (1-xs(i,1))^0.04 / (Ps(i,1) / Psc)^0.38 ;
            if xs(i,1) >=1
                fac = 8 ;
            end
            hc(i,1)   = hc(i,1) * fac ;
            Ux(i,1)   = 1. / (1./hc(i,1) + Rpipe) ;
            Nu(i,1)   = hc(i,1) * dp / ktl ;
            St(i,1)   = Nu(i,1) / (Re(i,1) * Pr0) ;
    
            % New enthalpy
            hs(i,1) = hs(i-1,1) + (Tp(i,1) - Ts(i,1)) * as / Gs(i,1) ;
            %hs(i,1) = hs(i-1,1) + (Tp(i,1) - Ts(i,1)) * as * Ux(i,1) / Gs(i,1) ;
            
            % Calculate steam properties
            if hs(i,1) > hv % steam
                xs(i,1)   = 1. ;
                Ts(i,1)   = RP1('HmassP_INPUTS',hs(i,1),Ps0,'T',steam);%Tsat + (hs(i,1)-hv) / csv ;
                rhos(i,1) = RP1('HmassP_INPUTS',hs(i,1),Ps0,'D',steam);% rhosv;
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
%    hp(:,2)   = hp(:,1) ;
    xp(:,2)   = xp(:,1) ;
    Tp(:,2)   = Tp(:,1) ;
    cp(:,2)   = cp(:,1) ;
    
    Ux(:,2)   = Ux(:,1) ;
    
    % Save data points every so often
    if mod(n,Nt/10)==0
        fprintf("TIMESTEP #%i COMPLETED.\n",n);
        GSsave(:,nout) = Gs(:,1) ;
        HSsave(:,nout) = hs(:,1) ;
        TSsave(:,nout) = Ts(:,1) ;
        XSsave(:,nout) = xs(:,1) ;
        RHOSsave(:,nout) = rhos(:,1) ;
        SSsave(:,nout) = ss(:,1) ;
        
        TPsave(:,nout) = Tp(:,1) ;
        XPsave(:,nout) = xp(:,1) ;
        RHOPsave(:,nout) = rhop(:,1) ;
        
        UXsave(:,nout) = Ux(:,1) ;
        
        nout = nout + 1 ;
    end
    
    time = time + dt ;
end

toc