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
% Author: Josh McTigue, Pau Farres-Antunez
% 17 April 2020
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
dp    = 0.1 ; % Pipe diameter
do    = 0.5 ; % Outer pipe diameter
Lp    = 1.0 ; % Pipe length
U     = 1e2 ; % Overall heat transfer coefficient, W/m2K

% STEAM
steam = fluid_class('Water','WF','CP','HEOS',2,30);
Psat  = 20e5 ; % Saturation pressure
xs0   = 1.0 ;  % Initial dryness fraction of steam
Tsat  = RP1('PQ_INPUTS',Psat,1.0,'T',steam) ; % Saturation temperature
Ts0   = Tsat ;
Ps0   = Psat ;
mdot0 = 5.0;  % Mass flow rate, kg/s

% PCM
Tpm   = Tsat - 15 ; % Melting point of PCM
Tp0   = Tsat - 16 ; % Initial temperature of PCM
xp0   = 0.0 ;
cpl   = 2.5e3 ;     % Specific heat capacity of liquid PCM
cps   = 1.5e3 ;     % Specific heat capacity of solid PCM
rhopl = 1e3 ;       % Density of liquid PCM
rhops = 0.8e3 ;     % Density of solid PCM
Lpcm  = 2e6 ;       % Latent heat of PCM
cp0   = cps;
rhop0 = rhops ;

dTm   = 1.00 ;           % Melting temperature range of PCM
cpm   = Lpcm / (2.*dTm) ; % specific heat capacity of PCM while melting/freezing
rhopm = 0.5 * (rhopl + rhops) ;

% GRID AND TIME STEPS
CFL   = 100. ; % Courant-Freidrich-Lewy number
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
csl   = RP1('PQ_INPUTS',Ps0,0.0,'CPMASS',steam) ; 
rhosl = RP1('PQ_INPUTS',Ps0,0.0,'D',steam) ; 
csv   = RP1('PQ_INPUTS',Ps0,1.0,'CPMASS',steam) ; 
rhosv = RP1('PQ_INPUTS',Ps0,1.0,'D',steam) ; 

mdot0 = 0.3*U * pi * dp * Lp * (Tsat - Tpm) / (hv-hl) ;

Apcm = mdot0 * tN * (hv - hl) / (rhops * Lp * Lpcm) ;
do   = sqrt(4.*Apcm/pi + dp^2) ;

% Areas
Ap    = 0.25 * pi * dp^2 ; % Pipe area
Ar    = 4.0 / dp ; % Heat transfer area
Apcm  = 4.0 / (dp * ((do/dp)^2 - 1.)) ;

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
Gs    = Gs0   * ones(Nx,2) ;
Ts    = Ts0   * ones(Nx,2) ;
Ps    = Ps0   * ones(Nx,2) ;
hs    = hs0   * ones(Nx,2) ;
ss    = ss0   * ones(Nx,2) ;
xs    = xs0   * ones(Nx,2) ;

% pcm properties
rhop  = rhop0 * ones(Nx,2) ; % 'p' stands for PCM
Tp    = Tp0   * ones(Nx,2) ; 
%hp    = hp0   * ones(Nx,2) ; 
%sp    = sp0   * ones(Nx,2) ;
xp    = xp0   * ones(Nx,2) ;
cp    = cp0   * ones(Nx,2) ; 
tic
%**** STEP FORWARDS IN TIME ****%
for n = 1 : Nt
    
    % constants
    ap = U * Apcm * dt ./ (rhop(:,2) .* cp(:,2)) ;
    as = U * Ar * dt;
    fuck = 0.0 ;
    
    % Steam conditions at the pipe inlet
    rhos(1,1) = rhos0 ;
    Gs(1,1)   = Gs0 ;
    hs(1,1)   = hs0 ;
    Ts(1,1)   = Ts0 ;
    xs(1,1)   = xs0 ;
    us(1,1)   = Gs(1,1) / rhos(1,1) ;
    
    % PCM conditions at pipe inlet
    expp      = exp(-ap(1)) ;
    Tp(1,1)   = Ts(1,1) * (1. - expp) + Tp(1,2) * expp ; % From analytical solution
    if Tp(1,1) < Tpm - dTm
        cp(1,1)   = cps ;
        rhop(1,1) = rhops ;
        xp(1,1)   = 0.0 ;
    elseif Tp(1,1) > Tpm + dTm
        cp(1,1)   = cpl ;
        rhop(1,1) = rhopl ;
        xp(1,1)   = 1.0 ;
    else
        cp(1,1)   = cpm ;
        rhop(1,1) = rhops ;
        xp(1,1)   = (Tp(1,1) - Tpm + dTm) / (2. * dTm) ;
    end
        
    
    %*** PREDICTOR STEP ***%
    for i = 2 : Nx
        
        % PCM energy equation
        Tsguess = 0.5 * (Ts(i,2) + Ts(i-1,2)) ; % Steam temp probably doesn't change quickly - particularly when changing phase!
        Tp(i,1) = (Tp(i,2) + ap(i) * Tsguess) / (1. + ap(i)) ;
        PP = Tp(i,1) - Tp(i,2) ; % Predictor-PCM
        
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
        
        % Steam energy equation
        s1 = U * Ar * (Tp(i,1) - Tsguess) * dx / Gs(i-1,1) ;
        s2 = fuck*(rhos(i-1,1) * hs(i-1,1) - rhos(i-1,2) * hs(i-1,2)) * dx / (dt * Gs(i-1,1)) ;
        PS      = s1 - s2 ; % Predictor-Steam
        hs(i,1) = hs(i-1,1) + PS ;
                
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
        
        % Now estimate G from mass continuity equation
        PC = (rhos(i,2) - rhos(i,1)) * dx / dt ; % Predictor-Continuity
        Gs(i,1) = Gs(i-1,1);% + PC ;
        
        
    end
    
    % Now have updated estimates of rhos, Ts, Tp. Update others:
    ap = U * Apcm * dt ./ (rhop(:,1) .* cp(:,1)) ;
                
    % *** CORRECTOR STEP *** %
    for i = 2 : Nx
        
        % PCM energy equation
        CP = (Tp(i,2) + ap(i) * Ts(i,1)) / (1. + ap(i)) - Tp(i,2) ; % Corrector-PCM
        Tp(i,1) = Tp(i,2) + 0.5 * (PP + CP) ;
        
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
        
        % Steam energy equation
        s1 = U * Ar * (Tp(i,1) - Ts(i,1)) * dx / Gs(i-1,1) ;
        s2 = fuck*(rhos(i-1,1) * hs(i-1,1) - rhos(i-1,2) * hs(i-1,2)) * dx / (dt * Gs(i-1,1)) ;
        CS = s1 - s2 ; % Corrector-Steam
        
        hs(i,1) =  hs(i-1,1) + 0.5 * (PS + CS) ;        
        
        % Calculate steam properties
        if hs(i,1) > hv % steam
            xs(i,1)   = 1. ;
            Ts(i,1)   = RP1('HmassP_INPUTS',hs(i,1),Ps0,'T',steam) ;%Tsat + (hs(i,1)-hv) / csv  ; 
            rhos(i,1) = RP1('HmassP_INPUTS',hs(i,1),Ps0,'D',steam) ;%rhosv; 
        elseif hs(i,1) < hl % water
            xs(i,1)   = 1. ;
            Ts(i,1)   = hs(i,1) / csl ;
            rhos(i,1) = rhosl ;
        else % steam-water mixture
            xs(i,1) = (hs(i,1) - hl) / (hv - hl) ; % Estimate dryness fraction from lever rule
            Ts(i,1) = Tsat ;
            vs = xs(i,1) * (vv - vl) + vl ; % Specific volume from lever rule
            rhos(i,1) = 1./vs ;
        end
        
        % Continuity equation
        CC      = (rhos(i,2) - rhos(i,1)) * dx / dt ; % Corrector-Continuity
        Gs(i,1) = Gs(i-1,1);% + 0.5 * (PC + CC) ;
        
        
    end
    
    if n == 1
        keyboard
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
    
    % Save data points every so often
    if mod(n,1e5)==0
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
        
        nout = nout + 1 ;
    end
    
    time = time + dt ;
end

toc