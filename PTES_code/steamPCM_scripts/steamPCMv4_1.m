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
% In v4.1, the PCM equation is written in terms of enthalpy which improves
% the stability of the solution. Have also added charge and discharge
% modes.
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

% STEAM CHARGING CONDITIONS
steam = fluid_class('Water','WF','CP','HEOS',2,30);
PsatC = 20e5 ; % Saturation pressure
xsC   = 0.95 ;  % Initial dryness fraction of steam
TsatC = RPN('PQ_INPUTS',PsatC,1.0,'T',steam) ; % Saturation temperature
TsC   = TsatC ;
PsC   = PsatC ;

% STEAM DISCHARGING CONDITIONS
PsatD = 12e5 ; % Saturation pressure
xsD   = 0.0 ;  % Initial dryness fraction of steam
TsatD = RPN('PQ_INPUTS',PsatD,0.0,'T',steam) ; % Saturation temperature
TsD   = TsatD ;
PsD   = PsatD ;

if PsatD > PsatC
    warning('Discharging steam pressure is larger than charging pressure. May cause problems');
end

% PCM
Tpm   = TsC - 10 ; % Melting point of PCM
TpC   = TsC - 5 ; % Charged temperature of PCM
xpC   = 1.0 ;
TpD   = TsC - 12 ; % Discharged temperature of PCM
xpD   = 0.0 ;
cpl   = 2.5e3 ;     % Specific heat capacity of liquid PCM
cps   = 1.5e3 ;     % Specific heat capacity of solid PCM
rhopl = 1e3 ;       % Density of liquid PCM
rhops = 0.8e3 ;     % Density of solid PCM
Lpcm  = 2e6 ;       % Latent heat of PCM
kpcm  = 5 ;         % Thermal conductivity of PCM

cpD   = cps;
rhopD = rhops ;

hps   = cps * Tpm ;
hpl   = hps + Lpcm ;
if TpD < Tpm
    hpD   = cpD * TpD ;
else
    error('Not implemented')
end

% GRID AND TIME STEPS
CFL   = 100 ;      % Courant-Freidrich-Lewy number
Nx    = 1000 ;    % Number of gridsteps
tN    = 4*3600;  % Nominal duration of charge, s
Nload = 2 ;      % Number of loads (charge and discharge are separate loads)
Iload = 1 ;      % Load number counter
Tload = tN * ones(Nload,1) ; % Array that contains the duration of each load
Load  = ["c";"d";"c"]; % 'c' - charge. 'd' - discharge, duh. 
PCMiterations = 2 ; % Number of times to iterate PCM equation. 2 gives good results. 1 is satisfactory.

%**** SET UP SOME PARAMETERS ****%

% Initial properties of steam 
hvC = RPN('PQ_INPUTS',PsatC,1.0,'H',steam) ; % Enthalpy of vapour
hlC = RPN('PQ_INPUTS',PsatC,0.0,'H',steam) ; % Enthalpy of liquid
vvC = 1./RPN('PQ_INPUTS',PsatC,1.0,'D',steam) ; % Specific volume of vapour
vlC = 1./RPN('PQ_INPUTS',PsatC,0.0,'D',steam) ; % Specific volume of liquid

hvD = RPN('PQ_INPUTS',PsatD,1.0,'H',steam) ; % Enthalpy of vapour
hlD = RPN('PQ_INPUTS',PsatD,0.0,'H',steam) ; % Enthalpy of liquid
vvD = 1./RPN('PQ_INPUTS',PsatD,1.0,'D',steam) ; % Specific volume of vapour
vlD = 1./RPN('PQ_INPUTS',PsatD,0.0,'D',steam) ; % Specific volume of liquid

% Going to do something a bit crappy to avoid having to use Coolprops (very slow)
Ntab  = 50 ;
htabC = linspace(hvC,2*hvC,Ntab)' ;
htabD = linspace(hvD,2*hvD,Ntab)' ;
TtabC = zeros(Ntab,1) ;
VtabC = zeros(Ntab,1) ;
TtabD = zeros(Ntab,1) ;
VtabD = zeros(Ntab,1) ;

for i = 1 : Ntab
   TtabC(i) = RPN('HmassP_INPUTS',htabC(i),PsatC,'T',steam);
   VtabC(i) = 1./RPN('HmassP_INPUTS',htabC(i),PsatC,'D',steam);
   
   TtabD(i) = RPN('HmassP_INPUTS',htabD(i),PsatD,'T',steam);
   VtabD(i) = 1./RPN('HmassP_INPUTS',htabD(i),PsatD,'D',steam);
end

% Correlation between H and T, and H and V (specific volume)
HH = [ones(Ntab,1) htabC] ;
TcoefC = HH \ TtabC ;
VcoefC = HH \ VtabC ;

HH = [ones(Ntab,1) htabD] ;
TcoefD = HH \ TtabD ;
VcoefD = HH \ VtabD ;

hsC   = RPN('PQ_INPUTS',PsC,xsC,'H',steam) ;
ssC   = RPN('PQ_INPUTS',PsC,xsC,'S',steam) ; 
rhosC = RPN('PQ_INPUTS',PsC,xsC,'D',steam) ; 
musC  = RPN('PQ_INPUTS',PsC,xsC,'VISCOSITY',steam) ; 

hsD   = RPN('PQ_INPUTS',PsD,xsD,'H',steam) ;
ssD   = RPN('PQ_INPUTS',PsD,xsD,'S',steam) ; 
rhosD = RPN('PQ_INPUTS',PsD,xsD,'D',steam) ; 
musD  = RPN('PQ_INPUTS',PsD,xsD,'VISCOSITY',steam) ; 

csl   = RPN('PQ_INPUTS',PsC,0.0,'CPMASS',steam) ; 
rhosl = RPN('PQ_INPUTS',PsC,0.0,'D',steam) ; 
musl  = RPN('PQ_INPUTS',PsC,0.0,'VISCOSITY',steam) ; 
ktl   = RPN('PQ_INPUTS',PsC,0.0,'L',steam) ;  % Thermal conductivity of liquid

csv   = RPN('PQ_INPUTS',PsC,1.0,'CPMASS',steam) ; 
rhosv = RPN('PQ_INPUTS',PsC,1.0,'D',steam) ; 
musv  = RPN('PQ_INPUTS',PsC,1.0,'VISCOSITY',steam) ; 
ktv   = RPN('PQ_INPUTS',PsC,1.0,'L',steam) ;  % Thermal conductivity of liquid

Psc   = RPN(0,0,0,'PCRIT',steam);      % Critical pressure of steam

% Areas
Ap    = 0.25 * pi * dp^2 ; % Pipe area
Ar    = 4.0 / dp ; % Heat transfer area

%%% NOW WE MUST ITERATE TO FIND do, mdot, and U %%%

do_lo = 1.01 * dp ; % Guess do
do_hi = 10 * dp ; % Guess do
PrC   = RPN('PQ_INPUTS',PsC,0.0,'PRANDTL',steam) ; % Prandtl number of liquid

% Now iterate to find do
cnt = 1 ;
cntN = 100 ;
err  = 1e10 ;
tol  = 0.01 ;
while err > tol && cnt <= cntN
    
   do_guess = 0.5 * (do_lo + do_hi) ; 
   
   Apcm  = pi * (do_guess^2 - dp^2) / 4 ; % This is actual PCM area per unit length
   mdotN = rhops * Lp * Lpcm * Apcm / (tN * (hvC - hlC)) ;
   U1    = mdotN * (hvC - hlC) / (pi * dp * Lp * (TsC - Tpm)) ;
   
   % Also calculate U from heat transfer terms
   Rpipe = 0.5 * dp * log (do_guess/dp) / kpcm ;
   GsC   = mdotN / Ap ; % Mass flux per unit area
   ReC   = GsC * dp / musl ; % Reynolds number based on all liquid flow
   if ReC < 3000
       NuC = 4.36 ;
   else
       NuC = 0.023 * ReC^0.8 * PrC^0.4 ; % Dittus-Boelter
   end
   fac   = 0.55 + 2.09 / (PsC / Psc)^0.38 ;
   hcC   = fac * NuC * ktl / dp ;
   U2    = 1. / (1./hcC + Rpipe) ;
   
   Uguess = U2 - U1 ;
   
   cnt = cnt + 1 ;
   err = abs(Uguess) ;
   
   if Uguess > 0
       do_lo = do_guess ;
   else
       do_hi = do_guess ;
   end
   
end

% Have now established correct values so assign them
do = do_guess ;
Apcm = pi * (do^2 - dp^2) / 4 ; % PCM area per unit length
mdotN = rhops * Lp * Lpcm * Apcm / (tN * (hvC - hlC)) ; % Nominal charging mass flow rate
mdot  = mdotN * ones(Nload,1) ; % Mass flow rate of each load number
mdot(2) = 2.5 * mdotN ;
Apcm  = 4.0 / (dp * ((do/dp)^2 - 1.)) ; % Normalized PCM area for computations
GsC   = mdotN / Ap ; % Mass flux per unit area
usC   = GsC / rhosC ;

% GRID
dx   = Lp / Nx ; % Grid step size
dt   = CFL * dx / usC ; % Time step size
Nt   = int64(tN / dt) ; % Number of timesteps

% Set up some property arrays. 1st column is new value, 2nd column is old.
% steam properties
rhos  = rhosC * ones(Nx,2) ; % 's' stands for steam
Gs    = GsC   * ones(Nx,2) ; % Mass flux
Ts    = TsC   * ones(Nx,2) ; % Temperature
Ps    = PsC   * ones(Nx,2) ; % Pressure
hs    = hsC   * ones(Nx,2) ; % Enthalpy
ss    = ssC   * ones(Nx,2) ; % Entropy
xs    = xsC   * ones(Nx,2) ; % Dryness fraction 
mus   = musC  * ones(Nx,2) ; % Viscosity

% Some dimensionless numbers
ReL   = GsC * dp / musl ; % Reynolds number based on all liquid flow
PrL   = RPN('PQ_INPUTS',PsC,0.0,'PRANDTL',steam) ; % Prandtl number of liquid
if ReL < 3000
    NuL = 4.36;
else
    NuL = 0.023 * ReL^0.8 * PrL^0.4 ; % Dittus-Boelter
end
hcL   = NuL * ktl / dp ;

ReV   = GsC * dp / musv ; % Reynolds number based on all liquid flow
PrV   = RPN('PQ_INPUTS',PsC,1.0,'PRANDTL',steam) ; % Prandtl number of liquid
NuV   = 0.023 * ReV^0.8 * PrV^0.4 ; % Dittus-Boelter
hcV   = NuV * ktv / dp ;

ReC   = ReV ;
PrC   = PrV ;
fac   = 0.55 + 2.09 / (PsC / Psc)^0.38 ;
hcC   = fac * NuL * ktl / dp ;
NuC   = hcC * dp / ktv ;
StC   = NuC / (ReC * PrC) ;
UC    = 1. / (1./hcC + Rpipe) ;

Re    = ReC * ones(Nx,2) ; % Reynolds number
Pr    = PrC * ones(Nx,2) ; % Prandtl number
Nu    = NuC * ones(Nx,2) ; % Nusselt number
St    = StC * ones(Nx,2) ; % Stanton number
hc    = hcC * ones(Nx,2) ; % Steam heat transfer coefficient

% pcm properties
rhop  = rhopD * ones(Nx,2) ; % 'p' stands for PCM
Tp    = TpD   * ones(Nx,2) ; 
hp    = hpD   * ones(Nx,2) ; 
%sp    = sp0   * ones(Nx,2) ;
xp    = xpD   * ones(Nx,2) ;
cp    = cpD   * ones(Nx,2) ; 

Ux    = UC   * ones(Nx,2) ; 

% Values to save
Nprof = 20 ;
GSsave   = zeros(Nx,Nprof,Nload) ;
HSsave   = zeros(Nx,Nprof,Nload) ;
TSsave   = zeros(Nx,Nprof,Nload) ;
XSsave   = zeros(Nx,Nprof,Nload) ;
RHOSsave = zeros(Nx,Nprof,Nload) ;
SSsave   = zeros(Nx,Nprof,Nload) ;

TPsave   = zeros(Nx,Nprof,Nload) ;
HPsave   = zeros(Nx,Nprof,Nload) ;
XPsave   = zeros(Nx,Nprof,Nload) ;
RHOPsave = zeros(Nx,Nprof,Nload) ;

UXsave   = zeros(Nx,Nprof,Nload) ;
    
tic

%**** EXECUTE EACH LOAD CYCLE ****%

while Iload <= Nload

    %**** STEP FORWARDS IN TIME ****%
    if strcmp(Load(Iload),'c')
        fprintf('\nLOAD CYCLE #%i. CHARGING.\n\n',Iload);
    else
        fprintf('\nLOAD CYCLE #%i. DISCHARGING.\n\n',Iload);
    end
    
    % Set steam properties as constant along pipe
    % Inlet properties of steam are known
    if strcmp(Load(Iload),'c')
        rhos  = rhosC * ones(Nx,2) ; % 's' stands for steam
        Gs    = mdot(Iload) * ones(Nx,2) / Ap ; % Mass flux
        Ts    = TsC   * ones(Nx,2) ; % Temperature
        Ps    = PsC   * ones(Nx,2) ; % Pressure
        hs    = hsC   * ones(Nx,2) ; % Enthalpy
        ss    = ssC   * ones(Nx,2) ; % Entropy
        xs    = xsC   * ones(Nx,2) ; % Dryness fraction
        mus   = musC  * ones(Nx,2) ; % Viscosity
        us    = Gs ./ rhos ;
        
        hl = hlC ;
        hv = hvC ;
        vl = vlC ;
        vv = vvC ;
        Tsat = TsC ;
                
        Tcoef = TcoefC ;
        Vcoef = VcoefC ;
    elseif strcmp(Load(Iload),'d')
        rhos  = rhosD * ones(Nx,2) ; % 's' stands for steam
        Gs    = mdot(Iload) * ones(Nx,2) / Ap ; % Mass flux
        Ts    = TsD   * ones(Nx,2) ; % Temperature
        Ps    = PsD   * ones(Nx,2) ; % Pressure
        hs    = hsD   * ones(Nx,2) ; % Enthalpy
        ss    = ssD   * ones(Nx,2) ; % Entropy
        xs    = xsD   * ones(Nx,2) ; % Dryness fraction
        mus   = musD  * ones(Nx,2) ; % Viscosity
        us    = Gs ./ rhos ;
        
        hl = hlD ;
        hv = hvD ;
        vl = vlD ;
        vv = vvD ;
        Tsat = TsD ;
        
        Tcoef = TcoefD ;
        Vcoef = VcoefD ;
    else
        error('You`ve mucked up')
    end
            
    n = 1; % Timestep number
    time = 0.;
    nout = 1 ;
    while time <= Tload(Iload)
        
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
        
        %** CALCULATE STEAM PROPERTIES **%
        % Inlet properties of steam are known
        Gs(1,1)   = mdot(Iload) / Ap ;
        if strcmp(Load(Iload),'c')
            hs(1,1)   = hsC ;
            Ts(1,1)   = TsC ;
            xs(1,1)   = xsC ;
            Ps(1,1)   = PsC ;
            rhos(1,1) = rhosC ;
            us(1,1)   = Gs(1,1) / rhos(1,1) ;
        elseif strcmp(Load(Iload),'d')
            hs(1,1)   = hsD ;
            Ts(1,1)   = TsD ;
            xs(1,1)   = xsD ;
            Ps(1,1)   = PsD ;
            rhos(1,1) = rhosD ;
            us(1,1)   = Gs(1,1) / rhos(1,1) ;
        else
            error('You`ve mucked up')
        end
        
        % Heat transfer conditions at pipe inlet
        if xs(1,1) >=1
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
            fac(1)  = (1-xs(1,1))^0.8 + 3.8 * xs(1,1)^0.76 * (1-xs(1,1))^0.04 / (Ps(1,1) / Psc)^0.38 ;
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
                    %Ts(i,1)   = RPN('HmassP_INPUTS',hs(i,1),Ps(i,1),'T',steam);
                    Ts(i,1)   = Tsat + (hs(i,1)-hv) / csv ;
                    %rhos(i,1) = RPN('HmassP_INPUTS',hs(i,1),Ps(i,1),'D',steam);% rhosv;
                    
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
    
    Iload = Iload + 1 ;

end

toc