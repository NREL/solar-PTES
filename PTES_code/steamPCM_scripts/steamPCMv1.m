%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCM with steam
% This code models a Phase Change Material (PCM) which surrounds a pipe 
% through which evaporating or condensing steam passes. The code uses 
% routines developed by Pau Farres-Antunez (Cambridge University) and
% Josh McTigue (NREL). The method is based partly on the work of Sharan et
% al. (2018) and partly on McTigue's PhD thesis.
%
% Author: Josh McTigue, Pau Farres-Antunez
% 31 March 2020
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

% Load CoolProp library (low-level interface)
load_coolprop

%**** INPUT PARAMETERS ****%
dp    = 0.1 ;  % Pipe diameter
Lp    = 1.0 ; % Pipe length
U     = 5e5 ; % Overall heat transfer coefficient, W/m2K
Psat  = 20e5 ; % Saturation pressure
x0    = 1.0 ;  % Initial dryness fraction
mdot0 = 1.0;  % Mass flow rate, kg/s
CFL   = 0.5 ; % Courant-Freidrich-Lewy number
Nx    = 100 ; % Number of gridsteps
tN    = 1*1000;  % Duration of simultion, s 

%**** SET UP SOME PARAMETERS ****%

% Set up steam class
steam = fluid_class('Water','WF','CP','HEOS',2,30);
Ap    = 0.25 * pi * dp^2 ; % Pipe area
Ar    = 4.0 / dp ;
Tsat = RP1('PQ_INPUTS',Psat,x0,'T',steam) ; % Saturation temperature

% Initial properties
T0 = Tsat ;
Tpcm = Tsat - 15 ;
P0 = Psat ;
h0 = RP1('PQ_INPUTS',Psat,x0,'H',steam) ; 
hv = RP1('PQ_INPUTS',Psat,1.0,'H',steam) ; % Enthalpy of vapour
hl = RP1('PQ_INPUTS',Psat,0.0,'H',steam) ; % Enthalpy of liquid

vv = 1./RP1('PQ_INPUTS',Psat,1.0,'D',steam) ; % Specific volume of vapour
vl = 1./RP1('PQ_INPUTS',Psat,0.0,'D',steam) ; % Specific volume of liquid

s0 = RP1('PQ_INPUTS',Psat,x0,'S',steam) ; 
rho0 = RP1('PQ_INPUTS',Psat,x0,'D',steam) ; 
G0   = mdot0 / Ap ; % Mass flux per unit area
u0   = G0 / rho0 ;

dx   = Lp / Nx ; % Grid step size
dt   = CFL * dx / u0 ; % Time step size
Nt   = int64(tN / dt) ; % Number of timesteps

% Set up some property arrays. 1st column is new value, 2nd column is old.
rho  = rho0 * ones(Nx,2) ;
G    = G0   * ones(Nx,2) ;
h    = h0   * ones(Nx,2) ;
s    = s0   * ones(Nx,2) ;
x    = x0   * ones(Nx,2) ;

time = 0.;
nout = 1 ;

%**** STEP FORWARDS IN TIME ****%
for n = 1 : Nt
    
    % Conditions at the pipe inlet
    rho(1,1) = rho0 ;
    G(1,1)   = G0 ;
    h(1,1)   = h0 ;
    x(1,1)   = x0 ;
    
    %** STEP FORWARD ALONG THE PIPE **%
    for i = 2 : Nx
        
        % Guess steam density
        rho(i,1) = 1.1*rho(i,2);%0.5 * (rho(i,2) + rho(i+1,2)) ;
        
        % Do this twice - basically refines guessed value of rho
        for k = 1:2
            % Continuity equation
            C1 = (rho(i,2) - rho(i,1)) * dx / dt ;
            G(i,1) = G(i-1,1) + C1 ;
            
            % Steam energy equation
            E1 = U * Ar * (Tpcm - Tsat) * dt ;
            E2 = G(i-1,1) * h(i-1,1) * dt / dx ;
            E3 = rho(i,2) * h(i,2) + E1 + E2 ;
            E4 = (rho(i,1) + G(i,1) * dt / dx) ;
            
            h(i,1) = E3 / E4 ;
            
            x(i,1) = (h(i,1) - hl) / (hv - hl) ; % Estimate dryness fraction from lever rule
            v = x(i,1) * (vv - vl) + vl ; % Specific volume from lever rule
            rho(i,1) = 1./v ;
            
            %x(i,1)   = RP1('HmassP_INPUTS',h(i,1),Psat,'Q',steam) ;
            %s(i,1)   = RP1('HmassP_INPUTS',h(i,1),Psat,'S',steam) ;
            %rhoaccurate = RP1('HmassP_INPUTS',h(i,1),Psat,'D',steam) ;
        end
    end
    
    rho(:,2) = rho(:,1) ;
    G(:,2)   = G(:,1) ;
    h(:,2)   = h(:,1) ;
    x(:,2)   = x(:,1) ;
    
    % Save data points every so often
    if mod(n,1e5)==0
        fprintf("TIMESTEP #%i COMPLETED.\n",n);
        Gsave(:,nout) = G(:,1) ;
        Hsave(:,nout) = h(:,1) ;
        Xsave(:,nout) = x(:,1) ;
        RHOsave(:,nout) = rho(:,1) ;
        Ssave(:,nout) = s(:,1) ;
        nout = nout + 1 ;
    end
    
    time = time + dt ;
end





