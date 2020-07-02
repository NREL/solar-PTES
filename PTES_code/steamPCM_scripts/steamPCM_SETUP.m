% Set up the geometry of the PCM tube based on charging conditions

% Areas
Ap    = 0.25 * pi * dp^2 ; % Pipe area
Ar    = 4.0 / dp ; % Heat transfer area

% Iterate to find do, mdot, and U %
do_lo = 1.01 * dp ; % Guess do
do_hi = 10 * dp ; % Guess do
PrC   = RP1('PQ_INPUTS',PsC,0.0,'PRANDTL',steam) ; % Prandtl number of liquid

% Now iterate to find do
cnt = 1 ;
cntN = 100 ;
err  = 1e10 ;
tol  = 0.01 ;
while err > tol && cnt <= cntN
    
   do_guess = 0.5 * (do_lo + do_hi) ; 
   
   Apcm  = pi * (do_guess^2 - dp^2) / 4 ; % This is actual PCM area per unit length
   mdotN = rhops * Lp * Lpcm * Apcm / (tN * (hsC - hlC)) ;
   U1    = mdotN * (hsC - hlC) / (pi * dp * Lp * (TsC - Tpm)) ;
   
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
do    = do_guess ;
Apcm  = pi * (do^2 - dp^2) / 4 ; % PCM area per unit length
mdotN = rhops * Lp * Lpcm * Apcm / (tN * (hsC - hlC)) ; % Nominal charging mass flow rate
mdot  = mdotN * ones(Nload,1) ; % Mass flow rate of each load number
Apcm  = 4.0 / (dp * ((do/dp)^2 - 1.)) ; % Normalized PCM area for computations
GsC   = mdotN / Ap ; % Mass flux per unit area
usC   = GsC / rhosC ;

if Lreadload
    Npipe = tes_mdot0 / mdotN ; % Number of steamPCM pipes required
end

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
PrL   = RP1('PQ_INPUTS',PsC,0.0,'PRANDTL',steam) ; % Prandtl number of liquid
if ReL < 3000
    NuL = 4.36;
else
    NuL = 0.023 * ReL^0.8 * PrL^0.4 ; % Dittus-Boelter
end
hcL   = NuL * ktl / dp ;

ReV   = GsC * dp / musv ; % Reynolds number based on all liquid flow
PrV   = RP1('PQ_INPUTS',PsC,1.0,'PRANDTL',steam) ; % Prandtl number of liquid
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
