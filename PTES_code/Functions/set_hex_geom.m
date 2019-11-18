function [HX] = set_hex_geom(fluidH, indH, fluidC, indC, eff_min, ploss_max, D1, t1_min, sigma, mode, var)
% Obtain geometric parameters based on performance objectives, using
% analytical solutions. It should be expected that objectives will be met
% accurately (only) when using fluids with small variations of
% thermophysical properties.

warning('*** further development needed for accurate computation - see notes below ***')
% Things to be improved:
% (1) Accurate computation of Nussel number for tube bundle inside
% developed_flow function, according to square/trilateral lattice
% configuration.
% (2) Minimum D2 value according to lattice configuration.

% Import fluid.state and fluid.stage
stateH = fluidH.state(indH(1),indH(2));
stateC = fluidC.state(indC(1),indC(2));

% Set inlet temperatures and pressures
THin = stateH.T;
pHin = stateH.p;
hHin = stateH.h;
TCin = stateC.T;
pCin = stateC.p;
hCin = stateC.h;
mH  = stateH.mdot;
mC  = stateC.mdot;

% Obtain temperature arrays as a function of the enthalpy arrays
[hvH,TvH] = get_h_T(fluidH,TCin-1,THin+1,pHin,100);
[hvC,TvC] = get_h_T(fluidC,TCin-1,THin+1,pCin,100);

% Obtain minimum and maximum enthalpy outlets and determine mean Cp
hHout_min = rtab1(TvH,hvH,TCin,1);
hCout_max = rtab1(TvC,hvC,THin,1);
CpHmean = (hHin - hHout_min)/(THin-TCin);
CpCmean = (hCout_max - hCin)/(THin-TCin);

% Determine mass flow rates
switch mode
    case 0
        % Both mass flow rates previously specified
        if any([mH,mC] == 0)
            error('mH and mC must be known if cond==0');
        end
        
    case 1
        % Only mass flow rate of hot fluid previously specified, compute
        % mass flow rate of cold fluid according to Crat
        if mH == 0, error('mH must be known if cond==1'); end
        Crat = var; %Crat = mH*CpH / (mC*CpC)
        mC = mH*CpHmean/(CpCmean*Crat);
        
    case 2
        % Only mass flow rate of cold fluid previously specified, compute
        % mass flow rate of hot fluid according to Crat        
        if mC == 0, error('mC must be known if cond==2'); end
        Crat = var; %Crat = mH*CpH / (mC*CpC)
        mH = Crat*mC*CpCmean/CpHmean;      
end

% Set the minimum number of transfer units that each stream should have to
% obtain the specified effectiveness
Ntu_min = 2/(1-eff_min);

% Declare the two fluid streams
S1 = stream;
S2 = stream;

% Determine which stream is the high pressure one (which flows inside the
% tubes - indicated 1) and which is the low pressure one (which flows on
% the shell side - indicated 2). Update the properties of each stream
% accordingly.
if pCin >= pHin
    % Cold stream is high pressure stream (tube side). Hot stream is low
    % pressure stream (shell side).
    S1.mdot = mC;
    S1.name = fluidC.name;
    S1.p    = pCin;
    S1.h    = 0.5*(hCout_max + hCin);
    S1 = stream_update(fluidC,S1,1);
    S2.mdot = mH;
    S2.name = fluidH.name;
    S2.p    = pHin;
    S2.h    = 0.5*(hHout_min + hHin);
    S2 = stream_update(fluidH,S2,1);
else
    % Hot stream is high pressure stream (tube side). Cold stream is low
    % pressure stream (shell side).
    S1.mdot = mH;
    S1.name = fluidH.name;
    S1.p    = pHin;
    S1.h    = 0.5*(hCout_max + hCin);
    S1 = stream_update(fluidH,S1,1);
    S2.mdot = mC;
    S2.name = fluidC.name;
    S2.p    = pCin;
    S2.h    = 0.5*(hHout_min + hHin);
    S2      = stream_update(fluidC,S2,1);
end

% Set pre-specified hydraulic diameter
S1.D = D1;

S1.Re = 100; %initial guess (assume laminar flow)
max_iter = 100;
tol = 1e-6;
for i=1:max_iter
    % Keep track of initial value (or value from previous iteration)
    Re1_0 = S1.Re;
    
    % Compute Cf and St
    [S1.Cf, S1.St] = developed_flow(S1.Re,S1.Pr,'circular');
    
    % Obtain the maximum mass flux that satisfies the Ntu and ploss
    % requirements
    S1.G = sqrt( (2*S1.p*S1.rho/Ntu_min) * (S1.St/S1.Cf) * ploss_max );
    
    % Compute the Reynolds number
    S1.Re = S1.G*S1.D/S1.mu;
    
    % Check convergence
    condition = abs((S1.Re - Re1_0)/Re1_0*100) < tol;
    if condition
        % Converged
        break
    end
end
if all([i>=max_iter,~condition])
    error('Re1 not converged')
end

% Compute heat transfer area and tube length
S1.A  = S1.mdot*Ntu_min/(S1.G*S1.St);
L     = S1.D*S1.A*S1.G/(4*S1.mdot);

% Compute flow area and number of tubes
S1.Af  = S1.mdot/S1.G; % total flow area (stream 1)
Af1one = pi/4*S1.D^2;    % flow area of one tube
N1     = round(S1.Af/Af1one); % number of tubes
S1.Af  = N1*Af1one;    % update

% Update G, Re, Cf and St
S1.G  = S1.mdot/S1.Af;
S1.Re = S1.G*S1.D/S1.mu;
[S1.Cf, S1.St] = developed_flow(S1.Re,S1.Pr,'circular');

% Set t1 (tube thickness)
t1_hoop = S1.p*S1.D/(2*sigma); %thickness for which sigma = hoop stress.
t1 = max([t1_min, 2*t1_hoop]);

% Compute A2
S2.A = L*N1*pi()*(S1.D+t1);

S2.Re = 100; %initial guess (assume laminar flow)
for i=1:max_iter
    % Keep track of initial value (or value from previous iteration)
    Re2_0 = S2.Re;
    
    % Compute Cf and St
    [S2.Cf, S2.St] = developed_flow(S2.Re,S2.Pr,'circular');
    
    % Obtain the maximum mass flux that satisfies the Ntu and ploss
    % requirements
    S2.G = sqrt( (2*S2.p*S2.rho/Ntu_min) * (S2.St/S2.Cf) * ploss_max );
    
    % Compute the hydraulic diameter
    S2.D = 4*L*S2.mdot / (S2.A*S2.G);
    
    % Ensure D2 is not smaller than D1 (space is needed between pipes)
    if S2.D < S1.D
        S2.D = S1.D;
        S2.G = 4*L*S2.mdot / (S2.A*S2.D);
    end
    
    % Ensure that Ntu2 and ploss2 satisfy conditions
    Ntu2   = 4*L/S2.D*S2.St
    ploss2 = Ntu2*S2.G^2*S2.v*(S2.Cf/S2.St)/(2*S2.p)
    if any([Ntu2<Ntu_min,ploss2>ploss_max])
        keyboard
    end
    
    
    % Update the Reynolds number
    S2.Re = S2.G*S2.D/S2.mu;
    
    % Check convergence
    condition1 = abs((S2.Re - Re2_0)/Re2_0*100) < tol;
    condition2 = i>=2;
    if all([condition1,condition2])
        % Converged
        break
    end
end
if all([i>=max_iter,~condition1])
    error('Re2 not converged')
end

% S2.D = 4*L*S2.St/Ntu_min;    
%     S2.G = 4*L*S2.mdot / (S2.A*S2.D);

% Find Af2
S2.Af = S2.mdot/S2.G;

Ntu1   = 4*L/S1.D*S1.St
Ntu2   = 4*L/S2.D*S2.St
ploss1 = Ntu1*S1.G^2*S1.v*(S1.Cf/S1.St)/(2*S1.p)
ploss2 = Ntu2*S2.G^2*S2.v*(S2.Cf/S2.St)/(2*S2.p)

keyboard

% Export results into HX structure
HX.shape     = 'circular';
HX.sigma     = sigma;      % Maximum allowable stress, Pa
HX.L         = L;          % Tube length, m
HX.D1        = S1.D;       % Tube diameter, m
HX.t1        = t1;         % Tube thickness, m
HX.AfT       = S1.Af + S2.Af; % Total flow area, m2
HX.AfR       = S2.Af/S1.Af;   % Ratio of flow areas, Af2/Af1, -

end


%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ hv, Tv ] = get_h_T( fluid, T1, T2, pressure, n )
% Obtain the hv and Tv arrays of a given fluid for the hex subroutines.
% Data ordered in regular intervals of h. T as a function of h.

if strcmp(fluid.read,'CP') %read from CoolProp
    
    h1  = CP1('PT_INPUTS',pressure,T1,'H',fluid.handle);
    h2  = CP1('PT_INPUTS',pressure,T2,'H',fluid.handle);
    hv  = linspace(h1,h2,n)';       % enthalpy array between TC1 and TH2
    pv  = ones(size(hv)).*pressure; % pressure array
    Tv  = CP1('HmassP_INPUTS',hv,pv,'T',fluid.handle); % temperature array
    
elseif strcmp(fluid.read,'TAB') %read from table
    
    Tx  = fluid.TAB(:,1);
    hy  = fluid.TAB(:,2);
    h1  = rtab1(Tx,hy,T1,0);
    h2  = rtab1(Tx,hy,T2,0);
    hv  = linspace(h1,h2,n)'; % enthalpy array between TC1 and TH2
    Tv  = rtab1(hy,Tx,hv,1);  % temperature array between TC1 and TH2
else
    error('not implemented')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%