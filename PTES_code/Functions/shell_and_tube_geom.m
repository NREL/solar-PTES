function [ C, H, HX] = shell_and_tube_geom( C, H, HX)
% SHELL-AND-TUBE GEOMETRY. Compute and assign derived geometrical
% parameters to each side of a shell-and-tube heat exchanger.
% INPUT parameters are the inlet pressures (C.pin and H.pin), length (L),
% total flow area (AfT), ratio of flow areas (AfR), tube diameter (D1) and
% tube thickness (t1).
% OUTPUT parameters are the number of tubes (N1), the hydraulic diameter of
% the shell-side (D2), the flow areas of each side, the mass fluxes (G1 and
% G2), the heat transfer areas (A1 and A2), and the combined volume of
% metal of the tubes (Vm).
% All parameters are extracted from and stored among three structures:
% C (cold stream), H (hot stream) and HX (heat exchanger).


% Determine which stream is the high pressure one (which flows inside the
% tubes - indicated 1) and which is the low pressure one (which flows on
% the shell side - indicated 2)
if C.pin >= H.pin
    % Cold stream is high pressure stream (tube side). Hot stream is low
    % pressure stream (shell side).
    p1 = C.pin;
    m1 = C.mdot;
    m2 = H.mdot;
elseif C.pin < H.pin
    % Hot stream is high pressure stream (tube side). Cold stream is low
    % pressure stream (shell side).
    p1 = H.pin;
    m1 = H.mdot;
    m2 = C.mdot;
else
    error('C.pin and H.pin must be defined')
end

% Extract parameters from HX structure
L     = HX.L;   % Length, m
AfT   = HX.AfT; % Total flow area = Af1 + Af2, m2
AfR   = HX.AfR; % Flow area ratio = Af2/Af1
D1    = HX.D1;  % Tube diameter, m
t1    = 0; % negligible thickness assumed at the moment
% sigma = HX.sigma; % Maximum allowable stress, Pa
% t1_min    = HX.t1_min;    % Minimum tube thickness, m
% t1_D1_min = HX.t1_D1_min; % Minimum tube thickness-to-diameter ratio
% 
% % Check that tube thickness is sufficient to withstand internal pressure
% t1 = max([t1_min,t1_D1_min*D1,D1*p1/sigma]);

% Obtain flow areas and number of tubes
Af1  = AfT/(AfR+1);
Af1one = pi/4*D1^2; % flow area of one tube
N1   = round(Af1/Af1one); %number of tubes
Af1  = N1*Af1one; %update
Af2  = AfT - Af1;

% Obtain mass fluxes
G1 = m1/Af1;
G2 = m2/Af2;

% Obtain heat transfer areas
A1 = pi*D1*L*N1;
A2 = pi*(D1+t1)*L*N1;

% Compute hydraulic diameter shell side (Dh = 4*Af/P by definition)
D2 = 4*Af2/(pi*(D1+t1)*N1); %must check if this is independent of lattice configuration

% Compute volume of metal of all tubes
Vm = (pi/4*(D1+t1)^2*N1 - Af1)*L;

% Export computed variables for output
HX.t1  = t1;
HX.N1  = N1;
HX.D2  = D2;
HX.G1  = G1;
HX.G2  = G2;
HX.Af1 = Af1;
HX.Af2 = Af2;
HX.A1  = A1;
HX.A2  = A2;
HX.Vm  = Vm;
if C.pin >= H.pin
    % Cold stream is high pressure stream
    C.D = D1;
    C.G = G1;
    C.A = A1;
    
    H.D = D2;
    H.G = G2;
    H.A = A2;
else
    % Hot stream is high pressure stream
    C.D = D2;
    C.G = G2;
    C.A = A2;
    
    H.D = D1;
    H.G = G1;
    H.A = A1;
end

end

