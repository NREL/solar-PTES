function [ A ] = liquid_tank_start( fluid_streams, ind, A, t, T0 )

% Computes the necessary initial fluid mass (and other initial conditions)
% for a tank A to supply the outgoing fluid_streams during a time t.
% The different fluid streams may have different mass flow rates,
% but they all depart with the same initial conditions.

% Obtain linear index of interest
idx = sub2ind(size(fluid_streams(1).state),ind(1),ind(2));

% Compute total mass flow rate of the fluid_streams
Mdot = 0;
n = numel(fluid_streams);
for i=1:n
    Mdot   = Mdot + fluid_streams(i).state(idx).mdot;
end

% Set conditions of source tank A
A.M = Mdot*t;
A.V = (1/fluid_streams(1).state(idx).rho)*A.M;
A.H = fluid_streams(1).state(idx).h*A.M;
A.S = fluid_streams(1).state(idx).s*A.M;
A.B = A.H - T0*A.S; %fluid_out.b*A.M;

end

