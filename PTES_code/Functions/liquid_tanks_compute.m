function [ A2, B1, B2, WL ] = liquid_tanks_compute( fluid_streams, i1, i2, A1, t, T0)
% Computes the conditions of the source and sink liquid storage tanks,
% according to the streams flowing from source to sink. Also computes the
% lost work of the mixing process at the inlet of the sink tank. A is the
% source tank, B is the sink tank. 1=start, 2=end. t is time.
% The subroutine assumes that the sink tank was initially empty (B1.M=0).
% i1 distinguishes the charge row to the discharge row of fluid states
% i2 distinguishes the inlet column (1) to the outlet column (2) of fluid states

% Give A2, B1 and B2 the same class as A1
A2 = A1; B1 = A1; B2 = A1;

% Compute total mass, enthalpy and entropy flows of the fluid_streams
Mdot = 0;
Hdot = 0;
Sdot = 0;
n = numel(fluid_streams);
for i=1:n
    Mdot   = Mdot + fluid_streams(i).state(i1,i2).mdot;
    Hdot   = Hdot + fluid_streams(i).state(i1,i2).h*fluid_streams(i).state(i1,i2).mdot;
    Sdot   = Sdot + fluid_streams(i).state(i1,i2).s*fluid_streams(i).state(i1,i2).mdot;
end

% Find specific enthalpy and pressure of combined stream
h_mix = Hdot/Mdot;
p_mix = fluid_streams(1).state(i1,i2).p;

% Set conditions of combined stream
mixed_stream_state = fluid_streams(1).state(i1,i2);
mixed_stream_state.mdot = Mdot;
mixed_stream_state.h = h_mix;
mixed_stream_state.p = p_mix;

% Update properties based on enthalpy and pressure
mixed_stream_state = update_state(mixed_stream_state,fluid_streams(1).handle,fluid_streams(1).read,fluid_streams(1).TAB,2);
T_mix = mixed_stream_state.T;
v_mix = 1/mixed_stream_state.rho;
s_mix = mixed_stream_state.s;

% Update end conditions of source tank A
A2.M = A1.M - Mdot*t;
A2.T = A1.T;
A2.V = A1.V*A2.M/A1.M;
A2.H = A1.H*A2.M/A1.M;
A2.S = A1.S*A2.M/A1.M;
A2.B = A1.B*A2.M/A1.M;

% Update conditions of sink tank B
% Start
B1.T = T_mix;
B1.M = 0;
B1.V = 0;
B1.H = 0;
B1.S = 0;
B1.B = 0;
% End
B2.M = Mdot*t;
B2.T = T_mix;
B2.V = v_mix*B2.M;
B2.H = h_mix*B2.M;
B2.S = s_mix*B2.M;
B2.B = B2.H - T0*B2.S;

% Compute entropy generation of mixing
Sdot_irr = s_mix*Mdot - Sdot;
WL = T0*Sdot_irr*t;

end