function [ENV_MIX,stream,iNew] = env_mix(ENV_MIX,stream,indM,T0,p0)

% One fluid stream is brought into thermal equilibrium with the environment
% iM is the index of the original stream

% Import fluid.state and fluid.stage
stateM = stream.state(indM(1),indM(2));
stageM = stream.stage(indM(1),indM(2));

% Extract initial values of enthalpy and entropy
h_st = stateM.h;
s_st = stateM.s;

% Check that pressures are the same
if abs((stateM.p - p0)/p0) > 1e-9
    warning('stream must have the same pressure as the environment')
    keyboard
    stateM.p = p0 ;
end

% Update states
stateM.T = T0;
[stateM] = update_state(stateM,stream,1);

% Compute stages
% Entropy change
Ds         = stateM.s - s_st;
Dh         = stateM.h - h_st;
% Main stream
stageM.Dh   = stateM.h - h_st;
stageM.sirr = Ds-Dh/T0;
stageM.q    = stageM.Dh;
stageM.w    = 0;
stageM.type = 'env_mix';

% Save losses into ENV_MIX object
ENV_MIX.Dh(indM(1),1)   = stageM.Dh;
ENV_MIX.sirr(indM(1),1) = stageM.sirr;
ENV_MIX.q(indM(1),1)    = ENV_MIX.Dh(indM(1),1);
ENV_MIX.w(indM(1),1)    = 0;
ENV_MIX.mdot(indM(1),1) = stateM.mdot;

% Cold stream
%ENV_MIX.Dh(indM(1),2)   = stateS.h - h2_st;
%ENV_MIX.sirr(indM(1),2) = (stateM.mdot*Ds + stateS.mdot*Ds2)/stateS.mdot;
%ENV_MIX.q(indM(1),2)    = ENV_MIX.Dh(indM(1),2);
%ENV_MIX.w(indM(1),2)    = 0;
%ENV_MIX.mdot(indM(1),2) = stateS.mdot;

% Reassing values to fluid.state and fluid.stage
stream.state(indM(1),indM(2)+1) = stateM;
stream.stage(indM(1),indM(2))   = stageM;

% Increase stage counter
iNew = indM(2) + 1;
end