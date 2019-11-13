function [stream,iM,iS] = mix_streams(stream,indM,indS)

% Two fluid streams are brought into thermal equilibrium. The second stream
% is incorporated into the main stream.
% stream is the class of the fluid stream
% iM is the index of the first stream (main)
% iS is the index of the second stream (second)

% Import fluid.state and fluid.stage
stateM = stream.state(indM(1),indM(2));
stageM = stream.stage(indM(1),indM(2));
stateS = stream.state(indS(1),indS(2));
stageS = stream.stage(indS(1),indS(2));

% Extract initial values of enthalpy and entropy
h1_st = stateM.h;
h2_st = stateS.h;
s1_st = stateM.s;
s2_st = stateS.s;

% Check that fluids and pressures are the same
if strcmp(stream.name,stream.name)
    % This is redundant now, because are the same by definition
else
    error('stream1 and stream2 must be the same fluid')
end
if abs((stateM.p - stateS.p)/stateM.p) > 1e-10
    error('stream1 and stream2 must have the same pressure')
end
p = stateM.p;

% Mass conservation
Mdot = stateM.mdot + stateS.mdot;

% Energy conservation
Hdot = stateM.mdot*stateM.h + stateS.mdot*stateS.h;
h    = Hdot/Mdot;

% Update states
stateM.h = h;
[stateM] = update_state(stateM,stream.handle,stream.read,stream.TAB,2);
stateS.h = h;
[stateS] = update_state(stateS,stream.handle,stream.read,stream.TAB,2);

% Compute stages
% Entropy change
Ds1         = stateM.s - s1_st;
Ds2         = stateS.s - s2_st;
% Main stream
stageM.Dh   = stateM.h - h1_st;
%stageM.Dh   = (stateM.mdot*(stateM.h - h1_st) + stateS.mdot*(stateS.h - h2_st))/Mdot;
stageM.sirr = (stateM.mdot*Ds1 + stateS.mdot*Ds2)/stateM.mdot;
stageM.q    = stageM.Dh;
stageM.w    = 0;
stageM.type = 'mixing';
% Cold stream
stageS.Dh   = stateS.h - h2_st;
stageS.sirr = 0; %avoid counting the lost work twice
stageS.q    = 0;
stageS.w    = 0;
stageS.type = 'mixing';

% Update mass flow rates (join streams after reaching equilibrium)
stateM.mdot = Mdot;
stateS.mdot = 0;

% Reassing values to fluid.state and fluid.stage
stream.state(indM(1),indM(2)+1) = stateM;
stream.stage(indM(1),indM(2))   = stageM;
stream.state(indS(1),indS(2)+1) = stateS;
stream.stage(indS(1),indS(2))   = stageS;

% Increase stage counter
iM = indM(2) + 1;
iS = indS(2) + 1;
end