function [main,second,iM,iS] = mix_streams(main,indM,second,indS)

% Two fluid streams are brought into thermal equilibrium. The second stream
% is incorporated into the main stream.

% Import fluid.state and fluid.stage
stateM = main.state(indM(1),indM(2));
stageM = main.stage(indM(1),indM(2));
stateS = second.state(indS(1),indS(2));
stageS = second.stage(indS(1),indS(2));

% Extract initial values of enthalpy and entropy
h1_st = stateM.h;
h2_st = stateS.h;
s1_st = stateM.s;
s2_st = stateS.s;

% Check that fluids and pressures are the same
if strcmp(main.name,second.name)
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
[stateM] = update_state(stateM,main.handle,main.read,main.TAB,2);
stateS.h = h;
[stateS] = update_state(stateS,second.handle,second.read,second.TAB,2);

% Compute stages
% Entropy change
Ds1         = stateM.s - s1_st;
Ds2         = stateS.s - s2_st;
% Main stream
stageM.Dh   = stateM.h - h1_st;
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
main.state(indM(1),indM(2)+1) = stateM;
main.stage(indM(1),indM(2))   = stageM;
second.state(indS(1),indS(2)+1) = stateS;
second.stage(indS(1),indS(2))   = stageS;

% Increase stage counter
iM = indM(2) + 1;
iS = indS(2) + 1;
end