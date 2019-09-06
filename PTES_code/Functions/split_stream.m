function [main,second,iM] = split_stream(main,indM,second,indS,fraction)

% Separate main stream into main and secondary stream
% The "stage" only happens for the main stream. The secondary stream is
% "born" from the main.

% Import fluid.state and fluid.stage for main stream
stateM = main.state(indM(1),indM(2));
stageM = main.stage(indM(1),indM(2));

% Second stream has same state as main stream
stateS = stateM;

% Update mass flow rates
mdot  = stateM.mdot;
mdotS = mdot*fraction;
stateS.mdot = mdotS;
stateM.mdot = mdot - mdotS;

%main.state(indM(1),indM(2)+1).mdot = mdot - mdotS;

% Compute stage energy flows
stageM.Dh   = 0;
stageM.type = 'split';
stageM.q    = 0;
stageM.w    = 0;
stageM.sirr = 0;

% Reassing values to fluid.state and fluid.stage
main.state(indM(1),indM(2)+1) = stateM;
main.stage(indM(1),indM(2))   = stageM;
second.state(indS(1),indS(2)) = stateS;

% Increase stage counter for main stream
iM  = indM(2) + 1;
end

