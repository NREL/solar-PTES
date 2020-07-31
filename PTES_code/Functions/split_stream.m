function [main,iM] = split_stream(main,iL,iM,iS,fraction)

% Separate main stream into main and secondary stream
% The "stage" only happens for the main stream. The secondary stream is
% "born" from the main.

% Import fluid.state and fluid.stage for main stream
stateM = main.state(iL,iM);
stageM = main.stage(iL,iM);

% Second stream has same state as main stream
stateS = stateM;

% Update mass flow rates
mdot  = stateM.mdot;
mdotS = mdot*fraction;
stateS.mdot = mdotS;
stateM.mdot = mdot - mdotS;

% Compute stage energy flows
stageM.Dh   = 0;
stageM.type = 'split';
stageM.q    = 0;
stageM.w    = 0;
stageM.sirr = 0;

% Reassing values to fluid.state and fluid.stage
main.state(iL,iM+1) = stateM;
main.stage(iL,iM)   = stageM;
main.state(iL,iS)   = stateS;

% Increase stage counter for main stream
iM  = iM + 1;
end

