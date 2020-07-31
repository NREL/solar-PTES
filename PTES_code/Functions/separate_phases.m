function [main,second,iM] = separate_phases(main,indM,second,indS)

% Separate 2-phase stream into liquid (main) and gas (second) streams
% The "stage" only happens for the main stream. The secondary stream is
% "born" from the main.

% Import fluid.state and fluid.stage for main stream
stateM = main.state(indM(1),indM(2));
stageM = main.stage(indM(1),indM(2));
stateS = stateM;

% Find gas fraction
fraction = RPN('HmassP_INPUTS',stateM.h,stateM.p,'Q',main);
if (fraction < 0 || fraction > 1)
    error('***Fluid is outside the saturation curve, separation not possible***')
end

% Update mass flow rates
mdot  = stateM.mdot;
mdotS = mdot*fraction;
stateS.mdot = mdotS;
stateM.mdot = mdot - mdotS;

% Update states
stateM.h = RPN('PQ_INPUTS',stateM.p,0.0,'H',main);
stateM   = update_state(stateM,main,2);
stateS.h = RPN('PQ_INPUTS',stateS.p,1.0,'H',second);
stateS   = update_state(stateS,second,2);

% Compute stage energy flows
stageM.Dh   = 0;
stageM.type = 'separate';
stageM.q    = 0;
stageM.w    = 0;
stageM.sirr = 0;

% Reassing values to fluid.state and fluid.stage
main.state(indM(1),indM(2)+1) = stateM;
main.stage(indM(1),indM(2))   = stageM;
second.state(indS(1),indS(2)) = stateS;

%Increase stage counter for main stream
iM  = indM(2) + 1;
end

