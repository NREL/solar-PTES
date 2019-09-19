function [fluidH, fluidC, iH, iC] = hex_TQ_cond(fluidH, indH, fluidC, indC, eff, Crat, ploss, stage_type, cond, var)
% RESOLVE HEX T-Q DIAGRAM FOR A GIVEN EFFECTIVENESS

% Set number of sections for hex_core algorithm
n = 100;

% Set inlet temperatures (nomenclature: cold inlet is position 1, hot inlet is position 2)
TH2 = fluidH.state(indH(1),indH(2)).T;
TC1 = fluidC.state(indC(1),indC(2)).T;
% Check which one is fluidH and which is fluidC and swap them if necessary
if TC1 > TH2 % swap needed
    swap = 1;
    fluidH0 = fluidH;
    fluidH  = fluidC;
    fluidC  = fluidH0;
    indH0 = indH;
    indH  = indC;
    indC  = indH0;
else
    swap = 0;
end

% Import fluid.state and fluid.stage
stateH = fluidH.state(indH(1),indH(2));
stageH = fluidH.stage(indH(1),indH(2));
stateC = fluidC.state(indC(1),indC(2));
stageC = fluidC.stage(indC(1),indC(2));

% Set stage type
stageH.type = stage_type;
stageC.type = stage_type;

% Set inlet pressures, enthalpies and entropies
TH2 = stateH.T;
pH2 = stateH.p;
hH2 = stateH.h;
sH2 = stateH.s;
TC1 = stateC.T;
pC1 = stateC.p;
hC1 = stateC.h;
sC1 = stateC.s;

% Set Cp and enthalpy arrays
[TvH,CpvH,hvH] = hex_get_prop(fluidH,TC1,TH2,pH2,n);
[TvC,CpvC,hvC] = hex_get_prop(fluidC,TC1,TH2,pC1,n);

% Set mass flow rates
if cond==4 %mass flow rates previously specified in any case
    mH = stateH.mdot;
    mC = stateC.mdot;
else
    if all(strcmp({fluidH.job,fluidC.job},{'SF','WF'}))
        CpHmean = mean(CpvH); CpCmean = mean(CpvC);
        mC  = stateC.mdot;
        mH  = mC*CpCmean*Crat / CpHmean; %Crat = mH*CpH / (mC*CpC)
        stateH.mdot = mH;
    elseif all(strcmp({fluidH.job,fluidC.job},{'WF','SF'}))
        CpHmean = mean(CpvH); CpCmean = mean(CpvC);
        mH  = stateH.mdot;
        mC  = mH*CpHmean / (Crat*CpCmean); %Crat = mH*CpH / (mC*CpC)
        stateC.mdot = mC;
    elseif all(strcmp({fluidH.job,fluidC.job},{'WF','WF'}))
        mH = stateH.mdot;
        mC = stateC.mdot;
    else
        error('not implemented')
    end
end
% Run HEX core algorithm
[QT] = hex_core_q(TvH, CpvH, hvH, TvC, CpvC, hvC, mH, mC, TC1, TH2, eff);

% Determine outlet enthalpies and temperatures
hC2 = hC1 + QT/mC;
TC2 = rtab1(hvC,TvC,hC2,1);
hH1 = hH2 - QT/mH;
TH1 = rtab1(hvH,TvH,hH1,1);

% -1 chosen since default value seems to be 0
if cond==-1 % var is minimum allowed TC2
    if TC2 < var && var < TH2
        TC2 = var;
        hC2 = rtab1(TvC,hvC,TC2,0);
        mC  = QT/(hC2 - hC1);
        stateC.mdot = mC;
    end
elseif cond==1 % var is maximum allowed TC2
    if TC2 > var
        TC2 = var;
        hC2 = rtab1(TvC,hvC,TC2,0);
        mC  = QT/(hC2 - hC1);
        stateC.mdot = mC;
    end
elseif cond==2 % var is minimum allowed TH1
    if TH1 < var
        TH1 = var;
        hH1 = rtab1(TvH,hvH,TH1,0);
        mH  = QT/(hH2 - hH1);
        stateH.mdot = mH;
    end
elseif cond==3 % var is maximum allowed TH1
    if TH1 > var && var > TC1
        TH1 = var;
        hH1 = rtab1(TvH,hvH,TH1,0);
        mH  = QT/(hH2 - hH1);
        stateH.mdot = mH;
    end
end

% Update states
stateH.h = hH1;    
stateH.p = pH2*(1-ploss);
stateH = update_state(stateH,fluidH.handle,fluidH.read,fluidH.TAB,2);
stateC.h = hC2;
stateC.p = pC1*(1-ploss);
stateC = update_state(stateC,fluidC.handle,fluidC.read,fluidC.TAB,2);

% Compute stages
% Entropy change
DsH         = stateH.s - sH2;
DsC         = stateC.s - sC1;
% Hot stream
stageH.Dh   = stateH.h - hH2;
stageH.sirr = (stateH.mdot*DsH + stateC.mdot*DsC)/stateH.mdot;
stageH.q    = stageH.Dh;
stageH.w    = 0;
stageH.type = stage_type;
% Cold stream
stageC.Dh   = stateC.h - hC1;
stageC.sirr = (stateC.mdot*DsC + stateH.mdot*DsH)/stateC.mdot;
stageC.q    = stageC.Dh;
stageC.w    = 0;
stageC.type = stage_type;
if strcmp(stage_type,'regen')
    stageH.sirr=0; %to avoid counting the lost work twice
end

% Export computed states and stages back into fluids
fluidH.state(indH(1),indH(2)+1) = stateH; % Result goes into next state
fluidH.stage(indH(1),indH(2))   = stageH; % Result stays in current stage
fluidC.state(indC(1),indC(2)+1) = stateC; % Result goes into next state
fluidC.stage(indC(1),indC(2))   = stageC; % Result stays in current stage

% Update mass flow rates for inlet state, if necessary
if cond~=4
    fluidH.state(indH(1),indH(2)).mdot = mH;
    fluidC.state(indC(1),indC(2)).mdot = mC;
end

% Reverse swap, if necessary
if swap == 1
    fluidH0 = fluidH;
    fluidH  = fluidC;
    fluidC  = fluidH0;
    indH0 = indH;
    indH  = indC;
    indC  = indH0;
end

% Increase stage counter
iH = indH(2) + 1;
iC = indC(2) + 1;

end