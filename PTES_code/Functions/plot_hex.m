function [ ] = plot_hex( fluidH, indH, fluidC, indC, n, fignum )
% Makes a TQ diagram of a heat exchanger

% Set inlet temperatures (nomenclature: cold inlet is position 1, hot inlet is position 2)
TH2 = fluidH.state(indH(1),indH(2)).T;
TC1 = fluidC.state(indC(1),indC(2)).T;

% Check which one is fluidH and which is fluidC and swap them if necessary
if TC1 > TH2 % swap needed
    fluidH0 = fluidH;
    fluidH  = fluidC;
    fluidC  = fluidH0;
    indH0 = indH;
    indH  = indC;
    indC  = indH0;
end

% Import fluid.state and fluid.stage
stateH = fluidH.state(indH(1),indH(2));
stageH = fluidH.stage(indH(1),indH(2));
stateC = fluidC.state(indC(1),indC(2));
stageC = fluidC.stage(indC(1),indC(2));

% Set inlet pressures, enthalpies and entropies
TH2 = stateH.T;
pH2 = stateH.p;
hH2 = stateH.h;
mH  = stateH.mdot;
TC1 = stateC.T;
pC1 = stateC.p;
hC1 = stateC.h;
mC  = stateC.mdot;

% Set Cp and enthalpy arrays
[TvH,~,hvH] = hex_get_prop(fluidH,TC1,TH2,pH2,n);
[TvC,~,hvC] = hex_get_prop(fluidC,TC1,TH2,pC1,n);

% Obtain total heat transfer
QTH = stageH.Dh*mH;
QT  = stageC.Dh*mC;
if abs((QTH+QT)/QT) > 1e-6
    error('incorrect energy balance')
end

% Pre-allocate arrays
TC = zeros(1,n);
TH = zeros(1,n);
hC = zeros(1,n);
hH = zeros(1,n);
QS = zeros(1,n);
dQ = QT/(n-1);

% Compute arrays
TC(1) = TC1;
TH(n) = TH2;
hC(1) = hC1;
hH(n) = hH2;
QS(1) = 0;
for i=1:n-1
    hC(i+1) = hC(i) + dQ/mC;
    TC(i+1) = rtab1(hvC,TvC,hC(i+1),1);
    QS(i+1) = QS(i) + dQ;
end
for i=n:-1:2
    hH(i-1) = hH(i) - dQ/mH;
    TH(i-1) = rtab1(hvH,TvH,hH(i-1),1);
end

% Set legends
LC = sprintf('%s, %.1f bar',fluidC.name,pC1/1e5);
LH = sprintf('%s, %.1f bar',fluidH.name,pH2/1e5);

% Make plots
figure(fignum)
if mod(fignum,2)==0
    plot(QS./QS(n),TC-273.15,QS./QS(n),TH-273.15)
    legend(LC,LH,'Location','SouthEast')
else
    plot(QS./QS(n),TH-273.15,QS./QS(n),TC-273.15)
    legend(LH,LC,'Location','SouthEast')
end
ylabel('Temperature [$^\circ$C]')
%ylabel('Temperature [K]')
xlabel('Normalised heat transfer')
set(gcf,'color','w')

end

