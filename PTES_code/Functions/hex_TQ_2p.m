function [fluidH, fluidC, iH, iC] = hex_TQ_2p(fluidH, indH, fluidC, indC, eff, ploss, stage_type, cond, var)
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

% Set inlet pressures, enthalpies, entropies and mass flow rates
TH2 = stateH.T;
pH2 = stateH.p;
hH2 = stateH.h;
sH2 = stateH.s;
TC1 = stateC.T;
pC1 = stateC.p;
hC1 = stateC.h;
sC1 = stateC.s;
mH = stateH.mdot;
mC = stateC.mdot; 

% Obtain temperature arrays as a function of the enthalpy arrays
[hvH,TvH] = hex_prop(fluidH,TC1-1,TH2+1,pH2,n);
[hvC,TvC] = hex_prop(fluidC,TC1-1,TH2+1,pC1,n);

% Obtain preliminary minimum and maximum enthalpy outlets (hot outlet
% cannot be colder than cold inlet, and vice-versa)
hH1_min = rtab1(TvH,hvH,TC1,1);
hC2_max = rtab1(TvC,hvC,TH2,1);

% Determine mass flow rates
switch cond
    case 0
        % Both mass flow rates previously specified
        if any([mH,mC] == 0), error('mH and mC must be known if cond==0'); end
        
    case 1
        % Only mass flow rate of hot fluid previously specified, compute
        % mass flow rate of cold fluid according to Crat
        Crat = var; %Crat = mH*CpH / (mC*CpC)
        CpHmean = (hH2 - hH1_min)/(TH2-TC1);
        CpCmean = (hC2_max - hC1)/(TH2-TC1);
        
        if mH == 0, error('mH must be known if cond==1'); end
        mC = mH*CpHmean/(CpCmean*Crat);
        stateC.mdot = mC;
        
    case 2
        % Only mass flow rate of cold fluid previously specified, compute
        % mass flow rate of hot fluid according to Crat
        Crat = var; %Crat = mH*CpH / (mC*CpC)
        CpHmean = (hH2 - hH1_min)/(TH2-TC1);
        CpCmean = (hC2_max - hC1)/(TH2-TC1);
        
        if mC == 0, error('mC must be known if cond==2'); end
        mH = Crat*mC*CpCmean/CpHmean;
        stateH.mdot = mH;                    
end

% Compute preliminary QMAX (hot outlet cannot be colder than cold inlet,
% and vice-versa) and update hH1_min accordingly
QMAX = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)]);
hH1_min = hH2 - QMAX/mH;

% Compute hH1 for which DT_min = 0
f = @(hH1) DT_area(hH1,mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,n);
[hH1,~,~,~,~] = golden_search(f,hH1_min,hH2,(hH2-hH1_min)/1e3,'Min',100);
% % To see how the golden_search evolves, comment line above and uncomment
% % lines below:
% [hH1,~,xv,yv,iter] = golden_search(f,hH1_min,hH2,(hH2-hH1_min)/1e3,'Min',100);
% figure(4)
% plot(1:iter,xv(1:iter))
% xlabel('Iterations')
% ylabel('hH1')
% figure(5)
% semilogy(1:iter,yv(1:iter))
% xlabel('Iterations')
% ylabel('DT area')

% Compute QMAX
QMAX = mH*(hH2 - hH1);

% Compute actual heat transfer based on heat exchanger effectiveness
QT   = QMAX*eff;

% Determine outlet enthalpies and temperatures
hC2 = hC1 + QT/mC;
hH1 = hH2 - QT/mH;

% To see the temperature distribution after applying the heat exchanger
% effectiveness, uncomment lines below:
[~,DT,TC,TH,QS] = DT_area(hH1,mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,n);
figure(10)
plot(QS./QS(end),TH,'r'); hold on;
plot(QS./QS(end),TC,'b'); hold off;
xlabel('Cumulative heat transfer')
ylabel('Temperature')
legend([fluidH.name,', ',sprintf('%.1f',pH2/1e5),' bar'],[fluidC.name,', ',sprintf('%.1f',pC1/1e5),' bar'],'Location','Best')
figure(11)
plot(QS./QS(end),DT,'r');
xlabel('Cumulative heat transfer')
ylabel('Temperature difference')

% Update states
stateH.h = hH1;    
stateH.p = pH2*(1-ploss);
stateH   = update_state(stateH,fluidH.handle,fluidH.read,fluidH.TAB,2);
stateC.h = hC2;
stateC.p = pC1*(1-ploss);
stateC   = update_state(stateC,fluidC.handle,fluidC.read,fluidC.TAB,2);

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
if any(cond==[1,2])
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

function [ hv, Tv ] = hex_prop( fluid, T1, T2, pressure, n )
% Obtain the hv and Tv arrays of a given fluid for the hex subroutines

if strcmp(fluid.read,'CP') %read from CoolProp
    
    h1  = CP1('PT_INPUTS',pressure,T1,'H',fluid.handle);
    h2  = CP1('PT_INPUTS',pressure,T2,'H',fluid.handle);
    hv  = linspace(h1,h2,n)';       % enthalpy array between TC1 and TH2
    pv  = ones(size(hv)).*pressure; % pressure array
    Tv  = CP1('HmassP_INPUTS',hv,pv,'T',fluid.handle); % temperature array
    
elseif strcmp(fluid.read,'TAB') %read from table
    %Tv  = linspace(TC1,TH2,n)'; %temperature array between TC1 and TH2
    Tx  = fluid.TAB(:,1);
    hy  = fluid.TAB(:,2);
    h1  = rtab1(Tx,hy,T1,0);
    h2  = rtab1(Tx,hy,T2,0);
    hv  = linspace(h1,h2,n)';       % enthalpy array between TC1 and TH2
    Tv  = rtab1(hy,Tx,hv,1);
else
    error('not implemented')
end

end

function varargout = DT_area(hH1,mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,n)
% Compute the temperature difference between the two streams inside the
% heat exchanger, for a given hH1 (enthalpy at hot stream outlet)

% Compute temperature distribution of hot stream
hH  = linspace(hH1,hH2,n)';
TH  = rtab1(hvH,TvH,hH,0);

% Compute temperature distribution of cold stream
QS  = (hH - hH1)*mH; % cummulative heat transfer
hC  = hC1 + QS/mC;
TC  = rtab1(hvC,TvC,hC,0);

% Compute temperature difference between the two streams
DT  = TH - TC;

% Integrate DT to find objective value (the DT 'area' on the T-Q diagram).
% Negative DT values (for non-physically possible hH1 values) are turned
% into positive and artificially increased, to avoid finding incorrect
% solutions when minimising the objective value.
DT_abs = DT;
neg    = DT_abs<0;
DT_abs(neg) = abs(DT_abs(neg))*1e6;

solution = sum(DT_abs);

% fprintf(1,'\nhH1 = %f, sol = %f',hH1,solution)
% figure(1)
% plot(QS./QS(end),TH,'r'); hold on;
% plot(QS./QS(end),TC,'b'); hold off;
% 
% figure(2)
% plot(QS./QS(end),DT,'r'); hold off;
% 
% figure(3)
% plot(QS./QS(end),DT_abs,'r'); hold off;
% 
% keyboard

if nargout == 1
    varargout{1} = solution;
else
    varargout{1} = solution;
    varargout{2} = DT;
    varargout{3} = TC;
    varargout{4} = TH;
    varargout{5} = QS;
    figure(1)
    plot(QS./QS(end),TH,'r'); hold on;
    plot(QS./QS(end),TC,'b'); hold off;
    
    figure(2)
    plot(QS./QS(end),abs(DT),'r'); hold off;
end

end