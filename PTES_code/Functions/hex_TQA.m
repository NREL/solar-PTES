function [fluidH, fluidC, iH, iC] = hex_TQA(fluidH, indH, fluidC, indC, HX, stage_type, mode, var)
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
switch mode
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

% Declare the two fluid streams and the HEX solid material
H  = stream;
C  = stream;
M  = solid;

% Adapt nomenclature of inputs
M.shape   = HX.shape;
M.V       = HX.Vm;
M.L       = HX.L;
C.D       = HX.D_2;
M.t_D_min = HX.t_D_min;
M.sigma   = HX.sigma;
M.VR      = HX.VR;
NX        = HX.NX;
% NI        = HX.NI;
% TOL       = HX.TOL;

% Compute derived geometric parameters and mass fluxes
[ H.D, C.G, H.G, C.Ax, H.Ax, M.Ax, C.A, H.A] = printed_circuit_geom(M,C.D,C.pin,H.pin,mC,mH);

% Compute preliminary QMAX (hot outlet cannot be colder than cold inlet,
% and vice-versa) and update hH1_min accordingly
QMAX = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)]);
hH1_min = hH2 - QMAX/mH;

% INITIAL GUESS
% Set initial conditions for iteration procedure
% Pressures
H.p = ones(NX+1,1)*pH2;
C.p = ones(NX+1,1)*pC1;

warning('pressure variation not implemented')
warning('printed_circuit_geom function requires upgrade')
warning('developed flow function requires upgrade')
warning('consider finding QMAX using separate, previous iteration loop? (Cp system?)')
% In golden_search routine, the objective function should have only one
% minimum/maximum, that's why doing abs(AC - C.A) seems necessary. I could
% also try using a built-in optimiser from Matlab instead, but I am not
% sure if that would be faster.


% Compute hH1 for computed area equals C.A
f1 = @(hH1) abs(compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,pH2,pC1,NX));
% [hH1,~,~,~,~] = golden_search(f,hH1_min,hH2,(hH2-hH1_min)/1e3,'Min',100);
% To see how the golden_search evolves, comment line above and uncomment
% lines below:
[hH1,~,xv,yv,iter] = golden_search(f1,hH1_min,hH2,(hH2-hH1_min)/1e6,'Min',100);

figure(4)
plot(1:iter,xv(1:iter))
xlabel('Iterations')
ylabel('hH1')
figure(5)
plot(1:iter,yv(1:iter),'ro')
xlabel('Iterations')
ylabel('Objective function')
figure(6)
plot(xv(1:iter),yv(1:iter),'o')
xlabel('hH1')
ylabel('Objective function')
% options = optimset('Display','final','PlotFcns','optimplotfval','TolX',(hH2-hH1_min)/1e6);
% [x,fval,exitflag,output] = fminbnd(f1,hH1_min,hH2,options);

keyboard

% Compute QT and hC2
QT  = mH*(hH2 - hH1);
hC2 = hC1 + QT/mC;

% To see the temperature distribution after applying the heat exchanger
% effectiveness, uncomment lines below:
[~,C,H,DT,QS] = compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,pH2,pC1,NX);

figure(10)
plot(QS./QS(end),H.T,'r'); hold on;
plot(QS./QS(end),C.T,'b'); hold off;
xlabel('Cumulative heat transfer')
ylabel('Temperature')
legend([fluidH.name,', ',sprintf('%.1f',pH2/1e5),' bar'],[fluidC.name,', ',sprintf('%.1f',pC1/1e5),' bar'],'Location','Best')
figure(11)
plot(QS./QS(end),DT,'r');
xlabel('Cumulative heat transfer')
ylabel('Temperature difference')

keyboard
error('not implemented')

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
if any(mode==[1,2])
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
% Data ordered in regular intervals of h.

if strcmp(fluid.read,'CP') %read from CoolProp
    
    h1  = CP1('PT_INPUTS',pressure,T1,'H',fluid.handle);
    h2  = CP1('PT_INPUTS',pressure,T2,'H',fluid.handle);
    hv  = linspace(h1,h2,n)';       % enthalpy array between TC1 and TH2
    pv  = ones(size(hv)).*pressure; % pressure array
    Tv  = CP1('HmassP_INPUTS',hv,pv,'T',fluid.handle); % temperature array
    
elseif strcmp(fluid.read,'TAB') %read from table
    
    Tx  = fluid.TAB(:,1);
    hy  = fluid.TAB(:,2);
    h1  = rtab1(Tx,hy,T1,0);
    h2  = rtab1(Tx,hy,T2,0);
    hv  = linspace(h1,h2,n)'; % enthalpy array between TC1 and TH2
    Tv  = rtab1(hy,Tx,hv,1);  % temperature array between TC1 and TH2
    
else
    error('not implemented')
end

end

function varargout = compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,pH2,pC1,NX)

% Compute enthalpy arrays from hH1 (outlet guess value) and hH2 and hC1
% (fixed inlet values)
H.h = linspace(hH1,hH2,NX+1)';
QS  = mH*(H.h - hH1);
C.h = hC1 + QS/mC;

% UPDATE PROPERTIES
% Cold stream
C = stream_update(fluidC,C,1);
% Hot stream
H = stream_update(fluidH,H,1);

% COMPUTE HEAT TRANSFER COEFFICIENTS
% Cold stream
C.Re = C.D*C.G./C.mu;
[C.Cf,C.St] = developed_flow(C.Re,C.Pr,'circular');
C.ht  = C.G*C.Cp.*C.St;
% Hot stream
H.Re = H.D*H.G./H.mu;
[H.Cf,H.St] = developed_flow(H.Re,H.Pr,'circular');
H.ht  = H.G*H.Cp.*H.St;
% Overall heat transfer coefficient (based on cold side heat transfer area).
% Neglects wall thermal resistance and axial conduction.
UlC  = 1./(C.A./(H.A*H.ht) + 1./C.ht);

%COMPUTE AVERAGED ARRAYS
H.T_AV = 0.5*(H.T(1:NX) + H.T(2:NX+1));
C.T_AV = 0.5*(C.T(1:NX) + C.T(2:NX+1));
UlC_AV = 0.5*(UlC(1:NX) + UlC(2:NX+1));
DT     = H.T - C.T;
DT_AV  = H.T_AV - C.T_AV;

% COMPUTE HEAT TRANSFER AREA (cold side)
dQ  = (H.h(2:NX+1) - H.h(1:NX))*mH;
dAC = dQ./(UlC_AV.*DT_AV);
AC  = sum(dAC);

% Next steps:
% (1) Compute pressure losses along each stream (starting from each inlet)
% (2) Update properties using new pressure array (perhaps in the following
% loop?)

solution = C.A - AC;
% Control physically impossible solutions
if any([DT_AV <= 0; abs(solution) > 1.5*C.A; isnan(solution)])
    solution = 1.5*C.A;
end

fprintf(1,'HEX AC = %.1f, computed AC = %.1f\n',C.A,AC);

% figure(10)
% plot(QS./QS(end),H.T,'r'); hold on;
% plot(QS./QS(end),C.T,'b'); hold off;
% xlabel('Cumulative heat transfer')
% ylabel('Temperature')
% legend([fluidH.name,', ',sprintf('%.1f',pH2/1e5),' bar'],[fluidC.name,', ',sprintf('%.1f',pC1/1e5),' bar'],'Location','Best')
% figure(11)
% plot(QS./QS(end),DT,'r');
% xlabel('Cumulative heat transfer')
% ylabel('Temperature difference')


if nargout == 1
    varargout{1} = solution;
else
    varargout{1} = solution;
    varargout{2} = C;
    varargout{3} = H;
    varargout{4} = DT;
    varargout{5} = QS;
    %     figure(10)
    %     plot(QS./QS(end),H.T,'r'); hold on;
    %     plot(QS./QS(end),C.T,'b'); hold off;
    %     xlabel('Cumulative heat transfer')
    %     ylabel('Temperature')
    %     legend([fluidH.name,', ',sprintf('%.1f',pH2/1e5),' bar'],[fluidC.name,', ',sprintf('%.1f',pC1/1e5),' bar'],'Location','Best')
    %     figure(11)
    %     plot(QS./QS(end),DT,'r');
    %     xlabel('Cumulative heat transfer')
    %     ylabel('Temperature difference')
end


end