function [fluidH, fluidC, iH, iC, HX, varargout] = hex_TQA(fluidH, indH, fluidC, indC, HX, stage_type, mode, var)
% RESOLVE HEX T-Q DIAGRAM FOR A GIVEN EFFECTIVENESS

% Set number of sections for hex_core algorithm
NX = HX.NX;

% Set inlet temperatures (nomenclature: cold inlet is position 1, hot inlet
% is position 2)
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
[hvH,TvH] = get_h_T(fluidH,TC1-1,TH2+1,pH2,NX);
[hvC,TvC] = get_h_T(fluidC,TC1-1,TH2+1,pC1,NX);

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

% Declare the two fluid streams
H = stream; H.mdot = mH; H.name = fluidH.name; H.pin = pH2;
C = stream; C.mdot = mC; C.name = fluidC.name; C.pin = pC1;

% Compute derived geometric parameters and mass fluxes
[C, H, HX] = shell_and_tube_geom(C, H, HX);

% Compute preliminary QMAX (hot outlet cannot be colder than cold inlet,
% and vice-versa) and update hH1_min accordingly
QMAX = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)]);
hH1_min = hH2 - QMAX/mH;

% Set initial conditions for iteration procedure
% Pressures
H.p = ones(NX+1,1)*H.pin;
C.p = ones(NX+1,1)*C.pin;

% Compute hH1 for computed area equals C.A
f1 = @(hH1) compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,HX);
%plot_function(f1,hH1_min,hH2,100,31); keyboard;
TolX = (hH2-hH1_min)/1e4; %tolerance
options = optimset('TolX',TolX,'Display','notify');
hH1  = fminbnd(f1,hH1_min,hH2,options);

% Obtain output parameters for converged solution
[~,C,H,QS,AS] = compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,HX);

% Update states
stateH.h = H.h(1);    
stateH.p = H.p(1);
stateH   = update_state(stateH,fluidH.handle,fluidH.read,fluidH.TAB,2);
stateC.h = C.h(NX+1);
stateC.p = C.p(NX+1);
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

% Compute QMAX. Find the value of hH1 for which DTmin = 0, using golden
% search method
f1   = @(hH1) DTmin(mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,length(hvH),'hH1',hH1);
TolX = (hH2-hH1_min)/1e4; %tolerance
options = optimset('TolX',TolX);%,'Display','iter');
hH1_0  = fminbnd(f1,hH1_min,hH2,options);
QMAX = mH*(hH2 - hH1_0);

% Save data for plots (use with plot_hex_TQA function)
HX.C  = C;
HX.H  = H;
HX.QS = QS;
HX.AS = AS;
HX.QMAX = QMAX;

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

if nargout == 5
    varargout = {};
elseif nargout == 6
    % Additionally return the effectiveness and pressure drops of each
    % channel
    eff = QS(NX+1)/QMAX;
    ploss1 = (C.pin - C.p(NX+1))/C.pin;
    ploss2 = (H.pin - H.p(1))/H.pin;
    varargout{1} = [1-eff,ploss1,ploss2];
end

end



%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ hv, Tv ] = get_h_T( fluid, T1, T2, pressure, n )
% Obtain the hv and Tv arrays of a given fluid for the hex subroutines.
% Data ordered in regular intervals of h. T as a function of h.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,HX)
% For a given hot fluit outlet enthalpy (hH1), compute the TQ diagram, the
% properties of the fluids at each point and the corresponding heat
% transfer area of the heat exchanger. Compare that to the reference heat
% transfer area and return the difference between the two.

% Extract parameters
NX = HX.NX;

% Compute enthalpy arrays from hH1 (outlet guess value) and hH2 and hC1
% (fixed inlet values)
H.h = linspace(hH1,hH2,NX+1)';
QS  = mH*(H.h - hH1);
C.h = hC1 + QS/mC;

% Create array to check convergence. First element is computed heat
% transfer area. Later come the pressure points along each stream
CON_0 = [0; H.p; C.p]; % initial value
NI    = 50;
RES   = zeros(1,NI); % residuals
TOL   = 1e-3;
impossible = false; %indicades impossible situation
for iI = 1:NI
    
    % UPDATE PROPERTIES
    % Cold stream
    C = stream_update(fluidC,C,1);
    % Hot stream
    H = stream_update(fluidH,H,1);
    
    % COMPUTE AVERAGED TEMPERATURE ARRAYS
    H.T_AV = 0.5*(H.T(1:NX) + H.T(2:NX+1));
    C.T_AV = 0.5*(C.T(1:NX) + C.T(2:NX+1));
    DT_AV  = H.T_AV - C.T_AV;
    
    % Break loop if H.T < C.T at any point
    if any(H.T <= C.T)
        impossible = true;
        AC = Inf;
        break
    end
    
    % COMPUTE HEAT TRANSFER COEFFICIENTS
    % Cold stream
    C.Re = C.D*C.G./C.mu;
    [C.Cf,C.St] = developed_flow(C.Re,C.Pr,HX.shape);
    C.ht  = C.G*C.Cp.*C.St;
    % Hot stream
    H.Re = H.D*H.G./H.mu;
    [H.Cf,H.St] = developed_flow(H.Re,H.Pr,HX.shape);
    H.ht  = H.G*H.Cp.*H.St;
    % Overall heat transfer coefficient (based on cold side heat transfer area).
    % Neglects wall thermal resistance and axial conduction.
    UlC  = 1./(C.A./(H.A*H.ht) + 1./C.ht);
    UlC_AV = 0.5*(UlC(1:NX) + UlC(2:NX+1));
    
    % COMPUTE HEAT TRANSFER AREA (cold side)
    dQ  = (H.h(2:NX+1) - H.h(1:NX))*mH;
    dAC = dQ./(UlC_AV.*DT_AV);
    AC  = sum(dAC);
    
    % COMPUTE PRESSURE PROFILES
    % Create averaged arrays of Cf and v
    Cf_H = 0.5*(H.Cf(1:NX) + H.Cf(2:NX+1));
    v_H  = 0.5*(H.v(1:NX)  + H.v(2:NX+1));
    Cf_C = 0.5*(C.Cf(1:NX) + C.Cf(2:NX+1));
    v_C  = 0.5*(C.v(1:NX)  + C.v(2:NX+1));
    % Obtain dL from dAC and AC
    dL = dAC/AC*HX.L;
    % Compute arrays of pressure loss
    Dp_H = - 2*H.G^2*Cf_H.*v_H.*dL./H.D;
    Dp_C = - 2*C.G^2*Cf_C.*v_C.*dL./C.D;
    % Update pressure profiles
    for i=NX+1:-1:2
        H.p(i-1) = H.p(i) + Dp_H(i-1);
    end
    for i=1:NX
        C.p(i+1) = C.p(i) + Dp_C(i);
    end
    % Artificially avoid pressures below 20% of p_in and set error flag if
    % needed
    cond1 = C.p < 0.2*C.pin;
    cond2 = H.p < 0.2*H.pin;
    C.p(cond1) = 0.2*C.pin;
    H.p(cond2) = 0.2*H.pin;
    %     if any(cond1)
    %         warning('DpC exceeds 20%!');
    %     end
    %     if any(cond2)
    %         warning('DpH exceeds 20%!');
    %     end
    
    % Update convergence array
    CON = [AC; H.p; C.p]; % initial value
    
    %     % Make plots (uncomment to manually check iteration procedure)
    %     figure(10)
    %     plot(QS./QS(end),H.T,'r'); hold on;
    %     plot(QS./QS(end),C.T,'b'); hold off;
    %     xlabel('Cumulative heat transfer')
    %     ylabel('Temperature')
    %     legend([fluidH.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],[fluidC.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')
    %     figure(11)
    %     plot(QS./QS(end),H.p/H.pin,'r-'); hold on
    %     plot(QS./QS(end),C.p/C.pin,'b-'); hold off
    %     ylim([0.99 1])
    %     xlabel('Cummulative heat transfer')
    %     ylabel('Relative pressure, p/p0')
    %     legend([fluidH.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],[fluidC.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')
    %     keyboard
    
    % Compute residual
    RES(iI) = max(abs((CON - CON_0)./CON));
    %fprintf(1,'\n iteration = %d, RES = %.6f',iI,RES(iI));
    
    if (RES(iI)>TOL)
        CON_0 = CON;        
    else
        %fprintf(1,'\n\n*** Successful convergence after %d iterations***\n',iI);
        break
    end
    
end
if all([iI>=NI,RES(iI)>TOL,~impossible])
    figure()
    semilogy(1:iI,RES(1:iI))
    xlabel('Iteration')
    ylabel('Convergence residual')
    error('Convergence not reached after %d iterations***\n',iI);
end

solution = C.A - AC;
limit    = 1.5*C.A;
% Control physically impossible solutions
if any([DT_AV <= 0; abs(solution) >= limit; solution < 0])
    solution = limit;
end

%fprintf(1,'HEX AC = %5.1f, computed AC = %5.1f, solution = %5.1f\n',C.A,AC,solution);


if nargout == 1
    varargout{1} = solution;
else
    % Compute cumulative heat transfer area (cold side)
    AS = zeros(size(QS));
    for i=1:(length(AS)-1)
        AS(i+1) = AS(i) + dAC(i);
    end
    varargout{1} = solution;
    varargout{2} = C;
    varargout{3} = H;
    varargout{4} = QS;
    varargout{5} = AS;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = DTmin(mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,n,mode,hout)
% Compute the temperature difference between the two streams inside the
% heat exchanger. In mode=0, hH1=hout. In mode=1, hC2=hout.

switch mode
    case 'hH1'
        
        % Set outlet enthalpy
        hH1 = hout;
        
        % Compute temperature distribution of hot stream
        hH  = linspace(hH1,hH2,n)';
        TH  = rtab1(hvH,TvH,hH,0);
        
        % Compute temperature distribution of cold stream
        QS  = (hH - hH1)*mH; % cummulative heat transfer
        hC  = hC1 + QS/mC;
        TC  = rtab1(hvC,TvC,hC,0);
        
    case 'hC2'
        
        % Set outlet enthalpy
        hC2 = hout;
        
        % Compute temperature distribution of cold stream
        hC  = linspace(hC1,hC2,n)';
        TC  = rtab1(hvC,TvC,hC,0);
        
        % Compute temperature distribution of hot stream
        QS  = (hC - hC1)*mC; % cummulative heat transfer
        hH1 = hH2 - QS(n)/mH;
        hH  = hH1 + QS/mH;
        TH  = rtab1(hvH,TvH,hH,0);
        
    otherwise
        error('not implemented')
end

% Compute temperature difference between the two streams
DT  = TH - TC;

solution = min(DT);
if solution < 0.0
    DTmax    = TH(n) - TC(1);
    solution = DTmax;
end

% % To visualise the temperature distribution every time the function is
% % called, uncomment the lines below
% figure(10)
% plot(QS./QS(end),TH,'r'); hold on;
% plot(QS./QS(end),TC,'b'); hold off;
% xlabel('Cumulative heat transfer')
% ylabel('Temperature')
% keyboard

if nargout == 1
    varargout{1} = solution;
elseif nargout == 5
    varargout{1} = solution;
    varargout{2} = DT;
    varargout{3} = TC;
    varargout{4} = TH;
    varargout{5} = QS;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%