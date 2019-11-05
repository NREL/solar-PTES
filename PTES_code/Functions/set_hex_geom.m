function [HX] = set_hex_geom(fluidH, indH, fluidC, indC, eff, ploss, mode, var)
% Obtain geometric parameters based on performance objectives, using
% analytical solutions. It should be expected that objectives will be met
% accurately (only) when using fluids with small variations of
% thermophysical properties.
        
% Import fluid.state and fluid.stage
stateH = fluidH.state(indH(1),indH(2));
stateC = fluidC.state(indC(1),indC(2));

% Set inlet temperatures and pressures
THin = stateH.T;
pHin = stateH.p;
TCin = stateC.T;
pCin = stateC.p;
mH  = stateH.mdot;
mC  = stateC.mdot;

% Obtain temperature arrays as a function of the enthalpy arrays
[hvH,TvH] = get_h_T(fluidH,TCin-1,THin+1,pHin,100);
[hvC,TvC] = get_h_T(fluidC,TCin-1,THin+1,pCin,100);

% Obtain minimum and maximum enthalpy outlets and determine mean Cp
hHout_min = rtab1(TvH,hvH,TCin,1);
hCout_max = rtab1(TvC,hvC,THin,1);
CpHmean = (hH2 - hHout_min)/(THin-TCin);
CpCmean = (hCout_max - hC1)/(THin-TCin);

% Determine mass flow rates
switch mode
    case 0
        % Both mass flow rates previously specified
        if any([mH,mC] == 0), error('mH and mC must be known if cond==0'); end
        
    case 1
        % Only mass flow rate of hot fluid previously specified, compute
        % mass flow rate of cold fluid according to Crat
        if mH == 0, error('mH must be known if cond==1'); end
        Crat = var; %Crat = mH*CpH / (mC*CpC)
        mC = mH*CpHmean/(CpCmean*Crat);
        
    case 2
        % Only mass flow rate of cold fluid previously specified, compute
        % mass flow rate of hot fluid according to Crat        
        if mC == 0, error('mC must be known if cond==2'); end
        Crat = var; %Crat = mH*CpH / (mC*CpC)
        mH = Crat*mC*CpCmean/CpHmean;      
end

% Set the minimum number of transfer units that each stream should have to
% obtain the specified effectiveness
Ntu_min = 2/(1-eff);

% Determine which stream is the high pressure one (which flows inside the
% tubes - indicated 1) and which is the low pressure one (which flows on
% the shell side - indicated 2)
if pCin >= pHin
    % Cold stream is high pressure stream (tube side). Hot stream is low
    % pressure stream (shell side).
    p1 = pCin;
    m1 = mC;
    m2 = mH;
else
    % Hot stream is high pressure stream (tube side). Cold stream is low
    % pressure stream (shell side).
    p1 = pHin;
    m1 = mH;
    m2 = mC;
end


warning('still in development')
keyboard



end

% % Declare the two fluid streams
% H = stream; H.mdot = mH; H.name = fluidH.name;
% C = stream; C.mdot = mC; C.name = fluidC.name;
% 
% % Extract parameters
% NX = HX.NX;
% 
% % Compute derived geometric parameters and mass fluxes
% [C, H, HX] = shell_and_tube_geom(C, H, HX);
% 
% % Compute preliminary QMAX (hot outlet cannot be colder than cold inlet,
% % and vice-versa) and update hH1_min accordingly
% QMAX = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)]);
% hH1_min = hH2 - QMAX/mH;
% 
% % Set initial conditions for iteration procedure
% % Pressures
% H.pin = pH2; C.pin = pC1;
% H.p = ones(NX+1,1)*H.pin;
% C.p = ones(NX+1,1)*C.pin;
% 
% % Compute hH1 for computed area equals C.A
% f1 = @(hH1) compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,HX);
% %plot_function(f1,hH1_min,hH2,100,31)
% TolX = (hH2-hH1_min)/1e4; %tolerance
% options = optimset('TolX',TolX,'Display','notify');
% hH1  = fminbnd(f1,hH1_min,hH2,options);
% 
% % Obtain output parameters for converged solution
% [~,C,H,~,~] = compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,HX);
% 
% % Update states
% stateH.h = H.h(1);    
% stateH.p = H.p(1);
% stateH   = update_state(stateH,fluidH.handle,fluidH.read,fluidH.TAB,2);
% stateC.h = C.h(NX+1);
% stateC.p = C.p(NX+1);
% stateC   = update_state(stateC,fluidC.handle,fluidC.read,fluidC.TAB,2);
% 
% % Compute stages
% % Entropy change
% DsH         = stateH.s - sH2;
% DsC         = stateC.s - sC1;
% % Hot stream
% stageH.Dh   = stateH.h - hH2;
% stageH.sirr = (stateH.mdot*DsH + stateC.mdot*DsC)/stateH.mdot;
% stageH.q    = stageH.Dh;
% stageH.w    = 0;
% stageH.type = stage_type;
% % Cold stream
% stageC.Dh   = stateC.h - hC1;
% stageC.sirr = (stateC.mdot*DsC + stateH.mdot*DsH)/stateC.mdot;
% stageC.q    = stageC.Dh;
% stageC.w    = 0;
% stageC.type = stage_type;
% if strcmp(stage_type,'regen')
%     stageH.sirr=0; %to avoid counting the lost work twice
% end
% 
% % Save data for plots (use with plot_hex_TQA function)
% [~,C,H,QS,AS] = compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,HX);
% HX.C  = C;
% HX.H  = H;
% HX.QS = QS;
% HX.AS = AS;
% 
% % Export computed states and stages back into fluids
% fluidH.state(indH(1),indH(2)+1) = stateH; % Result goes into next state
% fluidH.stage(indH(1),indH(2))   = stageH; % Result stays in current stage
% fluidC.state(indC(1),indC(2)+1) = stateC; % Result goes into next state
% fluidC.stage(indC(1),indC(2))   = stageC; % Result stays in current stage
% 
% % Update mass flow rates for inlet state, if necessary
% if any(mode==[1,2])
%     fluidH.state(indH(1),indH(2)).mdot = mH;
%     fluidC.state(indC(1),indC(2)).mdot = mC;
% end
% 
% % Reverse swap, if necessary
% if swap == 1
%     fluidH0 = fluidH;
%     fluidH  = fluidC;
%     fluidC  = fluidH0;
%     indH0 = indH;
%     indH  = indC;
%     indC  = indH0;
% end
% 
% % Increase stage counter
% iH = indH(2) + 1;
% iC = indC(2) + 1;
% 
% 
% 
% end
% 
% 
% 
% %%% SUPPORT FUNCTIONS %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [ hv, Tv ] = get_h_T( fluid, T1, T2, pressure, n )
% % Obtain the hv and Tv arrays of a given fluid for the hex subroutines.
% % Data ordered in regular intervals of h. T as a function of h.
% 
% if strcmp(fluid.read,'CP') %read from CoolProp
%     
%     h1  = CP1('PT_INPUTS',pressure,T1,'H',fluid.handle);
%     h2  = CP1('PT_INPUTS',pressure,T2,'H',fluid.handle);
%     hv  = linspace(h1,h2,n)';       % enthalpy array between TC1 and TH2
%     pv  = ones(size(hv)).*pressure; % pressure array
%     Tv  = CP1('HmassP_INPUTS',hv,pv,'T',fluid.handle); % temperature array
%     
% elseif strcmp(fluid.read,'TAB') %read from table
%     
%     Tx  = fluid.TAB(:,1);
%     hy  = fluid.TAB(:,2);
%     h1  = rtab1(Tx,hy,T1,0);
%     h2  = rtab1(Tx,hy,T2,0);
%     hv  = linspace(h1,h2,n)'; % enthalpy array between TC1 and TH2
%     Tv  = rtab1(hy,Tx,hv,1);  % temperature array between TC1 and TH2
% else
%     error('not implemented')
% end
% 
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function varargout = compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,HX)
% % For a given hot fluit outlet enthalpy (hH1), compute the TQ diagram, the
% % properties of the fluids at each point and the corresponding heat
% % transfer area of the heat exchanger. Compare that to the reference heat
% % transfer area and return the difference between the two.
% 
% % Extract parameters
% NX = HX.NX;
% 
% % Compute enthalpy arrays from hH1 (outlet guess value) and hH2 and hC1
% % (fixed inlet values)
% H.h = linspace(hH1,hH2,NX+1)';
% QS  = mH*(H.h - hH1);
% C.h = hC1 + QS/mC;
% 
% % Create array to check convergence. First element is computed heat
% % transfer area. Later come the pressure points along each stream
% CON_0 = [0; H.p; C.p]; % initial value
% NI    = 10;
% RES   = zeros(1,NI); % residuals
% TOL   = 1e-4;
% 
% for iI = 1:NI
%     
%     % UPDATE PROPERTIES
%     % Cold stream
%     C = stream_update(fluidC,C,1);
%     % Hot stream
%     H = stream_update(fluidH,H,1);
%     
%     % COMPUTE HEAT TRANSFER COEFFICIENTS
%     % Cold stream
%     C.Re = C.D*C.G./C.mu;
%     [C.Cf,C.St] = developed_flow(C.Re,C.Pr,HX.shape);
%     C.ht  = C.G*C.Cp.*C.St;
%     % Hot stream
%     H.Re = H.D*H.G./H.mu;
%     [H.Cf,H.St] = developed_flow(H.Re,H.Pr,HX.shape);
%     H.ht  = H.G*H.Cp.*H.St;
%     % Overall heat transfer coefficient (based on cold side heat transfer area).
%     % Neglects wall thermal resistance and axial conduction.
%     UlC  = 1./(C.A./(H.A*H.ht) + 1./C.ht);
%     
%     %COMPUTE AVERAGED ARRAYS
%     H.T_AV = 0.5*(H.T(1:NX) + H.T(2:NX+1));
%     C.T_AV = 0.5*(C.T(1:NX) + C.T(2:NX+1));
%     UlC_AV = 0.5*(UlC(1:NX) + UlC(2:NX+1));
%     DT_AV  = H.T_AV - C.T_AV;
%     
%     % COMPUTE HEAT TRANSFER AREA (cold side)
%     dQ  = (H.h(2:NX+1) - H.h(1:NX))*mH;
%     dAC = dQ./(UlC_AV.*DT_AV);
%     AC  = sum(dAC);
%     
%     % COMPUTE PRESSURE PROFILES
%     % Create averaged arrays of Cf and rho
%     Cf_H  = 0.5*(H.Cf(1:NX)  + H.Cf(2:NX+1));
%     rho_H = 0.5*(H.rho(1:NX) + H.rho(2:NX+1));
%     Cf_C  = 0.5*(C.Cf(1:NX)  + C.Cf(2:NX+1));
%     rho_C = 0.5*(C.rho(1:NX) + C.rho(2:NX+1));
%     % Obtain dL from dAC and AC
%     dL = dAC/AC;
%     % Compute arrays of pressure loss
%     Dp_H = - 2*H.G^2*Cf_H.*dL./(rho_H*H.D);
%     Dp_C = - 2*C.G^2*Cf_C.*dL./(rho_C*C.D);
%     % Update pressure profiles
%     for i=NX+1:-1:2
%         H.p(i-1) = H.p(i) + Dp_H(i-1);
%     end
%     for i=1:NX
%         C.p(i+1) = C.p(i) + Dp_C(i);
%     end
%     % Artificially avoid pressures below 20% of p_in and set error flag if
%     % needed
%     cond1 = C.p < 0.2*C.pin;
%     cond2 = H.p < 0.2*H.pin;
%     C.p(cond1) = 0.2*C.pin;
%     H.p(cond2) = 0.2*H.pin;
%     if any(cond1),warning('DpC exceeds 20%!');end
%     if any(cond2),warning('DpH exceeds 20%!');end
%     
%     % Update convergence array
%     CON = [AC; H.p; C.p]; % initial value
%     
%     % Compute residual
%     RES(iI) = max(abs((CON - CON_0)./CON));
%     %fprintf(1,'\n iteration = %d, RES = %.6f',iI,RES(iI));
%     
%     if (RES(iI)>TOL)
%         CON_0 = CON;        
%     else
%         %fprintf(1,'\n\n*** Successful convergence after %d iterations***\n',iI);
%         break
%     end
%     
% end
% if all([iI>=NI,RES(iI)>TOL])
%     error('Convergence not reached after %d iterations***\n',iI);
% end
% 
% solution = C.A - AC;
% limit    = 1.5*C.A;
% % Control physically impossible solutions
% if any([DT_AV <= 0; abs(solution) > limit; solution < 0])
%     solution = limit;
% end
% 
% %fprintf(1,'HEX AC = %5.1f, computed AC = %5.1f, solution = %5.1f\n',C.A,AC,solution);
% 
% % figure(10)
% % plot(QS./QS(end),H.T,'r'); hold on;
% % plot(QS./QS(end),C.T,'b'); hold off;
% % xlabel('Cumulative heat transfer')
% % ylabel('Temperature')
% % legend([fluidH.name,', ',sprintf('%.1f',pH2/1e5),' bar'],[fluidC.name,', ',sprintf('%.1f',pC1/1e5),' bar'],'Location','Best')
% % figure(11)
% % plot(QS./QS(end),DT,'r');
% % xlabel('Cumulative heat transfer')
% % ylabel('Temperature difference')
% % keyboard
% 
% 
% if nargout == 1
%     varargout{1} = solution;
% else
%     % Compute cumulative heat transfer area (cold side)
%     AS = zeros(size(QS));
%     for i=1:(length(AS)-1)
%         AS(i+1) = AS(i) + dAC(i);
%     end
%     varargout{1} = solution;
%     varargout{2} = C;
%     varargout{3} = H;
%     varargout{4} = QS;
%     varargout{5} = AS;
% end
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%