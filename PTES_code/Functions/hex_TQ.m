function [fluidH, fluidC, iH, iC] = hex_TQ(fluidH, indH, fluidC, indC, eff, ploss, stage_type, mode, par)
% RESOLVE HEX T-Q DIAGRAM FOR A GIVEN EFFECTIVENESS

% DESCRIPTION
% TC1 and TH2 are the cold and hot temperature inlets (known)
% TC2 and TH1 are the cold and hot temperature outlets (unknown)
%
% There is five modes of operation:
% In mode == 0, the two mass flow rates (mH and mC) must be known, and the
% parameter "par" is unused
% In mode == 1, mC (unknown) is computed from par = Crat = mH*CpH/(mC*CpC)
% In mode == 2, mH (unknown) is computed from par = Crat = mH*CpH/(mC*CpC)
% In mode == 3, TC2 is specified (TC2=par) and mC (unknown) is computed
% In mode == 4, TH1 is specified (TH1=par) and mH (unknown) is computed
%
% There is two possible algorithms available for modes 0, 1 and 2. Both
% algorithms start by finding the conditions under which DT_pinch = 0;
% 'DT_min' does so by iterating over TH1, while 'Equal_Cp' solves
% mH*CpH=mC*CpC (at the pinch point temperature). Both algorithms produce
% essentially equal outputs, although the latter can be slightly faster.

switch mode
    case {0,1,2}
        % Select solution algorithm (either 'DT_min' or 'Equal_Cp')
        algorithm = 'Equal_Cp';
end

% Set number of sections to solve algorithm
n = 100;

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
[hvH,TvH] = get_h_T(fluidH,TC1-5,TH2+5,pH2,n);
[hvC,TvC] = get_h_T(fluidC,TC1-5,TH2+5,pC1,n);

% Obtain preliminary minimum and maximum enthalpy outlets (hot outlet
% cannot be colder than cold inlet, and vice-versa)
hH1_min = rtab1(TvH,hvH,TC1,1);
hC2_max = rtab1(TvC,hvC,TH2,1);

% Compute average 'overall' specific heat capacities
CpHmean = (hH2 - hH1_min)/(TH2-TC1);
CpCmean = (hC2_max - hC1)/(TH2-TC1);

% Determine mass flow rates
switch mode
    case 0
        % Both mass flow rates previously specified
        if any([mH,mC] == 0), error('mH and mC must be known'); end
        
    case 1
        % Only mass flow rate of hot fluid previously specified, compute
        % mass flow rate of cold fluid according to Crat
        if mH == 0, error('mH must be known in mode==1'); end
        Crat = par; %Crat = mH*CpH / (mC*CpC)        
        mC = mH*CpHmean/(CpCmean*Crat);
        stateC.mdot = mC;
        
    case 2
        % Only mass flow rate of cold fluid previously specified, compute
        % mass flow rate of hot fluid according to Crat
        if mC == 0, error('mC must be known in mode==2'); end
        Crat = par; %Crat = mH*CpH / (mC*CpC)        
        mH = Crat*mC*CpCmean/CpHmean;
        stateH.mdot = mH;
        
    case 3
        % Set TC2 = par, and compute mC. Mass flow rate of hot fluid must
        % be previously specified
        if mH == 0, error('mH must be known in mode==3'); end
        if any([par<=TC1,par>=TH2]), error('par must be TC1<par<TH2'); end
        
    case 4
        % Set TH1 = par, and compute mH. Mass flow rate of cold fluid must
        % be previously specified
        if mC == 0, error('mC must be known in mode==4'); end
        if any([par<=TC1,par>=TH2]), error('par must be TC1<par<TH2'); end
end

switch mode
    case {0,1,2}
        
        % Compute preliminary QMAX (hot outlet cannot be colder than cold inlet,
        % and vice-versa) and update hH1_min accordingly
        QMAX0 = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)]);
        hH1_min = hH2 - QMAX0/mH;
        
        switch algorithm
            case 'DT_min'
                
                % Compute hH1 for which DTmin = 0, using golden search method
                f1   = @(hH1) DTmin(mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,n,'hH1',hH1);
                TolX = (hH2-hH1_min)/1e4; %tolerance
                options = optimset('TolX',TolX);%,'Display','iter');
                hH1  = fminbnd(f1,hH1_min,hH2,options);
                
                % Compute QMAX
                QMAX = mH*(hH2 - hH1);
                
            case 'Equal_Cp'
                
                % Compute QMAX using the fact that mH*CpH = mC*CpC for any intermediate
                % pinch point
                [~, CpH] = get_Cp( fluidH, TC1, TH2, pH2, n );
                [Tv,CpC] = get_Cp( fluidC, TC1, TH2, pC1, n );
                CH  = mH*CpH;
                CC  = mC*CpC;
                fx  = CH - CC;
                fx0 = fx(1);
                np  = 0; %number of pinch points found
                ip  = zeros(n-1,1); %pinch point index
                for i0=2:n-1
                    if fx(i0)*fx0<0
                        np = np + 1;
                        ip(np) = i0; %pinch point index
                        fx0 = fx(i0);
                    end
                end
                Tp = zeros(np,1);
                Qp = zeros(size(Tp));
                for i0=1:np
                    Tp(i0) = Tv(ip(i0));
                    QC = mC*(rtab1(TvC,hvC,Tp(i0),1)-hC1); % mC*DhC from TC1 to Tp
                    QH = mH*(hH2-rtab1(TvH,hvH,Tp(i0),1)); % mH*DhH from Tp to TH2
                    Qp(i0) = QC + QH;
                end
                % QMAX has to be the minimum between Qp, QC and QH:
                QMAX = min([Qp;QMAX0]);
        end
        
        % Compute actual heat transfer based on heat exchanger effectiveness
        QT  = QMAX*eff;
        
        % Determine outlet enthalpies and temperatures
        hC2 = hC1 + QT/mC;
        hH1 = hH2 - QT/mH;
        
    case 3
        
        % Set outlet conditions of cold fluid
        TC2 = par;
        hC2 = rtab1(TvC,hvC,TC2,1);
        
        % Compute preliminary QMAX (hot outlet cannot be colder than cold
        % inlet) and set boundaries accordingly
        QMAX0 = mH*(hH2 - hH1_min);
        mCmin = 0;
        mCmax = QMAX0/(hC2 - hC1);
        
        % Find value of mC for which DTmin = 0, using golden search method
        f2    = @(mC) DTmin(mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,n,'hC2',hC2);
        TolX  = (mCmax - mCmin)/1e4; %tolerance
        options = optimset('TolX',TolX);%,'Display','iter');
        mC   = fminbnd(f2,mCmin,mCmax,options);
        
        % Compute QMAX and QT (actual heat transfer based on heat exchanger
        % effectiveness)
        QMAX = mC*(hC2 - hC1);
        QT   = QMAX*eff;
        
        % Update mC and outlet conditions of hot fluid
        mC = QT/(hC2 - hC1);
        stateC.mdot = mC;
        hH1 = hH2 - QT/mH;
        %keyboard % Ok to ignore this?
        
    case 4
        
        % Set outlet conditions of hot fluid
        TH1 = par;
        hH1 = rtab1(TvH,hvH,TH1,1);
        
        % Compute preliminary QMAX (cold outlet cannot be hotter than hot
        % inlet) and set boundaries accordingly
        QMAX0 = mC*(hC2_max - hC1);
        mHmin = 0;
        mHmax = QMAX0/(hH2 - hH1);
        
        % Find value of mH for which DTmin = 0, using golden search method
        f2    = @(mH) DTmin(mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,n,'hH1',hH1);
        TolX  = (mHmax - mHmin)/1e4; %tolerance
        options = optimset('TolX',TolX);%,'Display','iter');
        mH   = fminbnd(f2,mHmin,mHmax,options);
        
        % Compute QMAX and QT (actual heat transfer based on heat exchanger
        % effectiveness)
        QMAX = mH*(hH2 - hH1);
        QT   = QMAX*eff;
        
        % Update mH and outlet conditions of cold fluid
        mH = QT/(hH2 - hH1);
        stateH.mdot = mH;
        hC2 = hC1 + QT/mC;
        
    otherwise
        error('not implemented')        
end

% % To see the temperature distribution after applying the heat exchanger
% % effectiveness, uncomment lines below:
% [~,~,TC,TH,QS] = DTmin(mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,n,'hH1',hH1);
% figure(10)
% plot(QS./QS(end),TH,'r'); hold on;
% plot(QS./QS(end),TC,'b'); hold off;
% xlabel('Cumulative heat transfer')
% ylabel('Temperature')
% legend([fluidH.name,', ',sprintf('%.1f',pH2/1e5),' bar'],...
%     [fluidC.name,', ',sprintf('%.1f',pC1/1e5),' bar'],'Location','Best')
% keyboard

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
if any(mode==[1,2,3,4])
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

function [ Tv, Cpv ] = get_Cp( fluid, T1, T2, pressure, n )
% Obtain the Tv and Cpv arrays of a given fluid for the hex subroutines.
% Data ordered in regular intervals of T. Cp as a function of T.

if strcmp(fluid.read,'CP') %read from CoolProp
    
    Tv   = linspace(T1,T2,n)'; % Temperature array between TC1 and TH2
    pv  = ones(size(Tv)).*pressure; % pressure array
    dT   = (T2-T1)/(n*2);
    h_up = CP1('PT_INPUTS',pv,Tv+dT,'H',fluid.handle);
    h_lo = CP1('PT_INPUTS',pv,Tv-dT,'H',fluid.handle);
    Cpv  = (h_up - h_lo)/(2*dT);
    
elseif strcmp(fluid.read,'TAB') %read from table
    
    T   = fluid.TAB(:,1);
    Cp  = fluid.TAB(:,5);    
    Tv  = linspace(T1,T2,n)'; % Temperature array between TC1 and TH2
    Cpv = rtab_1D_1out(T,Cp,Tv,0);
else
    error('not implemented')
end

end

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
else
    varargout{1} = solution;
    varargout{2} = DT;
    varargout{3} = TC;
    varargout{4} = TH;
    varargout{5} = QS;
end

end