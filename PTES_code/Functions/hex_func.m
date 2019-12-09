function [HX, fluidH, fluidC, iH, iC] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par)
% COMPUTE HEAT EXCHANGER OUTLET CONDITIONS
%   Description
%   TC1 and TH2 are the cold and hot temperature inlets (known)
%   TC2 and TH1 are the cold and hot temperature outlets (unknown)
%   
%   There is five modes of operation:
%   In mode == 0, the two mass flow rates (mH and mC) must be known, and
%   the parameter "par" is unused
%   In mode == 1, mC (unknown) is computed from par = Crat = mH*CpH/(mC*CpC)
%   In mode == 2, mH (unknown) is computed from par = Crat = mH*CpH/(mC*CpC)
%   In mode == 3, TC2 is specified (TC2=par) and mC (unknown) is computed
%   In mode == 4, TH1 is specified (TH1=par) and mH (unknown) is computed
%   In mode == 5, TC2 is specified (TC2=par) and mH (unknown) is computed
%   In mode == 6, TH1 is specified (TH1=par) and mC (unknown) is computed
%   
%   The HX object/structure can operate according to three different models
%   If model = 'eff', the heat exchanger effectiveness and pressure loss
%   are specified
%   If model = 'UA', the overall heat transfer coefficient and pressure
%   loss are specified
%   If model = 'geom, the heat exchanger geometry is specified


% Extract parameters from HX structure according to selected model
model = HX.model;
switch model
    case 'eff'
        eff   = HX.eff;
        ploss = HX.ploss;
        stage_type = HX.stage_type;
        NX  = HX.NX;
        
    case 'UA'
        UA_ref = HX.UA;
        ploss = HX.ploss;
        stage_type = HX.stage_type;
        NX  = HX.NX;
        
    case 'geom'
        stage_type = HX.stage_type;
        NX  = HX.NX;
        
    otherwise
        error('Invalid heat exchanger model')
end

% Set inlet temperatures (nomenclature: cold inlet is position 1, hot inlet
% is position 2)
TH2 = fluidH.state(iL,iH).T;
TC1 = fluidC.state(iL,iC).T;

% Check which one is fluidH and which is fluidC and swap them if necessary
if TC1 > TH2 % swap needed
    swap = 1;
    fluidH0 = fluidH;
    fluidH  = fluidC;
    fluidC  = fluidH0;
    iH0 = iH;
    iH  = iC;
    iC  = iH0;
else
    swap = 0;
end

% Import fluid.state and fluid.stage
stateH = fluidH.state(iL,iH);
stageH = fluidH.stage(iL,iH);
stateC = fluidC.state(iL,iC);
stageC = fluidC.stage(iL,iC);

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

% Declare the two fluid streams
H = stream; H.name = fluidH.name; H.pin = pH2;
C = stream; C.name = fluidC.name; C.pin = pC1;

% Obtain minimum hot fluid temperature, and maximum cold fluid temperature
if strcmp(fluidH.read,'CP')
    THmin = max([CP1(0,0,0,'Tmin',fluidH.handle),TC1]);
else
    THmin = TC1;
end
if strcmp(fluidC.read,'CP')
    TCmax = min([CP1(0,0,0,'Tmax',fluidC.handle),TH2]);
else
    TCmax = TH2;
end

% Obtain preliminary minimum and maximum enthalpy outlets (hot outlet
% cannot be colder than cold inlet, and vice-versa)
hH1_min = RP1('PT_INPUTS',pH2,THmin,'H',fluidH);
hC2_max = RP1('PT_INPUTS',pC1,TCmax,'H',fluidC);

% Compute average 'overall' specific heat capacities
CpHmean = (hH2 - hH1_min)/(TH2-THmin);
CpCmean = (hC2_max - hC1)/(TCmax-TC1);

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
        if any([par<=TC1,par>=TCmax]), error('par must be TC1<par<TCmax'); end
        
    case 4
        % Set TH1 = par, and compute mH. Mass flow rate of cold fluid must
        % be previously specified
        if mC == 0, error('mC must be known in mode==4'); end
        if any([par<=THmin,par>=TH2])
            warning(strcat('Condition THmin<par<TH2 must be true in mode==4. ',...
                'Changing to mode==2 with Crat=1'));
            mode=2;
            Crat = 1.0; %Crat = mH*CpH / (mC*CpC)
            mH = Crat*mC*CpCmean/CpHmean;
            stateH.mdot = mH;
        end
        
    case 5
        % Set TC2 = par, and compute mH. Mass flow rate of cold fluid must
        % be previously specified
        if mC == 0, error('mC must be known in mode==5'); end
        if any([par<=TC1,par>=TCmax]), error('par must be TC1<par<TCmax'); end
        
    case 6
        % Set TH1 = par, and compute mC. Mass flow rate of hot fluid must
        % be previously specified
        if mH == 0, error('mH must be known in mode==6'); end
        if any([par<=THmin,par>=TH2]), error('par must be THmin<par<TH2'); end
        
    otherwise
        error('Invalid operation mode')
end

% Run algorithm according to the different models and operation modes
switch model
    case {'eff','UA'}
        
        % Set options for matlab root-finders (if needed)
        options = []; %optimset('Display','iter');
        
        if strcmp(model,'eff')
            compare = 'DTmin';
            ref = 0;
        elseif strcmp(model,'UA')
            compare = 'UA';
            ref = UA_ref;
        end
        
        switch mode
            case {0,1,2}
                % Compute preliminary QMAX (hot outlet cannot be colder than cold inlet,
                % and vice-versa) and update hH1_min accordingly
                QMAX0 = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)]);
                hH1_min = hH2 - QMAX0/mH;
                
                % Find value of hH1 for which DTmin=ref or UA=ref
                f1  = @(hH1) compute_TQ(fluidH,fluidC,mH,mC,hH2,hC1,pH2,pC1,NX,'hH1',hH1,compare,ref);
                %plot_function(f1,hH1_min,hH2,1000,11);
                %symlog(gca,'y')
                hH1 = fzero(f1,[hH1_min,hH2],options);
                
                % Compute total heat transfer
                if strcmp(model,'eff')
                    QMAX = mH*(hH2 - hH1);
                    QT   = QMAX*eff;
                elseif strcmp(model,'UA')
                    QT   = mH*(hH2 - hH1);
                end
                
                % Determine outlet enthalpies and pressures
                hC2 = hC1 + QT/mC;
                hH1 = hH2 - QT/mH;
                pH1 = pH2*(1-ploss);
                pC2 = pC1*(1-ploss);
                
            case 3
                % Set outlet conditions of cold fluid
                TC2 = par;
                hC2 = RP1('PT_INPUTS',pC1,TC2,'H',fluidC);
                pC2 = pC1*(1-ploss);
                
                % Compute preliminary QMAX (hot outlet cannot be colder than cold
                % inlet) and set boundaries accordingly
                QMAX0 = mH*(hH2 - hH1_min);
                mCmin = mH*eps;
                mCmax = QMAX0/(hC2 - hC1)*(1+1e-3); %necessary to find root
                
                % Find value of mC for which DTmin=ref or UA=ref
                f1 = @(mC) compute_TQ(fluidH,fluidC,mH,mC,hH2,hC1,pH2,pC1,NX,'hC2',hC2,compare,ref);
                %plot_function(f1,mCmin,mCmax,100,11);
                mC = fzero(f1,[mCmin,mCmax],options);
                
                % Compute total heat transfer (and update mC if necessary)
                if strcmp(model,'eff')
                    QMAX = mC*(hC2 - hC1);
                    QT   = QMAX*eff;
                    mC = QT/(hC2 - hC1);
                elseif strcmp(model,'UA')
                    QT   = mC*(hC2 - hC1);
                end
                stateC.mdot = mC;
                
                % Update outlet conditions of hot fluid
                hH1 = hH2 - QT/mH;
                pH1 = pH2*(1-ploss);
                
            case 4
                % Set outlet conditions of hot fluid
                TH1 = par;
                hH1 = RP1('PT_INPUTS',pH2,TH1,'H',fluidH);
                pH1 = pH2*(1-ploss);
                
                % Compute preliminary QMAX (cold outlet cannot be hotter than hot
                % inlet) and set boundaries accordingly
                QMAX0 = mC*(hC2_max - hC1);
                mHmin = mC*eps;
                mHmax = QMAX0/(hH2 - hH1)*(1+1e-3); %necessary to find root
                
                % Find value of mH for which DTmin=ref or UA=ref
                f1 = @(mH) compute_TQ(fluidH,fluidC,mH,mC,hH2,hC1,pH2,pC1,NX,'hH1',hH1,compare,ref);
                mH = fzero(f1,[mHmin,mHmax],options);
                
                % Compute total heat transfer (and update mH if necessary)
                if strcmp(model,'eff')
                    QMAX = mH*(hH2 - hH1);
                    QT   = QMAX*eff;
                    mH   = QT/(hH2 - hH1);                    
                elseif strcmp(model,'UA')
                    QT   = mH*(hH2 - hH1);
                end
                stateH.mdot = mH;
                
                % Update outlet conditions of cold fluid
                hC2 = hC1 + QT/mC;
                pC2 = pC1*(1-ploss);
                
            case 5
                error('not implemented yet')
                
            case 6
                % Set outlet conditions of hot fluid
                TH1 = par;
                hH1 = RP1('PT_INPUTS',pH2,TH1,'H',fluidH);
                pH1 = pH2*(1-ploss);
                
%                 % Compute total heat transfer
%                 QT  = mH*(hH2-hH1);
% 
%                 if strcmp(model,'eff')
%                     % Compute QMAX according to eff and obtain outlet
%                     % conditions in ideal case (DTmin=0)
%                     QMAX = QT/eff;
%                     hH1_id = hH2 - QMAX/mH;
%                     
%                     % Compute mCmin and mCmax according to QMAX
%                     mCmin = QMAX/(hC2_max - hC1)*(0.99);
%                     mCmax = mCmin*1e3; %necessary to find root
%                     
%                     % Using hH1_id, find value of mC for which DTmin=0
%                     f1 = @(mC) compute_TQ(fluidH,fluidC,mH,mC,hH2,hC1,pH2,pC1,NX,'hH1',hH1_id,'DTmin',0);
%                     %plot_function(f1,mCmin,mCmax,100,11);
%                     %symlog(gca,'y')
%                     %keyboard
%                     mC = fzero(f1,[mCmin,mCmax],options);
%                     
%                 elseif strcmp(model,'UA')
%                     % Compute mCmin and mCmax according to QT
%                     mCmin = QT/(hC2_max - hC1)*(0.99);
%                     mCmax = mCmin*1e3; %necessary to find root
%                     
%                     % Find value of mC for which UA=ref
%                     f1 = @(mC) compute_TQ(fluidH,fluidC,mH,mC,hH2,hC1,pH2,pC1,NX,'hH1',hH1,'UA',ref);
%                     %plot_function(f1,mCmin,mCmax,100,11);
%                     %symlog(gca,'y')
%                     %keyboard
%                     mC = fzero(f1,[mCmin,mCmax],options);
%                 end
%                 stateC.mdot = mC;
                
                % Compute total heat transfer and compute mCmin and mCmax
                % accordingly
                QT  = mH*(hH2-hH1);
                mCmin = QT/(hC2_max - hC1)*(0.99);
                mCmax = mCmin*1e3; %necessary to find root
                
                % Find value of mC for which DTmin=ref or UA=ref
                f1 = @(mC) compute_TQ(fluidH,fluidC,mH,mC,hH2,hC1,pH2,pC1,NX,'hH1',hH1,compare,ref);
                %plot_function(f1,mCmin,mCmax,100,11);
                %symlog(gca,'y')
                %keyboard
                mC = fzero(f1,[mCmin,mCmax],options);                
                
                % Update value of mC according to effectiveness value (this
                % will also affect the cold enthalpy outlet)
                if strcmp(model,'eff')
                    mC = mC/eff;
                end
                stateC.mdot = mC;
                                
                % Update outlet conditions of cold fluid
                hC2 = hC1 + QT/mC;
                pC2 = pC1*(1-ploss);
        end
        
        % Save variables into C and H stream objects
        [TC,TH,hC,hH,QS] = compute_TQ(fluidH,fluidC,mH,mC,hH2,hC1,pH2,pC1,NX,'hH1',hH1);
        C.mdot = mC;
        H.mdot = mH;
        C.T = TC;
        H.T = TH;
        C.h = hC;
        H.h = hH;
        HX.C = C;
        HX.H = H;
        HX.QS  = QS;
        HX.AS  = [];
        
        
    case 'geom'
        
        if any(mode==[3,4,5,6])
            error('Operation modes 3, 4, 5 and 6 not implemented yet for geometry-based heat exchanger model');
        end
        
        % Import mH and mC into stream objects
        H.mdot = mH;
        C.mdot = mC;
        
        % Set initial conditions for iteration procedure
        % Pressures
        H.p = ones(NX+1,1)*H.pin;
        C.p = ones(NX+1,1)*C.pin;
        
        % Compute derived geometric parameters and mass fluxes
        [C, H, HX] = shell_and_tube_geom(C, H, HX);
        
        % Find value of hH1 for which computed area equals specified area
        f1 = @(hH1) compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,HX);
        %plot_function(f1,hH1_min,hH2,100,31);
        opt = optimset('TolX',(hH2-hH1_min)/1e12,'Display','notify');
        hH1 = fzero(f1,[hH1_min,hH2],opt);
        
        % Obtain output parameters for converged solution
        [C,H,QS,AS] = compute_area(hH1,fluidH,fluidC,H,C,mH,mC,hH2,hC1,HX);
        
        % Save outlet conditions
        hH1 = H.h(1);
        pH1 = H.p(1);
        hC2 = C.h(NX+1);
        pC2 = C.p(NX+1);
        
        % Save variables into HX structure
        HX.C  = C;
        HX.H  = H;
        HX.QS = QS;
        HX.AS = AS;
        
end

% Update states
stateH.h = hH1;
stateH.p = pH1;
stateH   = update_state(stateH,fluidH.handle,fluidH.read,fluidH.TAB,2);
stateC.h = hC2;
stateC.p = pC2;
stateC   = update_state(stateC,fluidC.handle,fluidC.read,fluidC.TAB,2);

% Update average specific heat capacities and compute Cmin and NTU. Save
% into HX structure. Also save DppC nd DppH.
TH1 = stateH.T;
TC2 = stateC.T;
CpHmean = (hH2 - hH1)/(TH2-TH1);
CpCmean = (hC2 - hC1)/(TC2-TC1);
Cmin  = min([mC*CpCmean,mH*CpHmean]);
dQ    = QS(2:NX+1)-QS(1:NX);
DT_AV = 0.5*(HX.H.T(1:NX)+HX.H.T(2:NX+1)) - 0.5*(HX.C.T(1:NX)+HX.C.T(2:NX+1));
UA    = sum(dQ./DT_AV);
NTU   = UA/Cmin;
DppH  = (pH2-pH1)/pH2;
DppC  = (pC1-pC2)/pC1;
HX.H.Cp_mean = CpHmean;
HX.C.Cp_mean = CpCmean;
HX.Cmin = Cmin;
HX.NTU  = NTU;
HX.DppH = DppH;
HX.DppC = DppC;

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
fluidH.state(iL,iH+1) = stateH; % Result goes into next state
fluidH.stage(iL,iH)   = stageH; % Result stays in current stage
fluidC.state(iL,iC+1) = stateC; % Result goes into next state
fluidC.stage(iL,iC)   = stageC; % Result stays in current stage

% Update mass flow rates for inlet state, if necessary
if any(mode==[1,2,3,4,5,6])
    fluidH.state(iL,iH).mdot = mH;
    fluidC.state(iL,iC).mdot = mC;
end

% Reverse swap, if necessary
if swap == 1
    fluidH0 = fluidH;
    fluidH  = fluidC;
    fluidC  = fluidH0;
    iH0 = iH;
    iH  = iC;
    iC  = iH0;
end

% Increase stage counter
iH = iH + 1;
iC = iC + 1;

end


%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = compute_TQ(fluidH,fluidC,mH,mC,hH2,hC1,pH2,pC1,n,mode,hout,varargin)
% Compute the TQ diagram of a two-stream counter-flow heat exchanger.
%   The "mode" string controls which enthalpy outlet "hout" is specified,
%   either mode='hH1' or mode='hC2'.
%
%   The optional inputs, "varargin", contain "compare" and "ref".
%   The "compare" string controls which parameter is to be compared to the
%   "ref" parameter, either compare='DTmin' (minimum temperature
%   difference) or compare='UA' (overall heat transfer coefficient).
%
%   If the number of optional outputs is 1, varargout contains "solution"
%   (the difference between the computed value and "ref").
%
%   If the number of optional outputs is 5, they contain the following
%   arrays, in order: TC, TH, hC, hH, QS.
%
%   Usage examples:
%   [TC,TH,hC,hH,QS] = compute_TQ(mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,NX,'hH1',hH1);
%   solution         = compute_TQ(mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,NX,'hH1',hH1,'UA',UA_ref);


% Check and assign variable inputs
if isempty(varargin)
    
elseif length(varargin)==2
    compare = varargin{1};
    ref = varargin{2};
    
else
    error('incorrect number of inputs')
end

% Compute TQ diagram based on either hH1 or hC2
switch mode
    case 'hH1'
        
        % Set outlet enthalpy
        hH1 = hout;
        
        % Compute temperature distribution of hot stream
        hH  = linspace(hH1,hH2,n+1)';
        pH  = pH2*ones(size(hH));
        TH  = RP1('HmassP_INPUTS',hH,pH,'T',fluidH);
        
        % Compute temperature distribution of cold stream
        QS  = (hH - hH1)*mH; % cummulative heat transfer
        hC  = hC1 + QS/mC;
        pC  = pC1*ones(size(hC));
        TC  = RP1('HmassP_INPUTS',hC,pC,'T',fluidC);
        
    case 'hC2'
        
        % Set outlet enthalpy
        hC2 = hout;
        
        % Compute temperature distribution of cold stream
        hC  = linspace(hC1,hC2,n+1)';
        pC  = pC1*ones(size(hC));
        TC  = RP1('HmassP_INPUTS',hC,pC,'T',fluidC);
        
        % Compute temperature distribution of hot stream
        QS  = (hC - hC1)*mC; % cummulative heat transfer
        hH1 = hH2 - QS(n+1)/mH;
        hH  = hH1 + QS/mH;
        pH  = pH2*ones(size(hH));
        TH  = RP1('HmassP_INPUTS',hH,pH,'T',fluidH);
        
    otherwise
        error('not implemented')
end

% % To visualise the temperature distribution every time the function is
% % called, uncomment the lines below
% figure(10)
% plot(QS./QS(end),TH,'r'); hold on;
% plot(QS./QS(end),TC,'b'); hold off;
% xlabel('Cumulative heat transfer')
% ylabel('Temperature')
% keyboard

% Compare computed 'DTmin' or 'UA' to reference values
if length(varargin) == 2
    
    % Compute temperature difference between the two streams
    DT = TH - TC;
    
    if strcmp(compare,'DTmin')
        
        solution = min(DT) - ref;
        
    elseif strcmp(compare,'UA')
        
        % Compute average temperature difference
        DT_AV = 0.5*(DT(2:n+1)+DT(1:n));
        
        % Compute heat transfer at each section
        dQ = QS(2:n+1)-QS(1:n);
        
        % Compute required heat transfer coefficient
        dUA = dQ./DT_AV;
        UA = sum(dUA);
        
        % Find difference between reference and computed
        solution = ref - UA;
        
        % Control physically impossible solutions
        if any(DT_AV <= 0)
            solution = - ref;
        end
        
    else
        error('not implemented')
    end
    
else
    
    solution = NaN;
    
end

% Assign output arguments
if nargout == 1
    varargout{1} = solution;
elseif nargout == 5
    varargout{1} = TC;
    varargout{2} = TH;
    varargout{3} = hC;
    varargout{4} = hH;
    varargout{5} = QS;
else
    error('Incorrect number of outputs')
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
    if any(cond1)
        warning('DpC exceeds 20%!');
    end
    if any(cond2)
        warning('DpH exceeds 20%!');
    end
    
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
% Control physically impossible solutions
if any(DT_AV <= 0)
    solution = - C.A;
end

if nargout == 1
    varargout{1} = solution;
else
    % Compute cumulative heat transfer area (cold side)
    AS = zeros(size(QS));
    for i=1:(length(AS)-1)
        AS(i+1) = AS(i) + dAC(i);
    end
    varargout{1} = C;
    varargout{2} = H;
    varargout{3} = QS;
    varargout{4} = AS;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%