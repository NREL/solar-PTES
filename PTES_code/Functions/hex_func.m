function [HX, fluidH, iH, fluidC, iC] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par)
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
%   In mode == 5, TH1 is specified (TH1=par) and mC (unknown) is computed
%   
%   The HX object/structure can operate according to three different models
%   If model = 'eff', the heat exchanger effectiveness and pressure loss
%   are specified
%   If model = 'UA', the overall heat transfer coefficient and pressure
%   loss are specified
%   If model = 'geom, the heat exchanger geometry is specified


% Extract parameters from HX structure according to selected model
model = HX.model;
stage_type = HX.stage_type;
NX  = HX.NX;
switch model
    case 'DT'
        DT     = HX.DT;
        ploss  = HX.ploss;
        
    case 'eff'
        eff    = HX.eff;
        ploss  = HX.ploss;
        
    case 'UA'
        UA_ref = HX.UA;
        ploss  = HX.ploss;
        
    case 'geom'
        % Set heat exchanger geometry (first time only)
        if ~HX.Lgeom_set
            [HX] = set_hex_geom2(HX, iL, fluidH, iH, fluidC, iC, mode, par);
        end
        
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
    error(strcat('swap not implemented for hx_class and set_hex_geom.',...
        'make sure that fluidH is fluidH and fluidC is fluidC when calling hex_func'))
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
H = stream; H.mdot = mH; H.name = fluidH.name;
C = stream; C.mdot = mC; C.name = fluidC.name;
H.read = fluidH.read; H.handle = fluidH.handle; H.HEOS = fluidH.HEOS;
C.read = fluidC.read; C.handle = fluidC.handle; C.HEOS = fluidC.HEOS;
H.pin = pH2;
C.pin = pC1;

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
        if any([par<=TC1,par>=TCmax])
            warning(strcat('Condition TC1<par<TCmax must be true in mode==3. ',...
                'Changing to mode==1 with Crat=1'));
            mode = 1;
            Crat = 1.0; %Crat = mH*CpH / (mC*CpC)
            mC = mH*CpHmean/(CpCmean*Crat);
            stateC.mdot = mC;
        end
        
    case 4
        % Set TH1 = par, and compute mH. Mass flow rate of cold fluid must
        % be previously specified
        if mC == 0, error('mC must be known in mode==4'); end
        if mH == 0, mH = []; end
        if any([par<=THmin,par>=TH2])
            warning(strcat('Condition THmin<par<TH2 must be true in mode==4. ',...
                'Changing to mode==2 with Crat=1'));
            mode = 2;
            Crat = 1.0; %Crat = mH*CpH / (mC*CpC)
            mH = Crat*mC*CpCmean/CpHmean;
            stateH.mdot = mH;
        end
        
    case 5
        % Set TH1 = par, and compute mC. Mass flow rate of hot fluid must
        % be previously specified
        if mH == 0, error('mH must be known in mode==5'); end
        if any([par<=THmin,par>=TH2]), error('par must be THmin<par<TH2'); end
        
    otherwise
        error('Invalid operation mode')
end

% Run algorithm according to the different models and operation modes
switch model
    case {'eff','UA','DT'}
        
        % Set options for matlab root-finders (if needed)
        options = []; %optimset('Display','iter');
        
        switch model
            case 'eff'
                compare = 'DTmin';
                ref = 0;
            case 'UA'
                compare = 'UA';
                ref = UA_ref;
            case 'DT'
                compare = 'DTmin';
                ref = DT;
        end
        
        % Set outlet pressures
        if ~isempty(HX.plossH0)
            plossH = HX.plossH0;
        else
            plossH = ploss;
        end
        if ~isempty(HX.plossC0)
            plossC = HX.plossC0;
        else
            plossC = ploss;
        end
        pH1 = pH2*(1-plossH);
        pC2 = pC1*(1-plossC);
        
        switch mode
            case {0,1,2}
                % Compute preliminary QMAX (hot outlet cannot be colder than cold inlet,
                % and vice-versa) and update hH1_min accordingly
                QMAX0 = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)])*(1.01); %necessary to find root
                hH1_min = hH2 - QMAX0/mH;
                
                % Find value of hH1 for which DTmin=ref or UA=ref
                f1  = @(hH1) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref);
                %{
                plot_function(f1,hH1_min,hH2,100,11);
                symlog(gca,'y')                
                f2  = @(hH1) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref,true);
                plot_function(f2,hH1_min,hH2,100,11);
                %}
                hH1 = fzero(f1,[hH1_min,hH2],options);
                
                % Compute total heat transfer
                switch model
                    case 'eff'
                        QMAX = mH*(hH2 - hH1);
                        QT   = QMAX*eff;
                    case {'DT','UA'}
                        QT   = mH*(hH2 - hH1);
                end
                
                % Determine outlet enthalpies
                hC2 = hC1 + QT/mC;
                hH1 = hH2 - QT/mH;                
                
            case 3
                % Set outlet conditions of cold fluid
                TC2 = par;
                hC2 = RP1('PT_INPUTS',pC2,TC2,'H',fluidC);   
                
                % Compute preliminary QMAX (hot outlet cannot be colder than cold
                % inlet) and set boundaries accordingly
                QMAX0 = mH*(hH2 - hH1_min);
                mCmin = mH*eps;
                mCmax = QMAX0/(hC2 - hC1)*(1.01); %necessary to find root
                
                % Find value of mC for which DTmin=ref or UA=ref
                f1 = @(mC) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hC2',hC2,compare,ref);
                %{
                plot_function(f1,mCmin,mCmax,100,11);
                f2 = @(mC) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hC2',hC2,compare,ref,true);
                plot_function(f2,mCmin,mCmax,100,11);
                %}
                % Check whether f1 changes sign over interval
                f1min = f1(mCmin) ;
                f1max = f1(mCmax) ;
                % If they don't change sign, choose mC that has minimum boundary
                if f1min*f1max >= 0
                    warning('Cold stream may not be heated to desired temperature');
                    if min(f1min,f1max) == f1min
                        mC = mCmin ;
                    else
                        mC = mCmax ;
                    end
                else
                    mC = fzero(f1,[mCmin,mCmax],options);
                end
                
                % Compute total heat transfer (and update mC if necessary)
                switch model
                    case 'eff'
                        QMAX = mC*(hC2 - hC1);
                        QT   = QMAX*eff;
                        mC = QT/(hC2 - hC1);
                    case {'DT','UA'}
                        QT   = mC*(hC2 - hC1);
                end
                stateC.mdot = mC;
                
                % Update outlet conditions of hot fluid
                hH1 = hH2 - QT/mH;
                
            case 4
                % Set outlet conditions of hot fluid
                TH1 = par;
                hH1 = RP1('PT_INPUTS',pH1,TH1,'H',fluidH);
                
                % Compute preliminary QMAX (cold outlet cannot be hotter than hot
                % inlet) and set boundaries accordingly
                QMAX0 = mC*(hC2_max - hC1);
                mHmin = eps(mC);
                mHmax = QMAX0/(hH2 - hH1)*(1.01); %necessary to find root
                
                % Find value of mH for which DTmin=ref or UA=ref
                f1 = @(mH) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref);
                %{
                plot_function(f1,mHmin,mHmax,100,11)
                f2 = @(mH) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref,true);
                plot_function(f2,mHmin,mHmax,100,11)
                %}
                mH = fzero(f1,[mHmin,mHmax],options);
                
                % Compute total heat transfer (and update mH if necessary)
                switch model
                    case 'eff'
                        QMAX = mH*(hH2 - hH1);
                        QT   = QMAX*eff;
                        mH   = QT/(hH2 - hH1);
                    case {'UA','DT'}
                        QT   = mH*(hH2 - hH1);
                end
                stateH.mdot = mH;
                
                % Update outlet enthalpy of cold fluid
                hC2 = hC1 + QT/mC;
                
            case 5
                % Set outlet conditions of hot fluid
                TH1 = par;
                hH1 = RP1('PT_INPUTS',pH1,TH1,'H',fluidH);                
                
                % Compute total heat transfer and compute mCmin and mCmax
                % accordingly
                QT  = mH*(hH2-hH1);
                mCmin = QT/(hC2_max - hC1)*(0.99);
                mCmax = mCmin*1e3; %necessary to find root
                
                % Find value of mC for which DTmin=ref or UA=ref
                f1 = @(mC) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref);
                %{
                plot_function(f1,mCmin,mCmax,100,11);
                symlog(gca,'y')
                keyboard
                %}
                mC = fzero(f1,[mCmin,mCmax],options);
                
                % Update value of mC according to effectiveness value (this
                % will also affect the cold enthalpy outlet)
                if strcmp(model,'eff')
                    mC = mC/eff;
                end
                stateC.mdot = mC;
                                
                % Update outlet conditions of cold fluid
                hC2 = hC1 + QT/mC;
        end
        
        % Save variables into C and H stream objects
        [TC,TH,hC,hH,QS] = compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1);
        C.mdot = mC;
        H.mdot = mH;
        C.T = TC;
        H.T = TH;
        C.h = hC;
        H.h = hH;
        HX.C(iL) = C;
        HX.H(iL) = H;
        HX.QS(iL,:) = QS;
        HX.AS  = [];
        HX.C(iL).pin = pC1;
        HX.H(iL).pin = pH2;
        
        
    case 'geom'
        
        % Import mH and mC into stream objects (one of these might still be
        % set to 0 at this stage -unknown-, depending on the mode)
        H.mdot = mH;
        C.mdot = mC;
        
        switch mode
            case {0,1,2}
                QMAX0 = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)])*(1.01); %necessary to find root
                hH1_min = hH2 - QMAX0/mH;
                
                % Find value of hH1 for which computed area equals specified area
                f1 = @(hH1) compute_area(HX,H,C,mH,mC,hH2,hC1,'hH1',hH1);
                %plot_function(f1,hH1_min,hH2,100,31);
                %keyboard
                opt = optimset('TolX',(hH2-hH1_min)/1e12,'Display','notify');
                hH1 = fzero(f1,[hH1_min,hH2],opt);
                
            case 4
                % Set outlet conditions of hot fluid. Assume small effect
                % of pressure loss when computing outlet enthalpy from the
                % objective outlet temperature.
                TH1 = par;
                hH1 = RP1('PT_INPUTS',pH2,TH1,'H',fluidH);
                
                % Compute preliminary QMAX (cold outlet cannot be hotter
                % than hot inlet) and set boundaries accordingly
                QMAX0 = mC*(hC2_max - hC1);
                mHmin = eps(mC);
                mHmax = QMAX0/(hH2 - hH1)*(1.01); %necessary to find root
                
                % Find value of mH for which computed area equals specified area
                f1  = @(mH) compute_area(HX,H,C,mH,mC,hH2,hC1,'hH1',hH1);
                opt = optimset('TolX',(mHmax-mHmin)/1e12,'Display','notify');
                mH  = fzero(f1,[mHmin,mHmax],opt);
                
            case 5
                % Set outlet conditions of hot fluid. Assume small effect
                % of pressure loss when computing outlet enthalpy from the
                % objective outlet temperature.
                TH1 = par;
                hH1 = RP1('PT_INPUTS',pH2,TH1,'H',fluidH);
                
                % Compute total heat transfer and compute mCmin and mCmax
                % accordingly
                QT    = mH*(hH2-hH1);
                mCmin = QT/(hC2_max - hC1)*1.01;
                mCmax = mCmin*100; %necessary to find root
                
                % Find value of mC for which computed area equals specified area
                %keyboard
                f1  = @(mC) compute_area(HX,H,C,mH,mC,hH2,hC1,'hH1',hH1);
                %plot_function(f1,mCmin,mCmax,5,15,'semilogx');
                opt = optimset('TolX',(mCmax-mCmin)/1e12,'Display','notify');
                mC = fzero(f1,[mCmin,mCmax],opt);
                
            case 3
                error(['Operation mode 3 has not been implemented yet',...
                    'for the geometry-based heat exchanger model']);
        end
        
        % Obtain output parameters for converged solution
        [C,H,QS,AS] = compute_area(HX,H,C,mH,mC,hH2,hC1,'hH1',hH1);
        
        % Save outlet conditions
        hH1 = H.h(1);
        pH1 = H.p(1);
        hC2 = C.h(NX+1);
        pC2 = C.p(NX+1);
        
        % Save variables into HX structure
        HX.C(iL)  = C;
        HX.H(iL)  = H;
        HX.QS(iL,:) = QS;
        HX.AS = AS;
end

% Update states
stateH.h = hH1;
stateH.p = pH1;
stateH   = update_state(stateH,fluidH,2);
stateC.h = hC2;
stateC.p = pC2;
stateC   = update_state(stateC,fluidC,2);

% Update average specific heat capacities and compute Cmin and NTU. Save
% into HX structure. Also save DppC nd DppH.
TH1 = stateH.T;
TC2 = stateC.T;
CpHmean = (hH2 - hH1)/(TH2-TH1);
CpCmean = (hC2 - hC1)/(TC2-TC1);
Cmin  = min([mC*CpCmean,mH*CpHmean]);
dQ    = QS(2:NX+1)-QS(1:NX);
DT_AV = 0.5*(HX.H(iL).T(1:NX)+HX.H(iL).T(2:NX+1)) - 0.5*(HX.C(iL).T(1:NX)+HX.C(iL).T(2:NX+1));
UA    = sum(dQ./DT_AV);
NTU   = UA/Cmin;
DppH  = (pH2-pH1)/pH2;
DppC  = (pC1-pC2)/pC1;

dTa = HX.H(iL).T(1) - HX.C(iL).T(1) ;
dTb = HX.H(iL).T(end) - HX.C(iL).T(end) ;

HX.H(iL).Cp_mean = CpHmean;
HX.C(iL).Cp_mean = CpCmean;
HX.Cmin(iL) = Cmin;
HX.NTU(iL)  = NTU;
HX.DppH(iL) = DppH;
HX.DppC(iL) = DppC;
HX.UA(iL)   = UA ;
HX.LMTD(iL) = (dTa - dTb) / log(dTa / dTb) ;

% If this is the first time that hex_func is called, save the initial
% values UA0, NTU0 and LMTD0
if isempty(HX.UA0)
    HX.UA0 = HX.UA(iL) ;
    HX.NTU0 = HX.NTU(iL) ;
    HX.LMTD0 = HX.LMTD(iL) ;
end

% *** DELETE EVENTUALLY >>>
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

% <<< DELETE EVENTUALLY ***

% Compute losses, but insert into hx_class
% Entropy change
DsH         = stateH.s - sH2;
DsC         = stateC.s - sC1;
% Hot stream
HX.Dh(iL,1)   = stateH.h - hH2;
HX.sirr(iL,1) = (stateH.mdot*DsH + stateC.mdot*DsC)/stateH.mdot;
HX.q(iL,1)    = stageH.Dh;
HX.w(iL,1)    = 0;

% Cold stream
HX.Dh(iL,2)   = stateC.h - hC1;
HX.sirr(iL,2) = (stateC.mdot*DsC + stateH.mdot*DsH)/stateC.mdot;
HX.q(iL,2)    = stageC.Dh;
HX.w(iL,2)    = 0;

if strcmp(HX.stage_type,'regen')
    HX.sirr(iL,1) = 0; %to avoid counting the lost work twice
end

% Update mass flow rates for inlet state, if necessary
if any(mode==[1,2,3,4,5])
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

function varargout = compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,n,mode,hout,varargin)
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

visualise = 0;
% Check and assign variable inputs
if isempty(varargin)
    
elseif length(varargin)==2
    compare = varargin{1};
    ref = varargin{2};
    
elseif length(varargin)==3
    compare = varargin{1};
    ref = varargin{2};
    visualise = varargin{3};
    
else
    error('incorrect number of inputs')
end

% Set pressure arrays
pH  = linspace(pH1,pH2,n+1)';
pC  = linspace(pC1,pC2,n+1)';

% Compute TQ diagram based on either hH1 or hC2
switch mode
    case 'hH1'
        
        % Set outlet enthalpy
        hH1 = hout;
        
        % Compute temperature distribution of hot stream
        hH  = linspace(hH1,hH2,n+1)';
        TH  = RP1('HmassP_INPUTS',hH,pH,'T',fluidH);
        
        % Compute temperature distribution of cold stream
        QS  = (hH - hH1)*mH; % cummulative heat transfer
        hC  = hC1 + QS/mC;
        TC  = RP1('HmassP_INPUTS',hC,pC,'T',fluidC);
        
    case 'hC2'
        
        % Set outlet enthalpy
        hC2 = hout;
        
        % Compute temperature distribution of cold stream
        hC  = linspace(hC1,hC2,n+1)';
        TC  = RP1('HmassP_INPUTS',hC,pC,'T',fluidC);
        
        % Compute temperature distribution of hot stream
        QS  = (hC - hC1)*mC; % cummulative heat transfer
        hH1 = hH2 - QS(n+1)/mH;
        hH  = hH1 + QS/mH;
        TH  = RP1('HmassP_INPUTS',hH,pH,'T',fluidH);
        
    otherwise
        error('not implemented')
end

% Compare computed 'DTmin' or 'UA' to reference values
if any(length(varargin) == [2,3])
    
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

if visualise
    % Visualise the temperature distribution every time the function is
    % called
    figure(10)
    plot(QS./QS(end),TH,'r'); hold on;
    plot(QS./QS(end),TC,'b'); hold off;
    xlabel('Cumulative heat transfer')
    ylabel('Temperature')
    keyboard
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

function varargout = compute_area(HX,H,C,mH,mC,hH2,hC1,mode,hout,varargin)
%COMPUTE_AREA Solve the TQ and TA diagrams diagrams of a two-stream
%counter-flow heat exchanger.
%
%   For a given fluid outlet enthalpy (hH1 or hC2), compute the TQ diagram,
%   the properties of the fluids at each point and the corresponding heat
%   transfer area of the heat exchanger. Compare that to the reference heat
%   transfer area and return the difference between the two.
%
%   The "mode" string controls which enthalpy outlet "hout" is specified,
%   either mode='hH1' or mode='hC2'.
%
%   The optional input, "varargin", contains the logical argument
%   "visualise". If "visualise" is set to true, the compute_area function
%   pauses at the end of the internal iteration procedure and plots the
%   result.

% Select mode and assing value of known enthalpy outlet
switch mode
    case 'hH1'
        hH1 = hout;
    case 'hC2'
        error('not implemented yet')
    otherwise
        error('not implemented')
end

% Check and assign variable inputs
if isempty(varargin)
    visualise = 0;
elseif length(varargin)==1
    visualise = varargin{1};
else
    error('not implemented')
end

% Extract parameters
NX = HX.NX;

% Compute mass fluxes
[C, H, HX] = shell_and_tube_geom(C, H, HX);

% Set initial conditions for iteration procedure
% Pressures
H.p = ones(NX+1,1)*H.pin;
C.p = ones(NX+1,1)*C.pin;

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
    C = stream_update(C,2);
    % Hot stream
    H = stream_update(H,2);
    
    % COMPUTE AVERAGED TEMPERATURE ARRAYS
    H.T_AV = 0.5*(H.T(1:NX) + H.T(2:NX+1));
    C.T_AV = 0.5*(C.T(1:NX) + C.T(2:NX+1));
    DT_AV  = H.T_AV - C.T_AV;
    
    % Break loop if H.T < C.T at any point
    if any(H.T <= C.T)
        impossible = true;
        AC = Inf;
        %warning(['impossible condition reached.',...
        %    'breaking loop and proceeding'])
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
    
    % In a situation where H.h(NX+1)==H.h(1) (i.e. no heat exchange), the
    % code above results in dAC=0, AC=0 and Dp_H=NaN, Dp_C=NaN. If so, set
    % Dp_H and Dp_C to zero and proceed.
    if any(isnan([Dp_H;Dp_C]))
        if abs(H.h(NX+1)-H.h(1)) < 10*eps(H.h(1))
            Dp_H = zeros(size(Dp_H));
            Dp_C = zeros(size(Dp_C));
        else
            error('NaN value found in hex_func!')
        end
    end
    
    % Update pressure profiles. Assume a linear profile (to avoid computing
    % a slow 'for loop') and limit max pressure loss to 80%.
    DppH = abs(sum(Dp_H))/H.pin;
    DppC = abs(sum(Dp_C))/C.pin;
    if DppH > 0.8
        DppH = 0.8;
        warning('DpH exceeds 80%!');
    end
    if DppC > 0.8
        DppC = 0.8;
        warning('DpC exceeds 80%!');
    end
    H.p  = linspace(H.pin*(1-DppH),H.pin,NX+1)';
    C.p  = linspace(C.pin,C.pin*(1-DppC),NX+1)';
    
    % Update convergence array
    CON = [AC; H.p; C.p]; % initial value
    
    % Compute residual
    RES(iI) = max(abs((CON - CON_0)./CON));
    
    if (RES(iI)>TOL)
        if visualise
            fprintf(1,'\n iteration = %d, RES = %.6f',iI,RES(iI));
        end
        CON_0 = CON;
    else
        if visualise
            % Make plots
            figure(10)
            plot(QS./QS(end),H.T,'r'); hold on;
            plot(QS./QS(end),C.T,'b'); hold off;
            xlabel('Cumulative heat transfer')
            ylabel('Temperature')
            legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],[C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')
            
            figure(11)
            plot(QS./QS(end),H.p/H.pin,'r-'); hold on
            plot(QS./QS(end),C.p/C.pin,'b-'); hold off
            ylim([0.90 1])
            xlabel('Cummulative heat transfer')
            ylabel('Relative pressure, p/p0')
            legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],[C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')
            
            fprintf(1,'\n\n*** Successful convergence after %d iterations***\n',iI);
            keyboard
        end
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

% If the value of 'solution' is negative, it means that the computed area
% is larger than the actual area (heat exchanger too small to achieve
% selected operating conditions). If the value of 'solution' is positive,
% the computed area is smaller than the actual area (heat exchanger too
% large for selected operating conditions).
solution = C.A - AC;
% Control physically impossible solutions
if impossible
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