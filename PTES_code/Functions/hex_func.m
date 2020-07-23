function [HX, fluidH, iH, fluidC, iC] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par)
% COMPUTE HEAT EXCHANGER OUTLET CONDITIONS
%   Description
%   TC1 and TH2 are the cold and hot temperature inlets (known)
%   TC2 and TH1 are the cold and hot temperature outlets (unknown in modes
%   0, 1 and 2)
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
%   If model = 'DT', the pinch point temperature difference and pressure
%   loss are specified
%   If model = 'geom, the heat exchanger geometry is specified


% Set inlet temperatures (nomenclature: cold inlet is position 1, hot inlet
% is position 2)
TH2 = fluidH.state(iL,iH).T;
TC1 = fluidC.state(iL,iC).T;

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
        
    case 'geom'
        % Set heat exchanger geometry (first time only)
        if ~HX.Lgeom_set
            [HX] = hex_set_geom(HX, iL, fluidH, iH, fluidC, iC, mode, par);
        end
        %fprintf('\nGH = %.2f, GC = %.2f, LH = %.3f, LC = %.3f,\n',HX.H.G,HX.C.G,HX.L2,HX.L1)
        %disp([HX.H.dpdL, HX.C.dpdL])
        %disp([HX.H.Cf, HX.C.Cf])
        %disp([HX.H.ht, HX.C.ht])
        %disp([HX.H.Re, HX.C.Re])
        %disp([HX.H.shape,', ',HX.C.shape])
        %keyboard
        
    otherwise
        error('Invalid heat exchanger model')
end

% Check which one is fluidH and which is fluidC and swap them if necessary
swap = 0;
if TC1 > TH2 % swap needed
    error(['Swap not implemented for hx_class and hex_set_geom. Make sure',...
        ' that fluidH is fluidH and fluidC is fluidC when calling hex_func'])
    %{
    swap = 1;
    fluidH0 = fluidH;
    fluidH  = fluidC;
    fluidC  = fluidH0;
    iH0 = iH;
    iH  = iC;
    iC  = iH0;
    %}
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
mH  = stateH.mdot;
mC  = stateC.mdot;

% Declare the two fluid streams
H = stream; H.mdot = mH; H.name = fluidH.name;
C = stream; C.mdot = mC; C.name = fluidC.name;
H.read = fluidH.read; H.handle = fluidH.handle; H.HEOS = fluidH.HEOS;
C.read = fluidC.read; C.handle = fluidC.handle; C.HEOS = fluidC.HEOS;
switch HX.shape
    case 'cross-flow'
        H.shape = 'circular';
        C.shape = 'cross-flow';
    otherwise
        H.shape = HX.shape;
        C.shape = HX.shape;
end
H.pin = pH2;
C.pin = pC1;
H.hin = hH2;
C.hin = hC1;

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
hH1_min = RPN('PT_INPUTS',pH2,THmin,'H',fluidH);
hC2_max = RPN('PT_INPUTS',pC1,TCmax,'H',fluidC);

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
    case {'eff','DT'}
        
        % Set options for matlab root-finders (if needed)
        options = []; %optimset('Display','iter');
        
        switch model
            case 'eff'
                compare = 'DTmin';
                ref     = (1-eff)*(TH2-TC1);
                
            case 'DT'
                compare = 'DTmin';
                ref     = DT;
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
        % Account for different pressure losses in the two channels. Set
        % the highest pressure loss equal to 'ploss' and estimate the
        % lowest pressure with the ratio of geometry-independent pressure
        % loss factors.
        %{
        if isempty(HX.plossH0) && isempty(HX.plossC0) && any(mode==[0,1,2])
            % Obtain average specific volumes for the two channels
            vH = 0.5*(1/RPN('PT_INPUTS',pH2,TH2,'D',fluidH) + 1/RPN('PT_INPUTS',pH2,TC1,'D',fluidH));
            vC = 0.5*(1/RPN('PT_INPUTS',pC1,TH2,'D',fluidC) + 1/RPN('PT_INPUTS',pC1,TC1,'D',fluidC));
            % Obtain average viscosities
            muH = 0.5*(RPN('PT_INPUTS',pH2,TH2,'VISCOSITY',fluidH) + RPN('PT_INPUTS',pH2,TC1,'VISCOSITY',fluidH));
            muC = 0.5*(RPN('PT_INPUTS',pC1,TH2,'VISCOSITY',fluidC) + RPN('PT_INPUTS',pC1,TC1,'VISCOSITY',fluidC));
            % Obtain geometry-independent pressure loss factors
            if any([muH,muC]==Inf)
                DppH_fac = mH^2.0*vH/pH2;
                DppC_fac = mC^2.0*vC/pC1;
            else
                DppH_fac = mH^1.67*muH^0.33*vH/pH2;
                DppC_fac = mC^1.67*muC^0.33*vC/pC1;
            end
            % Set pressure losses
            if DppH_fac > DppC_fac
                plossH = ploss;
                plossC = ploss*DppC_fac/DppH_fac;
            else
                plossC = ploss;
                plossH = ploss*DppH_fac/DppC_fac;
            end
            %{
            fprintf(1,'\n');
            fprintf(1,'Hot:  %10s, %5.1f bar, ploss = %6.4f\n',valid_name(fluidH.name,2),pH2/1e5,plossH)
            fprintf(1,'Cold: %10s, %5.1f bar, ploss = %6.4f\n',valid_name(fluidC.name,2),pC1/1e5,plossC)
            keyboard
            %}
        end
        %}
        pH1 = pH2*(1-plossH);
        pC2 = pC1*(1-plossC);
        
        switch mode
            case {0,1,2}
                % Compute preliminary QMAX (hot outlet cannot be colder than cold inlet,
                % and vice-versa) and update hH1_min accordingly
                QMAX0 = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)])*(1.01); %necessary to find root
                hH1_min = hH2 - QMAX0/mH;
                
                % Find value of hH1 for which DTmin=ref
                f1  = @(hH1) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref);
                %{
                plot_function(f1,hH1_min,hH2,100,11);
                symlog(gca,'y')
                f2  = @(hH1) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref,true);
                plot_function(f2,hH1_min,hH2,100,11);
                %}
                hH1 = fzero(f1,[hH1_min,hH2],options);
                
                % Compute total heat transfer
                QT  = mH*(hH2 - hH1);
                
                % Determine outlet enthalpies
                hC2 = hC1 + QT/mC;
                hH1 = hH2 - QT/mH;
                
            case 3
                % Set outlet conditions of cold fluid
                TC2 = par;
                hC2 = RPN('PT_INPUTS',pC2,TC2,'H',fluidC);
                
                % Compute preliminary QMAX (hot outlet cannot be colder than cold
                % inlet) and set boundaries accordingly
                QMAX0 = mH*(hH2 - hH1_min);
                mCmin = mH*eps;
                mCmax = QMAX0/(hC2 - hC1)*(1.01); %necessary to find root
                
                % Find value of mC for which DTmin=ref
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
                
                % Store new mC value into stateC structure
                stateC.mdot = mC;
                
                % Compute total heat transfer
                QT  = mC*(hC2 - hC1);
                
                % Update outlet conditions of hot fluid
                hH1 = hH2 - QT/mH;
                
            case 4
                % Set outlet conditions of hot fluid
                TH1 = par;
                hH1 = RPN('PT_INPUTS',pH1,TH1,'H',fluidH);
                
                % Compute preliminary QMAX (cold outlet cannot be hotter than hot
                % inlet) and set boundaries accordingly
                QMAX0 = mC*(hC2_max - hC1);
                mHmin = eps(mC);
                mHmax = QMAX0/(hH2 - hH1)*(1.01); %necessary to find root
                
                % Find value of mH for which DTmin=ref
                f1 = @(mH) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref);
                %{
                plot_function(f1,mHmin,mHmax,100,11)
                f2 = @(mH) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref,true);
                plot_function(f2,mHmin,mHmax,100,11)
                %}
                mH = fzero(f1,[mHmin,mHmax],options);
                
                % Store new mH value into stateH structure
                stateH.mdot = mH;
                
                % Compute total heat transfer
                QT  = mH*(hH2 - hH1);
                
                % Update outlet enthalpy of cold fluid
                hC2 = hC1 + QT/mC;
                
            case 5
                % Set outlet conditions of hot fluid
                TH1 = par;
                hH1 = RPN('PT_INPUTS',pH1,TH1,'H',fluidH);
                
                % Compute total heat transfer and compute mCmin and mCmax
                % accordingly
                QT    = mH*(hH2-hH1);
                mCmin = QT/(hC2_max - hC1)*(0.98);
                mCmax = mCmin*100; %necessary to find root
                
                % Find value of mC for which DTmin=ref
                f1 = @(mC) compute_TQ(fluidH,fluidC,mH,mC,hH2,pH2,pH1,hC1,pC1,pC2,NX,'hH1',hH1,compare,ref);
                %{
                plot_function(f1,mCmin,mCmax,1000,31,'semilogx');
                symlog(gca,'y')
                keyboard
                %}
                % Check whether f1 changes sign over interval. If they
                % don't change sign, choose equal heat capacity ratios
                f1min = f1(mCmin) ;
                f1max = f1(mCmax) ;
                if f1min*f1max >= 0
                    warning('Selected outlet conditions might not agree with HX performance');
                    mC = mH*CpHmean/CpCmean;
                else
                    mC = fzero(f1,[mCmin,mCmax],options);
                end
                
                % Store new mC value into stateC structure
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
        
        % Import streams
        HX.H(iL) = H;
        HX.C(iL) = C;
        
        % Import mH and mC into stream objects (one of these might still be
        % set to 0 at this stage -unknown-, depending on the mode)
        HX.H(iL).mdot = mH;
        HX.C(iL).mdot = mC;
        
        switch mode
            case {0,1,2}
                QMAX0 = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)])*(1.01); %necessary to find root
                hH1_min = hH2 - QMAX0/mH;
                hH1_max = hH2 - (hH2-hH1_min)*HX.eff*0.90;
                
                % Find value of hH1 for which computed area equals specified area
                switch HX.shape
                    case {'circular','PCHE'}
                        f1 = @(hH1) hex_compute_area(HX,iL,'hH1',hH1);
                    case 'cross-flow'
                        f1 = @(hH1) hex_compute_Xflow(HX,iL,'hH1',hH1);
                    otherwise
                        error('not implemented')
                end
                %plot_function(f1,hH1_min,hH1_max,100,31);
                %keyboard
                opt = optimset('TolX',(hH2-hH1_min)/1e12,'Display','notify');
                hH1 = fzero(f1,[hH1_min,hH1_max],opt);
                
            case 3
                % Set outlet conditions of cold fluid. Assume small effect
                % of pressure loss when computing outlet enthalpy from the
                % objective outlet temperature.
                TC2 = par;
                hC2 = RPN('PT_INPUTS',pC1,TC2,'H',fluidC);
                
                % Compute preliminary QMAX (hot outlet cannot be colder
                % than cold inlet) and set boundaries accordingly
                QMAX0 = mH*(hH2 - hH1_min);
                mCmax = QMAX0/(hC2 - hC1)*(1.01); %necessary to find root
                mCmin = mCmax/1e6;
                
                % Find value of mC for which computed area equals specified area
                switch HX.shape
                    case {'circular','PCHE'}
                        f1  = @(mC) hex_compute_area(HX,iL,'hC2',hC2,'mC',mC);
                    case 'cross-flow'
                        f1  = @(mC) hex_compute_Xflow(HX,iL,'hC2',hC2,'mC',mC);
                    otherwise
                        error('not implemented')
                end
                %plot_function(f1,mCmin,mCmax,20,30)
                %keyboard
                opt = optimset('TolX',(mCmax-mCmin)/1e12,'Display','notify');
                mC  = fzero(f1,[mCmin,mCmax],opt);
                
                % Store new mC value into HX and stateC structures
                HX.C(iL).mdot = mC;
                stateC.mdot   = mC;
                
            case 4
                % Set outlet conditions of hot fluid. Assume small effect
                % of pressure loss when computing outlet enthalpy from the
                % objective outlet temperature.
                TH1 = par;
                hH1 = RPN('PT_INPUTS',pH2,TH1,'H',fluidH);
                
                % Compute preliminary QMAX (cold outlet cannot be hotter
                % than hot inlet) and set boundaries accordingly
                QMAX0 = mC*(hC2_max - hC1);
                mHmax = QMAX0/(hH2 - hH1)*(1.01); %necessary to find root
                mHmin = mHmax/1e6;
                
                % Find value of mH for which computed area equals specified area
                switch HX.shape
                    case {'circular','PCHE'}
                        f1  = @(mH) hex_compute_area(HX,iL,'hH1',hH1,'mH',mH);
                    case 'cross-flow'
                        f1  = @(mH) hex_compute_Xflow(HX,iL,'hH1',hH1,'mH',mH);
                    otherwise
                        error('not implemented')
                end
                opt = optimset('TolX',(mHmax-mHmin)/1e12,'Display','notify');
                mH  = fzero(f1,[mHmin,mHmax],opt);
                
                % Store new mH value into HX and stateH structures
                HX.H(iL).mdot = mH;
                stateH.mdot   = mH;
                
            case 5
                % Set outlet conditions of hot fluid. Assume small effect
                % of pressure loss when computing outlet enthalpy from the
                % objective outlet temperature.
                TH1 = par;
                hH1 = RPN('PT_INPUTS',pH2,TH1,'H',fluidH);
                
                % Compute total heat transfer and compute mCmin and mCmax
                % accordingly
                QT    = mH*(hH2-hH1);
                mCmin = QT/(hC2_max - hC1)*0.98;
                mCmax = mCmin*100; %necessary to find root
                
                % Find value of mC for which computed area equals specified area
                %keyboard
                switch HX.shape
                    case {'circular','PCHE'}
                        f1  = @(mC) hex_compute_area(HX,iL,'hH1',hH1,'mC',mC);
                    case 'cross-flow'
                        f1  = @(mC) hex_compute_Xflow(HX,iL,'hH1',hH1,'mC',mC);
                    otherwise
                        error('not implemented')
                end
                %plot_function(f1,mCmin,mCmax,10,15,'semilogx');
                opt = optimset('TolX',(mCmax-mCmin)/1e12,'Display','notify');
                mC = fzero(f1,[mCmin,mCmax],opt);
                
                % Store new mC value into HX and stateC structures
                HX.C(iL).mdot = mC;
                stateC.mdot   = mC;
                
        end
        
        % Obtain output parameters for converged solution
        switch mode
            case {0,1,2,4,5}
                switch HX.shape
                    case {'circular','PCHE'}
                        [~,HX] = hex_compute_area(HX,iL,'hH1',hH1);
                    case 'cross-flow'
                        [~,HX] = hex_compute_Xflow(HX,iL,'hH1',hH1);
                end
            case 3
                switch HX.shape
                    case {'circular','PCHE'}
                        [~,HX] = hex_compute_area(HX,iL,'hC2',hC2);
                    case 'cross-flow'
                        [~,HX] = hex_compute_Xflow(HX,iL,'hC2',hC2);
                end
        end
        
        % Extract outlet conditions
        hH1 = HX.H(iL).h(1);
        pH1 = HX.H(iL).p(1);
        hC2 = HX.C(iL).h(NX+1);
        pC2 = HX.C(iL).p(NX+1);
end

% Update states
stateH.h = hH1;
stateH.p = pH1;
stateH   = update_state(stateH,fluidH,2);
stateC.h = hC2;
stateC.p = pC2;
stateC   = update_state(stateC,fluidC,2);

% *** DELETE EVENTUALLY >>>
% Compute stages
% Entropy change
DsH         = stateH.s - sH2;
DsC         = stateC.s - sC1;
% Hot stream
stageH.Dh   = stateH.h - hH2;
stageH.sirr = (stateH.mdot*DsH + stateC.mdot*DsC)/stateH.mdot;
%stageH.sirr = (vpa(stateH.mdot*DsH) + vpa(stateC.mdot*DsC))/stateH.mdot;
stageH.q    = stageH.Dh;
stageH.w    = 0;
stageH.type = stage_type;
% Cold stream
stageC.Dh   = stateC.h - hC1;
stageC.sirr = (stateC.mdot*DsC + stateH.mdot*DsH)/stateC.mdot;
%stageC.sirr = (vpa(stateC.mdot*DsC) + vpa(stateH.mdot*DsH))/stateC.mdot;
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
        TH  = RPN('HmassP_INPUTS',hH,pH,'T',fluidH);
        
        % Compute temperature distribution of cold stream
        QS  = (hH - hH1)*mH; % cummulative heat transfer
        hC  = hC1 + QS/mC;
        TC  = RPN('HmassP_INPUTS',hC,pC,'T',fluidC);
        
    case 'hC2'
        
        % Set outlet enthalpy
        hC2 = hout;
        
        % Compute temperature distribution of cold stream
        hC  = linspace(hC1,hC2,n+1)';
        TC  = RPN('HmassP_INPUTS',hC,pC,'T',fluidC);
        
        % Compute temperature distribution of hot stream
        QS  = (hC - hC1)*mC; % cummulative heat transfer
        hH1 = hH2 - QS(n+1)/mH;
        hH  = hH1 + QS/mH;
        TH  = RPN('HmassP_INPUTS',hH,pH,'T',fluidH);
        
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%