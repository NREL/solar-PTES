function varargout = hex_compute_Xflow(HX,iL,mode,par,varargin)
%HEX_COMPUTE_XFLOW Solve the TQ and TA diagrams of a two-stream cross-flow
%heat exchanger.
%
%   There are three modes of operation, controled by the 'mode' argument.
%
%   In modes 'hH1' and 'hC2', the geometry of the heat exchanger is known
%   and the outlet enthalpy of one of the two fluids is specified (in the
%   'par' argument). The outlet enthalpy is used to compute the TQ diagram,
%   the properties of the fluids at each point and the required heat
%   transfer area of the heat exchanger. The computed heat transfer area is
%   compared to the reference heat transfer area and the difference between
%   the two is returned as output argument. This disagreement between areas
%   is employed by the iteration procedure of the calling function to find
%   the correct value of 'hH1' or 'hC2' for that geometry.
%
%   In mode 'Af', the geometry of the heat exchanger is not know yet, and
%   is to be created according to the performance objectives (effectiveness
%   and pressure loss). Because the effectiveness is specified, the TQ
%   diagram is known by the calling function. The calling function calls
%   HEX_COMPUTE_XFLOW while iterating over different 'Af1' values. The
%   value of 'Af1' (stored in the 'par' argument) is used to compute the
%   mass flux, the Re, St and Cf values at each point and generate a
%   geometry that satisfies the given TQ diagram. Then, the predicted
%   pressure losses are compared to the specified pressure loss and the
%   difference is returned as output argument (to be used by the iteration
%   procedure of the calling function) untill the correct value of 'Af1' is
%   found.
%
%   The optional inputs, "varargin", contain 'mC', 'mH' and the logical
%   argument "visualise". If "visualise" is set to true, the function
%   pauses at the end of the internal iteration procedure and plots the
%   result.
%
%   USAGE, E.G.:
%   hex_compute_Xflow(HX,iL,'hH1',hH1,true)
%   hex_compute_Xflow(HX,iL,'hH1',hH1,'mH',mH)
%   hex_compute_Xflow(HX,iL,'hH1',hH1,'mC',mC,true)
%   hex_compute_Xflow(HX,iL,'hC2',hC2,'mC',mC)
%   hex_compute_Xflow(HX,iL,'Af',Af1,0)

% Check that this is the correct algorithm to be called
switch HX.shape
    case 'cross-flow'
    case {'circular','PCHE'}
        error(['The hex_compute_area function must be used for',...
            ' counter-flow heat exchangers.'])
    otherwise
        error('HX.shape not recognised')
end

% Select mode and assing value of par
switch mode
    case 'hH1'
        hH1 = par;
    case 'hC2'
        hC2 = par;
    case 'Af'
        HX.Af1 = par;
    otherwise
        error('not implemented')
end

% Check and assign variable inputs
visualise = 0;
switch length(varargin)
    case 0
        
    case 1
        visualise = varargin{1};
        
    case {2,3}
        switch varargin{1}
            case 'mC'
                HX.C(iL).mdot = varargin{2};
            case 'mH'
                HX.H(iL).mdot = varargin{2};
            otherwise
                error('not implemented')
        end
        if length(varargin)==3
            visualise = varargin{3};
        end
        
    otherwise
        error('not implemented')
end

% Extract parameters
H   = HX.H(iL);
C   = HX.C(iL);
mH  = H.mdot;
mC  = C.mdot;
NX  = HX.NX;
hH2 = H.hin;
hC1 = C.hin;

% Compute enthalpy arrays from guess outlet values
switch mode
    case 'hH1'
        H.h = linspace(hH1,hH2,NX+1)';
        C.h = hC1 + (H.h - hH1)*mH/mC;
    case 'hC2'
        C.h = linspace(hC1,hC2,NX+1)';
        H.h = hH2 - (hC2 - C.h)*mC/mH;
end

% Set pressure arrays (guess)
if isempty(H.p)
    H.p = ones(NX+1,1)*H.pin;
    C.p = ones(NX+1,1)*C.pin;
end

% Heat fluxes
switch mode
    case {'hH1','hC2'}
        H.D  = HX.D2;
        H.Af = HX.Af2;
        H.A  = HX.A2;
        C.D  = HX.D1;
        C.Af = HX.Af1;
        C.A  = HX.A1;
        q  = (H.h(NX+1)-H.h(1))/H.A;
        %fprintf(1,'\nHX.A1 = %.3e, HX.A2 = %.3e\n',HX.A1,HX.A2);
        %keyboard
    case 'Af'
        q = 100; %guess
end
H.q = q;
C.q = q*C.A/H.A;

% Convention for heat rejection unit: fluidC (air) flows on channel 1,
% while fluidH (working fluid) flows on channel 2.
%
% In the 'Af' mode, define cross-flow geometry according to description of
% surface 8.0_3/8T from Kays&London and Nellis&Klein.
%
% Here we assume that the air side has a square shape, with equal width and
% height (W1 = H1). The length of the air side is equal to the width of the
% working fluid side (L = W2). The length of the working fluid side is
% equal to the width of the air side (L2 = W1). The height of the air side
% is qual to the height of the working fluid side (H1 = H2).
%
% Geometry-specific characteristics:
% The air side hydraulic diameter is      DC =  3.63 mm
% The working fluid hydraulic diameter is DH = 10.21 mm
% The air's free-flow area to frontal area ratio is sigma1 = 0.534.
% The definition of sigma means that: W1 = H1 = sqrt(Af1/sigma1);
% The working fluid's free-flow area ratio was computed as sigma2 = 0.146

% Check that the cold stream is environmental air
cond1 = any(strcmp(C.name,{'Air','Nitrogen'}));      % check name
cond2 = C.pin <= 1.2e5 && C.pin >= 0.8e5;            % check pressure range
C_Tin = RPN('HmassP_INPUTS',C.hin,C.pin,'T',C);
cond3 = C_Tin <= 340 && C_Tin >= 240;                % check temperature range
if ~all([cond1,cond2,cond3])
    error('It has not been possible to verify that fluidC is ambient air.')
end

% Set shape of each channel
H.shape = 'circular';
C.shape = 'cross-flow';

switch mode
    case 'Af'
        % Set fixed geometrical parameters
        DC =  3.63e-3;
        DH = 10.21e-3;
        sigma1 = 0.534;
        sigma2 = 0.146;
        
        HX.D1  = DC;            % Hydraulic diameter (air side)
        HX.D2  = DH;            % Hydraulic diameter (working fluid side)
        
        % Geometry initial guess
        L1  = 0.1; % heat exchanger length (air side)
        WHr = 1;   % Width-to-height ratio (air side)
end

% UPDATE PROPERTIES
% Cold stream
C   = stream_update(C,2);
% Hot stream
H   = stream_update(H,2);

% Check if the HEX is a condenser
if any(H.x >= 0.0 & H.x <= 1.0)
    Lcondenser = true;
else
    Lcondenser = false;
end

% Heat capacity rates
CpC  = ( C.h(NX+1)-C.h(1) )/( C.T(NX+1)-C.T(1) );
CpH  = ( H.h(NX+1)-H.h(1) )/( H.T(NX+1)-H.T(1) );
%{
if Lcondenser
    Cmax = mH*CpH;
    Cmin = mC*CpC;
else
    Cmax = max([mC*CpC,mH*CpH]);
    Cmin = min([mC*CpC,mH*CpH]);
end
%}
%%{
Cmax = max([mC*CpC,mH*CpH]);
Cmin = min([mC*CpC,mH*CpH]);
%%}
Cr   = Cmin/Cmax;

% Cummulative and total heat transfer
QS = mH*(H.h - H.h(1));
QT = mH*(H.h(NX+1) - H.h(1));

% Compute heat exchanger effectiveness for given hH1
hH1min = RPN('PT_INPUTS',H.pin,C.T(1),'H',H);
hC2max = RPN('PT_INPUTS',C.pin,H.T(NX+1),'H',C);
%%{
if Lcondenser
    QMAX = mC*(hC2max - hC1);
    %QMAX = mH*(hH2 - hH1min);
else
    QMAX   = min([mC*(hC2max - hC1),mH*(hH2 - hH1min)]);
end
%%}
%QMAX   = min([mC*(hC2max - hC1),mH*(hH2 - hH1min)]);
eff    = QT/QMAX;

% Create array to check convergence. First element is computed NTU.
% Later come the pressure points along each stream
CON_0 = [0; H.p; C.p]; % initial value
NI  = 50;
RES = zeros(1,NI); % residuals
TOL = 1e-4;
impossible_p = false; %indicates impossible situations (regarding Dp loss)
for iI = 1:NI
        
    switch mode
        case 'Af'
            % Compute NTU using analytical solutions
            if Lcondenser || Cr < 1e-3
                fun =@(NTU) 1 - exp(-NTU) - eff;
            else
                fun =@(NTU) 1 - exp( (1/Cr)*NTU.^0.22.*(exp(-Cr*NTU.^0.78) - 1) ) - eff;
            end
            NTU_needed = fzero(fun,[0 1e3]);
            
            % Main geometric parameters
            W1     = WHr*sqrt(HX.Af1/sigma1);   % Width  (air side)
            H1     = sqrt(HX.Af1/sigma1)/WHr;   % Height (air side)
            L2     = W1;                        % Length (working fluid side)
            H2     = H1;                        % Height (working fluid side)
            HX.Af2 = L1*H2*sigma2;              % Free-flow area (working fluid side)
            HX.A1  = 4*HX.Af1*L1/HX.D1;         % Heat transfer area (air side)
            HX.A2  = 4*HX.Af2*L2/HX.D2;         % Heat transfer area (wf side)
            
            % Store information inside each stream structure
            H.D  = HX.D2;
            H.Af = HX.Af2;
            H.A  = HX.A2;
            C.D  = HX.D1;
            C.Af = HX.Af1;
            C.A  = HX.A1;
    end
    
    % Compute mass fluxes
    H.G = mH/H.Af;
    C.G = mC/C.Af;
    
    % COMPUTE HEAT TRANSFER AND FRICTION COEFFICIENTS
    % Cold stream
    [C] = developed_flow(C,'heating');
    % Hot stream
    [H] = developed_flow(H,'cooling');
    
    % Compute average heat transfer coefficients
    htC = mean(C.ht);
    htH = mean(H.ht);
    UA  = 1/(1/(htC*C.A) + 1/(htH*H.A));
    
    % Compute number of transfer units
    NTU = UA/Cmin;
    
    switch mode
        case 'Af'
            % Compute heat transfer areas and channel length
            % By definition of hydraulic diameter, D = 4*Af/P = 4*Af*L/A
            AC = C.A*NTU_needed/NTU;
            AH = H.A*NTU_needed/NTU;
            L1 = C.D*AC/(4*C.Af);
            %L1 = max([L1,0.1]);
            
            % Heat flux
            H.q = QT/AH;
            C.q = q*C.A/H.A;
            
        case {'hH1','hC2'}
            if Lcondenser || Cr < 1e-3
                eff_an = 1 - exp(-NTU);
            else
                eff_an = 1 - exp( (1/Cr)*NTU.^0.22.*(exp(-Cr*NTU.^0.78) - 1) );
            end
            L1 = C.D*C.A/(4*C.Af);
            L2 = H.D*H.A/(4*H.Af);
    end
    
    % COMPUTE PRESSURE PROFILES
    % Compute arrays of pressure loss due to flow friction
    dL2 = L2/NX*ones(NX,1);
    dL1 = L1/NX*ones(NX,1);
    Dpf_H = 0.5*(H.dpdL(1:NX) + H.dpdL(2:NX+1)).*dL2;
    Dpf_C = 0.5*(C.dpdL(1:NX) + C.dpdL(2:NX+1)).*dL1;
    
    % In a situation where H.h(NX+1)==H.h(1) (i.e. no heat exchange), the
    % code above results in dAC=0, AC=0 and Dp_H=NaN, Dp_C=NaN. If so, set
    % Dp_H and Dp_C to zero and proceed.
    if any(isnan([Dpf_H;Dpf_C]))
        if abs(H.h(NX+1)-H.h(1)) < 10*eps(H.h(1))
            Dpf_H = zeros(size(Dpf_H));
            Dpf_C = zeros(size(Dpf_C));
        else
            error('NaN value found in hex_func!')
        end
    end
    
    % Compute arrays of pressure loss due to flow acceleration
    %Dpa_H = -( H.G.^2.*H.v(1:NX) - H.G.^2.*H.v(2:NX+1) );
    %Dpa_C = -( C.G.^2.*C.v(2:NX+1) - C.G.^2.*C.v(1:NX) );
    Dpa_H = 0; %ignore flow acceleration
    Dpa_C = 0; %ignore flow acceleration
    
    % Update pressure profiles. Assume a linear profile (to avoid computing
    % a slow 'for loop') and limit max pressure loss.
    DppH = abs(sum(Dpf_H + Dpa_H))/H.pin;
    DppC = abs(sum(Dpf_C + Dpa_C))/C.pin;
    switch mode
        case 'Af'
            if DppH > 0.02
                WHr  = WHr/max([2,sqrt(DppH)]);
                DppH = 0.02;
            end
            if DppC > 0.2
                DppC = 0.2;
                impossible_p = true;
                break
            end
        case {'hH1','hC2'}
            if DppH > 0.8
                DppH = 0.8;
            end
            if DppC > 0.8
                DppC = 0.8;
            end
    end
    H.p  = linspace(H.pin*(1-DppH),H.pin,NX+1)';
    C.p  = linspace(C.pin,C.pin*(1-DppC),NX+1)';
    
    % Update convergence array
    CON = [NTU; H.p; C.p];
    
    % Compute residual
    RES(iI) = max(abs((CON - CON_0)./CON));
    
    if RES(iI)>TOL
        if visualise
            fprintf(1,'\n iteration = %d, RES = %.6f',iI,RES(iI));
        end
        CON_0 = CON;
    else
        if visualise
            % Make TQ plot
            figure(10)
            plot(QS./QS(end),H.T,'r'); hold on;
            plot(QS./QS(end),C.T,'b'); hold off;
            xlabel('Cumulative heat transfer')
            ylabel('Temperature')
            legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],...
                [C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],...
                'Location','Best')
            
            % Make pQ plot
            figure(11)
            plot(QS./QS(end),H.p/H.pin,'r-'); hold on
            plot(QS./QS(end),C.p/C.pin,'b-'); hold off
            ylim([0.90 1])
            xlabel('Cummulative heat transfer')
            ylabel('Relative pressure, p/p0')
            legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],...
                [C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],...
                'Location','Best')
            
            formatSpec = ['\n\n*** Successful convergence in ',...
                'hex_compute_Xflow inner loop after %d iterations***\n'];
            fprintf(formatSpec,iI);
            keyboard
        end
        break
    end
    
end
if all([iI>=NI,RES(iI)>TOL])
    figure()
    semilogy(1:iI,RES(1:iI))
    xlabel('Iteration')
    ylabel('Convergence residual')
    error('Convergence not reached after %d iterations***\n',iI);
end

switch mode
    case {'hH1','hC2'}
        % If the value of 'solution' is negative, it means that the
        % computed area is larger than the actual area (heat exchanger too
        % small to achieve selected operating conditions). If the value of
        % 'solution' is positive, the computed area is smaller than the
        % actual area (heat exchanger too large for selected operating
        % conditions).
        solution = eff_an - eff;
        %formatSpec = 'hH1 = %9.4g, hC2 = %9.4g, mC = %6.0f, mH = %4.1f, UA = %8.3g, NTU = %5.2f, eff_num = %5.3f, eff_an = %5.3f, solution = %9.3g, iter = %2d\n';
        %fprintf(1,formatSpec,H.h(1),C.h(end),mC,mH,UA,NTU,eff,eff_an,solution,iI)
        %keyboard
        
    case 'Af'
        % Compare ploss and max_ploss and return the difference. If the
        % value of 'solution' is positive, it means that the pressure loss
        % in one of the two streams is larger than the set objective
        % (indicating that the flow area might be too small), and
        % vice-versa.
        if impossible_p
            max_ploss = 1.0;
        else
            max_ploss = DppC;%max([DppH,DppC]);
        end
        solution = max_ploss - HX.ploss;
        %formatSpec = ['Af1 = %8.3g, L1 = %6.3f, WHr = %5.3f, UA = %8.3g, ',...
        %    'NTU = %5.2f, NTU2 = %5.2f, eff = %5.3f, htC*AC = %8.3g, htH*AH = %8.3g, solution = %9.3g, iter = %2d\n'];
        %fprintf(1,formatSpec,HX.Af1,L1,WHr,UA,NTU,NTU_needed,eff,htC*AC,htH*AH,solution,iI)
        %keyboard
end

%{
% Make TQ plot
figure(10)
plot(QS./QS(end),H.T,'r'); hold on;
plot(QS./QS(end),C.T,'b'); hold off;
xlabel('Cumulative heat transfer')
ylabel('Temperature')
legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],...
    [C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],...
    'Location','Best')
keyboard
%}

if nargout == 1
    varargout{1} = solution;
else
    % Extract inlet and outlet conditions
    TH1 = H.T(1);
    TH2 = H.T(NX+1);
    TC1 = C.T(1);
    TC2 = C.T(NX+1);
    hH1 = H.h(1);
    hH2 = H.h(NX+1);
    hC1 = C.h(1);
    hC2 = C.h(NX+1);
    pH1 = H.p(1);
    pH2 = H.p(NX+1);
    pC1 = C.p(1);
    pC2 = C.p(NX+1);
    
    % Compute performance indicators
    H.Cp_mean = (hH2 - hH1)/(TH2-TH1);
    C.Cp_mean = (hC2 - hC1)/(TC2-TC1);
    Cmin  = min([mC*C.Cp_mean,mH*H.Cp_mean]);
    NTU   = UA/Cmin;
    DppH  = (pH2-pH1)/pH2;
    DppC  = (pC1-pC2)/pC1;
    DTmin = min(H.T - C.T);
    effDT = 1 - DTmin/(TH2 - TC1);
    
    % Save parameters into HX structure
    HX.H(iL) = H;
    HX.C(iL) = C;
    HX.Cmin(iL)  = Cmin;
    HX.NTU(iL)   = NTU;
    HX.UA(iL)    = UA ;
    HX.DppH(iL) = DppH;
    HX.DppC(iL) = DppC;
    HX.DTmin(iL) = DTmin;
    HX.effDT(iL) = effDT;
    HX.QS(iL,:) = QS';
    % Save parameters into HX structure
    HX.L1    = L1;
    HX.L2    = L2;
    switch mode
        case 'Af'
            HX.A1 = AC;
            HX.A2 = AH;
    end
    
    % If this is the first time that hex_func is called, save the initial
    % values UA0, NTU0 and LMTD0. Also save iL0, corresponding to the load
    % period in which it is called for the first time.
    if isempty(HX.UA0)
        HX.UA0   = HX.UA(iL) ;
        HX.NTU0  = HX.NTU(iL) ;
        HX.iL0   = iL;
    end
    
    % Set varargout fields
    varargout{1} = solution;
    varargout{2} = HX;
end

end