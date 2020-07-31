function varargout = hex_compute_area(HX,iL,mode,par,varargin)
%HEX_COMPUTE_AREA Solve the TQ and TA diagrams of a two-stream counter-flow
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
%   HEX_COMPUTE_AREA while iterating over different 'Af' values. The value
%   of 'Af' (stored in the 'par' argument) is used to compute the mass
%   flux, the Re, St and Cf values at each point and generate a geometry
%   that satisfies the given TQ diagram. Then, the predicted pressure
%   losses are compared to the specified pressure loss and the difference
%   is returned as output argument (to be used by the iteration procedure
%   of the calling function) untill the correct value of 'Af' is found.
%
%   The optional inputs, "varargin", contain 'mC', 'mH' and the logical
%   argument "visualise". If "visualise" is set to true, the compute_area
%   function pauses at the end of the internal iteration procedure and
%   plots the result.
%
%   USAGE, E.G.:
%   hex_compute_area(HX,iL,'hH1',hH1,true)
%   hex_compute_area(HX,iL,'hH1',hH1,'mH',mH)
%   hex_compute_area(HX,iL,'hH1',hH1,'mC',mC,true)
%   hex_compute_area(HX,iL,'hC2',hC2,'mC',mC)
%   hex_compute_area(HX,iL,'Af',Af,0)

% Check that this is the correct algorithm to be called
switch HX.shape
    case {'circular','PCHE'}
    case 'cross-flow'
        error(['The hex_compute_Xflow function must be used for',...
            ' cross-flow heat exchangers.'])
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
        Af = par;
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

% In the 'Af' mode, assume equal hydraulic diameters and flow areas
switch mode
    case 'Af'
        HX.D2  = HX.D1;
        HX.Af1 = Af/2;
        HX.Af2 = Af/2;
        HX.AfT = HX.Af1+HX.Af2;
        HX.AfR = HX.Af2/HX.Af1;
end

% Determine which fluid flows inside the pipes and set hydraulic diameters,
% flow areas and heat transfer areas (the latter are still unknown in 'Af'
% mode)
if H.pin > C.pin
    H.D  = HX.D1;
    H.Af = HX.Af1;
    H.A  = HX.A1;
    C.D  = HX.D2;
    C.Af = HX.Af2;
    C.A  = HX.A2;
else
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

% Compute enthalpy arrays from guess outlet values
switch mode
    case 'hH1'
        H.h = linspace(hH1,hH2,NX+1)';
        C.h = hC1 + (H.h - hH1)*mH/mC;
    case 'hC2'
        C.h = linspace(hC1,hC2,NX+1)';
        H.h = hH2 - (hC2 - C.h)*mC/mH;
end

% Set initial conditions for iteration procedure
% Pressures
H.p = ones(NX+1,1)*H.pin;
C.p = ones(NX+1,1)*C.pin;
% Heat flux
switch mode
    case {'hH1','hC2'}
        q = ones(NX+1,1)*(C.h(NX+1)-C.h(1))/C.A;
    case 'Af'
        q = ones(NX+1,1)*100; %guess
end

% Create array to check convergence. First element is computed heat
% transfer area. Later come the pressure points along each stream
CON_0 = [0; H.p; C.p]; % initial value
NI    = 50;
RES   = zeros(1,NI); % residuals
TOL   = 1e-3;
impossible_h = false; %indicates impossible situations (regarding h guess)
impossible_p = false; %indicates impossible situations (regarding Dp loss)
for iI = 1:NI
    
    % UPDATE PROPERTIES
    % Cold stream
    C   = stream_update(C,2);
    C.q = q;
    % Hot stream
    H   = stream_update(H,2);
    H.q = q;
    
    % COMPUTE AVERAGED TEMPERATURE ARRAYS
    H.T_AV = 0.5*(H.T(1:NX) + H.T(2:NX+1));
    C.T_AV = 0.5*(C.T(1:NX) + C.T(2:NX+1));
    DT_AV  = H.T_AV - C.T_AV;
    
    % Break loop if H.T < C.T at any point
    if any(H.T <= C.T)
        impossible_h = true;
        A = Inf;
        break
    end
    
    % COMPUTE HEAT TRANSFER COEFFICIENTS
    % Cold stream
    [C] = developed_flow(C,'heating');
    % Hot stream
    [H] = developed_flow(H,'cooling');
    % Local (total) heat transfer coefficient.
    % Neglects wall thermal resistance and axial conduction, and assumes
    % that the heat transfer area is the same for the hot and cold sides.
    Ul    = 1./(1./H.ht + 1./C.ht);
    Ul_AV = 0.5*(Ul(1:NX) + Ul(2:NX+1));
    
    % COMPUTE HEAT TRANSFER AREA and CUMMULATIVE HEAT TRANSFER
    dQ = (H.h(2:NX+1) - H.h(1:NX))*mH;
    QS = mH*(H.h - H.h(1));
    dA = dQ./(Ul_AV.*DT_AV);
    A  = sum(dA);
    
    % Update heat flux
    q(1:NX) = dQ./dA;
    q(NX+1) = q(NX);
    
    switch mode
        case 'Af'
            % COMPUTE HEAT EXCHANGER LENGTH
            % By definition of hydraulic diameter, D = 4*Af/P = 4*Af*L/A
            L = C.D*A/(4*C.Af);
            HX.L1 = L;
            HX.L2 = L;
    end
    
    % COMPUTE PRESSURE PROFILES
    % Obtain dL from dAC and AC
    dL   = dA/A*HX.L1;
    % Compute cummulative sums of the dL array (used to account for entry
    % effects)
    LSH  = cumsum(dL,'reverse');
    H.LS = [LSH(1);LSH];
    LSC  = cumsum(dL);
    C.LS = [LSC;LSC(end)];
    % Compute arrays of pressure loss due to flow friction
    Dpf_H = 0.5*(H.dpdL(1:NX) + H.dpdL(2:NX+1)).*dL;
    Dpf_C = 0.5*(C.dpdL(1:NX) + C.dpdL(2:NX+1)).*dL;
    
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
    Dpa_H = -( H.G.^2.*H.v(1:NX) - H.G.^2.*H.v(2:NX+1) );
    Dpa_C = -( C.G.^2.*C.v(2:NX+1) - C.G.^2.*C.v(1:NX) );
    %Dpa_H = 0;
    %Dpa_C = 0;
    
    % Update pressure profiles. Assume a linear profile (to avoid computing
    % a slow 'for loop') and limit max pressure loss to 80%.
    DppH = abs(sum(Dpf_H + Dpa_H))/H.pin;
    DppC = abs(sum(Dpf_C + Dpa_C))/C.pin;
    if DppH > 0.8
        DppH = 0.8;
        impossible_p = true;
        break
    end
    if DppC > 0.8
        DppC = 0.8;
        impossible_p = true;
        break
    end
    H.p  = linspace(H.pin*(1-DppH),H.pin,NX+1)';
    C.p  = linspace(C.pin,C.pin*(1-DppC),NX+1)';
    
    % Update convergence array
    CON = [A; H.p; C.p]; % initial value
    
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
                'hex_compute_area inner loop after %d iterations***\n'];
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
        solution = C.A - A;
        % Control physically impossible solutions
        if impossible_h
            solution = - C.A;
        end
        
    case 'Af'
        % Compare ploss and max_ploss and return the difference. If the
        % value of 'solution' is positive, it means that the pressure loss
        % in one of the two streams is larger than the set objective
        % (indicating that the flow area might be too small), and
        % vice-versa.
        if impossible_p
            max_ploss = 1.0;
        else
            max_ploss = max([DppH,DppC]);
        end
        solution = max_ploss - HX.ploss;
end

if nargout == 1
    varargout{1} = solution;
    
else
    % Compute cummulative heat transfer area
    AS = zeros(size(QS));
    for i=1:(length(AS)-1)
        AS(i+1) = AS(i) + dA(i);
    end
    
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
    
    % Compute several performance indicators
    H.Cp_mean = (hH2 - hH1)/(TH2-TH1);
    C.Cp_mean = (hC2 - hC1)/(TC2-TC1);
    Cmin  = min([mC*C.Cp_mean,mH*H.Cp_mean]);
    DTmin = min(H.T - C.T);
    effDT = 1 - DTmin/(TH2 - TC1);
    UA    = sum(dQ./DT_AV);
    NTU   = UA/Cmin;
    DppH  = (pH2-pH1)/pH2;
    DppC  = (pC1-pC2)/pC1;    
    dTa   = TH1 - TC1 ;
    dTb   = TH2 - TC2 ;
    LMTD  = (dTa - dTb) / log(dTa / dTb);
    
    % Save parameters into HX structure
    HX.H(iL) = H;
    HX.C(iL) = C;
    HX.Cmin(iL)  = Cmin;
    HX.NTU(iL)   = NTU;
    HX.UA(iL)    = UA ;
    HX.LMTD(iL)  = LMTD;
    HX.DTmin(iL) = DTmin;
    HX.effDT(iL) = effDT;
    HX.DppH(iL) = DppH;
    HX.DppC(iL) = DppC;
    HX.AS(iL,:) = AS';
    HX.QS(iL,:) = QS';
    HX.Ul(iL,:) = Ul';
    HX.dL(iL,:) = dL';
    HX.A1 = A;
    HX.A2 = A;
    
    % If this is the first time that hex_func is called, save the initial
    % values UA0, NTU0 and LMTD0. Also save iL0, corresponding to the load
    % period in which it is called for the first time.
    if isempty(HX.UA0)
        HX.UA0   = HX.UA(iL) ;
        HX.NTU0  = HX.NTU(iL) ;
        HX.LMTD0 = HX.LMTD(iL) ;
        HX.iL0   = iL;
    end
    
    % Set varargout fields
    varargout{1} = solution;
    varargout{2} = HX;
end

end