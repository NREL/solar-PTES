function [HX] = set_hex_geom2(HX, iL, fluidH, iH, fluidC, iC, mode, par)
%SET_HEX Determine the geometry of a heat exchanger
%
%   There are two scenarios when this is required:
%
%   1) The heat exchangers are defined in 'geom' mode so geometry is
%   required. The geometry is obtained to satisfy the performance
%   objectives set by DT and ploss.
%
%   2) The heat exchanger was already solved using the 'eff' or 'DT' modes
%   but the geometry should now be estimated for economic calculations.
%
%   Note that SET_HEX has the same inputs as HEX_FUNC, which is used in the
%   first scenario.
%
%   The sizing procedure follows the steps below:
%
%   (a) To simplify things, assume that both sides have the same hydraulic
%   diameter and cross-sectional area (this leads to slightly sub-optimal
%   designs and conservative cost estimates, but seems necessary to speed
%   things up at run-time while obtaining geometries that accurately match
%   the performance objectives)
%
%   (b) For a given hydraulic diameter, iterate over different values of
%   Af.
%
%   (c) For each value of Af, obtain in the following order: G1 and G2, Re1
%   and Re2, St2 and St2, Cf1 and Cf2, and the overall heat transfer
%   coefficient at each heat exchanger section.
%
%   (d) Compute required heat transfer area and total pressures loss.


% This function should only be called if the geometry has not been set.
if HX.Lgeom_set
    error('The heat exchanger geometry has already been set!')
end

% Set geometry according to scenario
switch HX.model
    case 'geom'
        
        for i0=1:2
            % Use hex_func to obtain temperature profiles. Temporarily set
            % HX.model to 'DT' in order to achieve this.
            HX.model = 'DT';
            [HX,~] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par);
            HX.model = 'geom';
            
            % Employ the compute_pressure function to determine for which
            % value of Af the pressure loss is the same as the objective.
            % Check whether f1 changes sign over interval. If it doesn't,
            % choose Af that is closest to zero.
            f1 = @(Af) compute_pressure(HX,iL,Af,0);
            Afmin = 1e-3;
            Afmax = 1e+3;
            f1min = f1(Afmin) ;
            f1max = f1(Afmax) ;
            if f1min*f1max >= 0
                if abs(f1min) < abs(f1max)
                    Af0 = Afmin ;
                else
                    Af0 = Afmax ;
                end
            else
                opt = optimset('TolX',Afmin/1e3,'Display','notify');
                Af0 = fzero(f1,[Afmin,Afmax],opt);
            end
            [~,HX] = f1(Af0);
            
            % Update the 'design' values of ploss (either plossH0 = ploss
            % or plossC0 = ploss, with the other one being lower) and
            % repeat design process once.
            HX.plossH0 = HX.DppH(iL);
            HX.plossC0 = HX.DppC(iL);
        end
        
    case {'eff','DT'}
        
        % Employ the compute_pressure function to determine for which
        % value of Af the pressure loss is the same as the objective.
        % Check whether f1 changes sign over interval. If it doesn't,
        % choose Af that is closest to zero.
        f1 = @(Af) compute_pressure(HX,iL,Af,0);
        Afmin = 1e-3;
        Afmax = 1e+3;
        f1min = f1(Afmin) ;
        f1max = f1(Afmax) ;
        if f1min*f1max >= 0
            if abs(f1min) < abs(f1max)
                Af0 = Afmin ;
            else
                Af0 = Afmax ;
            end
        else
            opt = optimset('TolX',Afmin/1e4,'Display','notify');
            Af0 = fzero(f1,[Afmin,Afmax],opt);
        end
        [~,HX] = f1(Af0);
        
    otherwise
        error('not implemented')
end

% Note that geometry is now set
HX.Lgeom_set = true;

end


%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = compute_pressure(HX,iL,Af,varargin)

% Check and assign variable inputs
if isempty(varargin)
    visualise = 0;
elseif length(varargin)==1
    visualise = varargin{1};
else
    error('not implemented')
end

% Extract parameters
H   = HX.H(iL);
C   = HX.C(iL);
mH  = H.mdot;
mC  = C.mdot;
NX  = HX.NX;

% Assume equal hydraulic diameters and flow areas
HX.D2  = HX.D1;
HX.Af1 = Af/2;
HX.Af2 = Af/2;
HX.AfT = HX.Af1+HX.Af2;
HX.AfR = HX.Af2/HX.Af1;

% Determine which fluid flows inside the pipes and set hydraulic diameters
% and flow areas
if H.pin > C.pin
    H.D  = HX.D1;
    H.Af = HX.Af1;
    C.D  = HX.D2;
    C.Af = HX.Af2;
else
    H.D  = HX.D2;
    H.Af = HX.Af2;
    C.D  = HX.D1;
    C.Af = HX.Af1;
end

% Compute mass fluxes
H.G = mH/H.Af;
C.G = mC/C.Af;

% Set pressure arrays (this is approximate because actual profile is not
% exactly linear, but the most important are the inlet and outlet points).
H.p = ones(NX+1,1)*H.pin;
C.p = ones(NX+1,1)*C.pin;

% Create array to check convergence. First element is computed heat
% transfer area. Later come the pressure points along each stream
CON_0 = [0; H.p; C.p]; % initial value
NI    = 50;
RES   = zeros(1,NI); % residuals
TOL   = 1e-3;
impossible = false;
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
    
    % COMPUTE HEAT TRANSFER COEFFICIENTS
    % Cold stream
    C.Re = C.D*C.G./C.mu;
    [C.Cf,C.St] = developed_flow(C.Re,C.Pr,HX.shape);
    C.ht  = C.G*C.Cp.*C.St;
    % Hot stream
    H.Re = H.D*H.G./H.mu;
    [H.Cf,H.St] = developed_flow(H.Re,H.Pr,HX.shape);
    H.ht  = H.G*H.Cp.*H.St;
    % Local (total) heat transfer coefficient.
    % Neglects wall thermal resistance and axial conduction, and assumes
    % that the heat transfer area is the same for the hot and cold sides.
    Ul  = 1./(1./H.ht + 1./C.ht);
    Ul_AV = 0.5*(Ul(1:NX) + Ul(2:NX+1));
    
    % COMPUTE HEAT TRANSFER AREA and CUMMULATIVE HEAT TRANSFER
    dQ = (H.h(2:NX+1) - H.h(1:NX))*mH;
    QS = mH*(H.h - H.h(1));
    dA = dQ./(Ul_AV.*DT_AV);
    A  = sum(dA);
    
    % COMPUTE HEAT EXCHANGER LENGTH
    % By definition of hydraulic diameter, D = 4*Af/P = 4*Af*L/A
    % So, L = D*A/(4*Af)
    HX.L = C.D*A/(4*C.Af);
    
    % Restrict designs to L<=50m
    %if HX.L > 50
    %    impossible = true;
    %    break
    %end
    
    % COMPUTE PRESSURE PROFILES
    % Create averaged arrays of Cf and v
    Cf_H = 0.5*(H.Cf(1:NX) + H.Cf(2:NX+1));
    v_H  = 0.5*(H.v(1:NX)  + H.v(2:NX+1));
    Cf_C = 0.5*(C.Cf(1:NX) + C.Cf(2:NX+1));
    v_C  = 0.5*(C.v(1:NX)  + C.v(2:NX+1));
    % Obtain dL from dAC and AC
    dL = dA/A*HX.L;
    % Compute arrays of pressure loss
    Dp_H = - 2*H.G^2*Cf_H.*v_H.*dL./H.D;
    Dp_C = - 2*C.G^2*Cf_C.*v_C.*dL./C.D;
    % Update pressure profiles. Assume a linear profile (to avoid computing
    % a slow 'for loop') and limit max pressure loss to 80%.
    DppH = abs(sum(Dp_H))/H.pin;
    DppC = abs(sum(Dp_C))/C.pin;
    if DppH > 0.8
        impossible = true;
        break
    end
    if DppC > 0.8
        impossible = true;
        break
    end
    H.p  = linspace(H.pin*(1-DppH),H.pin,NX+1)';
    C.p  = linspace(C.pin,C.pin*(1-DppC),NX+1)';
    
    % Update convergence array
    CON = [A; H.p; C.p]; % initial value
    
    % Compute residual
    RES(iI) = max(abs((CON - CON_0)./CON));
    
    if all([RES(iI)>TOL])
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
            legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],...
                [C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')
            
            figure(11)
            plot(QS./QS(end),H.p/H.pin,'r-'); hold on
            plot(QS./QS(end),C.p/C.pin,'b-'); hold off
            ylim([0.90 1])
            xlabel('Cummulative heat transfer')
            ylabel('Relative pressure, p/p0')
            legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],...
                [C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')
            
            fprintf(1,'\n\n*** Successful convergence after %d iterations***\n',iI);
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

if impossible
    max_ploss = 1.0;
else
    max_ploss = max([DppH,DppC]);
end

% Compare ploss and max_ploss and return the difference. If the value of
% 'solution' is positive, it means that the pressure loss in one of the two
% streams is larger than the set objective (indicating that the flow area
% might be too small), and vice-versa.
solution = max_ploss - HX.ploss;

if nargout == 1
    varargout{1} = solution;
else
    % Compute cummulative heat transfer area
    AS = zeros(size(QS));
    for i=1:(length(AS)-1)
        AS(i+1) = AS(i) + dA(i);
    end
    % Save parameters into HX structure
    HX.H(iL) = H;
    HX.C(iL) = C;
    HX.DppH(iL) = DppH;
    HX.DppC(iL) = DppC;
    HX.AS    = AS;
    HX.A1    = A;
    HX.A2    = A;
    varargout{1} = solution;
    varargout{2} = HX;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%