function [HX] = set_hex_geom(HX, iL, fluidH, iH, fluidC, iC, mode, par, NTUmin, ploss_max, D)
% Obtain the heat exchanger geometry based on the performance objectives
% specified by NTUmin and ploss_max.

% Things to be improved:
% (1) Accurate computation of Nussel number for tube bundle (see
% developed_flow function) according to square/trilateral lattice
% configuration.
% (2) Minimum D2 value according to lattice configuration.
% (3) Introduction of tube thickness if metal volume is wanted

if isempty(HX.UA0)
    % Use hex_func to obtain the thermal profiles for the specified NTU.
    % First, run with 'eff'=1.0 to obtain approximate Cmin.
    model0 = HX.model;
    eff0   = HX.eff;
    UA0    = HX.UA;
    HX.model = 'eff';
    HX.eff   = 1.0;
    HX.ploss = ploss_max;
    [HX,~,~,~,~] = hex_func(HX,iL,fluidH,iH,fluidC,iC,mode,par);
    % Now use Cmin to compute UA and update profiles.
    HX.model = 'UA';
    HX.UA    = HX.Cmin*NTUmin;
    [HX,~,~,~,~] = hex_func(HX,iL,fluidH,iH,fluidC,iC,mode,par);
    % Reset parameters
    HX.model = model0; HX.eff=eff0; HX.UA=UA0;
    %plot_hex(HX,1,'C')
end

% Declare the two fluid streams
SH = stream;
SC = stream;

% Import properties from HX structure
SH.h = HX.H(iL).h;
SH.p = HX.H(iL).pin*ones(size(SH.h));
SH.mdot = HX.H(iL).mdot;
SH = stream_update(fluidH,SH,1);
SC.h = HX.C(iL).h;
SC.p = HX.C(iL).pin*ones(size(SC.h));
SC.mdot = HX.C(iL).mdot;
%keyboard
SC = stream_update(fluidC,SC,1);
%keyboard

% Determine which one is the "weak" stream (the one likely to present
% higher thermo-hydraulic losses for a given value of G). Sizing will start
% from that stream. Thermo-hydraulic losses are found to be proportional to
% the factor v/p. More specifically, (Dp/p)/Ntu = G^2/2 * Cf/St * v/p, with
% Cf/St ~ 2*Pr^2/3.
LF_H = mean(SH.v./SH.p.*SH.Pr.^(2/3)); %"loss factor" hot stream
LF_C = mean(SC.v./SC.p.*SC.Pr.^(2/3)); %"loss factor" cold stream
if LF_H > LF_C
    % Hot stream is "weak" stream. Cold stream is "strong" stream
    SW = SH;
    SS = SC;
else
    % Cold stream is "weak" stream. Hot stream is "strong" stream
    SW = SC;
    SS = SH;
end

% Set the minimum number of transfer units that each stream should have to
% obtain the specified overall NTU
Ntu_min = 2*NTUmin;


%%% DESIGN FIRST STREAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set pre-specified hydraulic diameter
SW.D = D;

% Obtain the maximum mass flux that satisfies the Ntu and ploss
% requirements. Make an initial guess assuming Cf/St=2*Pr^2/3;
SW.G = sqrt(mean( 2*SW.p./SW.v * ploss_max/Ntu_min .* 1./(2*SW.Pr.^(2/3)) ));

max_iter = 100;
tol = 1e-6;
for i=1:max_iter
    % Keep track of initial value (or value from previous iteration)
    G1_0 = SW.G;
    
    % Compute Cf and St
    SW.Re = SW.G*SW.D./SW.mu;
    [SW.Cf, SW.St] = developed_flow(SW.Re,SW.Pr,'circular');
    
    % Update the value of G
    SW.G = sqrt(mean( 2*SW.p./SW.v * ploss_max/Ntu_min .* SW.St./SW.Cf ));
    
    % Check convergence
    condition = abs((SW.G - G1_0)/G1_0*100) < tol;
    if condition
        % Converged
        break
    end
end
if all([i>=max_iter,~condition])
    error('Convergence not found')
end

% Compute heat transfer area, tube length and flow area
SW.A  = SW.mdot*Ntu_min/(SW.G*mean(SW.St));
SW.L  = SW.D*SW.A*SW.G/(4*SW.mdot);
SW.Af = SW.mdot/SW.G;


%%% DESIGN SECOND STREAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set L2=L1 (i.e. counter-flow) and A2=A1 (i.e. thin tubes)
SS.L = SW.L;
SS.A = SW.A;

% Assume an initial value of D2 and update until it satisfies
% thermo-hydraulic performance requirements
% SS.D = SW.D;

% Find the optimal hydraulic diameter that minimises losses on the second
% stream
fun = @(D) stream_losses(D,SS);
xmin = SW.D/100;
xmax = SW.D*100;
%plot_function(fun,xmin,xmax,100,30,'loglog')
%options = optimset('Display','notify','TolX',0.01*xmin,...
%    'FunValCheck','on','AlwaysHonorConstraints','bounds');
options = optimset('Display','notify','TolX',0.01*xmin);
SS.D = fminbnd(fun,xmin,xmax,options);

% Compute mass flux
SS.G = 4*SS.L*SS.mdot / (SS.A*SS.D);

% Compute Cf and St
SS.Re = SS.G*SS.D./SS.mu;
[SS.Cf, SS.St] = developed_flow(SS.Re,SS.Pr,'circular');

% Compute flow area
SS.Af = SS.mdot/SS.G;

% Check results
% Ntu1   = 4*SW.L/SW.D*mean(SW.St)
% Ntu2   = 4*SS.L/SS.D*mean(SS.St)
% ploss1 = mean(Ntu1*SW.G^2*SW.v.*(SW.Cf./SW.St)./(2*SW.p))
% ploss2 = mean(Ntu2*SS.G^2*SS.v.*(SS.Cf./SS.St)./(2*SS.p))
% keyboard

% Export results into HX structure
HX.shape = 'circular';
HX.L     = SW.L;       % Tube length, m
HX.AfT   = SW.Af + SS.Af; % Total flow area, m2
if mean(SW.p) > mean(SS.p)
    HX.D1    = SW.D;        % Tube diameter, m
    HX.AfR   = SS.Af/SW.Af; % Ratio of flow areas (shell/tube)
else
    HX.D1    = SS.D;        % Tube diameter, m
    HX.AfR   = SW.Af/SS.Af; % Ratio of flow areas (shell/tube)
end

end


%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = stream_losses(D,S)
% For a given value of hydraulic diameter (D), compute the addition of the
% non-dimensional thermal loss and non-dimensional pressure loss
% corresponding to one heat exchanger stream. The S structure must contain
% the following geometrical parameters: length (L), mass flow rate (mdot)
% and heat transfer area (A); and the following thermophysical properties:
% dynamic viscosity (mu), Prandtl number (Pr), specific volume (v) and
% pressure (p).

% Compute mass flux
S.G = 4*S.L*S.mdot / (S.A*D);

% Compute the Reynolds number, friction coefficient and Stanton number
S.Re = S.G*D./S.mu;
[S.Cf, S.St] = developed_flow(S.Re,S.Pr,'circular');

% Compute thermal loss factor (TL = 1/Ntu)
Ntu = 4*S.L/D*mean(S.St);
TL  = 1/Ntu;

% Compute pressure loss factor (PL = Dp/p)
PL = 2*S.G^2*mean(S.Cf.*S.v./S.p)*S.L/D;

% Add the thermal and the pressure loss factors
out = PL + TL;
%keyboard

% if any([TL > 1/Ntu_min, PL > ploss_max])
%     out = 1000*out;
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%