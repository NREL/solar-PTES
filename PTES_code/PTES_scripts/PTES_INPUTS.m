%%% SET INPUT VARIABLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set atmospheric conditions and cycle parameters
T0      = 27 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
TC_0    = T0;           % temperature of discharged cold fluid, K
TH_0    = 250 + 273.15; % temperature of discharged hot fluid, K
mdot    = 10;           % working fluid mass flow rate, kg/s
t_ch    = 10 * 3600;    % charge time, s
pmax    = 100e5;        % top pressure, Pa
PRch    = 3.5;          % charge pressure ratio
PRr     = 1.4;          % discharge pressure ratio: PRdis = PRch*PRr
PRr_min = 0.1;          % minimum PRr for optimisation
PRr_max = 3.0;          % maximum PRr for optimisation
setTmax = 1;            % set Tmax? (this option substitutes PRch)
Tmax    = 570 + 273.15; % maximum temp at compressor outlet, K

% Set operation modes
mode       = 0; % cycle mode: 0=PTES, 1=Heat pump, 2=Heat engine
Nc_ch      = 1; % number of compressions during charge
Ne_ch      = 2; % number of expansions during charge

% Set working fluids, storage fluids, and heat rejection streams. 'WF' or
% 'SF' indicates working fluid or storage fluid. 'CP' or 'TAB' indicate
% CoolProp or Tabular reading modes. 'backend' is used by CoolProp (either
% 'HEOS', 'TTSE' or 'BICUBIC&HEOS'). 'num' indicates number of preallocated
% elements in state arrays.
% Working fluid
gas = fluid_class('Nitrogen','WF','CP','TTSE',30);

% Set double tanks. Double tanks have four states:
% 1=begin charge, 2=end charge, 3=start discharge, 4=end discharge
HT = double_tank_class(4);  %hot double tank
CT = double_tank_class(4);  %cold double tank

% Storage fluids
fluidH(1:Nc_ch) = fluid_class('SolarSalt','SF','TAB',0,2);
fluidC(1:Ne_ch) = fluid_class('INCOMP::MEG2[0.56]','SF','TAB',0,2);
% Heat rejection streams
environ = environment_class(T0,p0,10);

% Set component parameters
eta    = 0.90;  % polytropic efficiency
eff    = 0.97;  % heat exchanger effectiveness
ploss  = 0.01;  % pressure loss in HEXs

multi_run  = 0; % run cycle several times with different parameters?
optimise   = 1; % optimise cycle?
make_plots = 1; % make plots?
save_figs  = 0; % save figures at the end?

% Set number of points to plot each stage
num = 10;

% Variables to run cycle multiple times and plot curves. The variables must
% have been defined in the PTES_SET_MULTI_RUN script
if multi_run==1
    Vpnt = 'TH_0';   % variable along curve
    Npnt = 5;       % points on curve
    pnt1 = 300+273.15;    % min value
    pnt2 = 400+273.15;    % max value
    Apnt = linspace(pnt1,pnt2,Npnt); % array
    Vcrv = 'Ne_ch';  % variable between curves
    Acrv = 1;%[1,2,3];
    Ncrv = numel(Acrv);
    %     Vcrv = 'eta';  % variable between curves
    %     Ncrv = 3;      % number of curves
    %     crv1 = 0.95;   % min value
    %     crv2 = 0.99;   % max value
    %     Acrv = linspace(crv1,crv2,Ncrv); % array
else
    Npnt=1; Ncrv=1;
    % Start new logfile
    delete ./Outputs/log.txt
    diary  ./Outputs/log.txt
end

% Save copy of input file in "Outputs" folder
copyfile(['./PTES_scripts/',mfilename,'.m'],'./Outputs/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%