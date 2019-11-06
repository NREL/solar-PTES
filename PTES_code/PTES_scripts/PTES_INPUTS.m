%%% SET INPUT VARIABLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set atmospheric conditions and cycle parameters
T0      = 30 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
pmax    = 80e5;        % top pressure, Pa
PRch    = 3.5;          % charge pressure ratio
PRr     = 1.05;          % discharge pressure ratio: PRdis = PRch*PRr
PRr_min = 0.1;          % minimum PRr for optimisation
PRr_max = 3.0;          % maximum PRr for optimisation
setTmax = 1;            % set Tmax? (this option substitutes PRch)
Tmax    = 570 + 273.15; % maximum temp at compressor outlet, K
% Number of intercooled/interheated compressions/expansions
Nc_ch = 1; % number of compressions during charge
Ne_ch = 2; % number of expansions during charge

% Hot storage tanks
fHname  = 'SolarSalt';  % fluid name
TH_dis0 = 230 + 273.15; % initial temperature of discharged hot fluid, K
MH_dis0 = 1e6;          % initial mass of discharged hot fluid, kg
TH_chg0 = 550 + 273.15; % initial temperature of charged hot fluid, K
MH_chg0 = 0.00*MH_dis0; % initial mass of charged hot fluid, kg
nH      = Nc_ch;        % number of hot fluid streams
Nhot = 1;               % number of hot tanks
% Cold storage tanks
fCname  = 'INCOMP::MEG2[0.56]'; % fluid name
TC_dis0 = T0;           % initial temperature of discharged cold fluid, K
MC_dis0 = 1e6;          % initial mass of discharged cold fluid, kg
TC_chg0 = T0-50;        % initial temperature of charged cold fluid, K
MC_chg0 = 0.00*MC_dis0; % initial mass of charged cold fluid, kg
nC      = Ne_ch;        % number of cold fluid streams
Ncld = 1;               % number of cold tanks

% The Load structure stores information about the duration, type of cycle
% (charge, storage or discharge) and mass flow rate of each time period.
Load.mode = 3;
switch Load.mode
    case 0 % PTES
        Load.time = [10;10;4;10].*3600;        % time spent in each load period, s
        Load.type = ["chg";"chg";"str";"dis"];    % type of load period
        Load.mdot = [6;4;0;10];              % working fluid mass flow rate, kg/s
    case 1 % Heat pump
        Load.time = 10.*3600;                  % time spent in each load period, s
        Load.type = "chg";                     % type of load period
        Load.mdot = 10;                        % working fluid mass flow rate, kg/s
    case 2 % Heat engine (no cold tanks)
        error('not implemented')
        Load.time = [0,10].*3600;                  % time spent in each load period, s
        Load.type = ["sol","dis"];                 % type of load period
        Load.mdot = [0,10];                        % working fluid mass flow rate, kg/s
    case 3 % JB charge, Rankine discharge
        Load.time = [10;4;10].*3600;        % time spent in each load period, s
        Load.type = ["chg";"str";"ran"];    % type of load period
        Load.mdot = [10;0;10];              % working fluid mass flow rate, kg/s
        nH = 2;
end
Load.num = numel(Load.time);
Load.ind = 1:Load.num;

% Set working fluids, storage fluids, and heat rejection streams. 'WF' or
% 'SF' indicates working fluid or storage fluid. 'CP' or 'TAB' indicate
% CoolProp or Tabular reading modes. 'backend' is used by CoolProp (either
% 'HEOS', 'TTSE' or 'BICUBIC&HEOS'). 'num' indicates number of preallocated
% elements in state arrays.
% Working fluid
gas = fluid_class('Nitrogen','WF','CP','TTSE',Load.num,30);
if Load.mode ==3
    steam = fluid_class('Water','WF','CP','TTSE',Load.num,30);
end
% Storage fluids
fluidH(1:nH) = fluid_class(fHname,'SF','TAB',NaN,Load.num,10);
fluidC(1:nC) = fluid_class(fCname,'SF','TAB',NaN,Load.num,10);
% Set double tanks
HT = double_tank_class(fluidH,TH_dis0,p0,MH_dis0,TH_chg0,p0,MH_chg0,T0,Load.num+1); %hot double tank
CT = double_tank_class(fluidC,TC_dis0,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0,Load.num+1); %cold double tank
% Heat rejection streams
environ = environment_class(T0,p0,10);

% Set component parameters
eta   = 0.90;  % polytropic efficiency
eff   = 0.97;  % heat exchanger effectiveness
ploss = 0.01;  % pressure loss in HEXs

multi_run  = 0; % run cycle several times with different parameters?
optimise   = 0; % optimise cycle?
make_plots = 1; % make plots?
save_figs  = 0; % save figures at the end?

% Set number of points to plot each stage
num = 25;

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