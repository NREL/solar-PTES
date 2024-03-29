%%% SET INPUT VARIABLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set atmospheric conditions and cycle parameters
T0      = 27 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
pmax    = 250e5;        % top pressure, Pa
PRch    = 3.0;          % charge pressure ratio
PRr     = 0.80;          % discharge pressure ratio: PRdis = PRch*PRr
PRr_min = 0.1;          % minimum PRr for optimisation
PRr_max = 3.0;          % maximum PRr for optimisation
setTmax = 0;            % set Tmax? (this option substitutes PRch)
Tmax    = 500 + 273.15; % maximum temp at compressor outlet, K
Lcld    = false ;       % Make cold store as cold as possible?

% Number of intercooled/interheated compressions/expansions
Nc_ch = 1; % number of compressions during charge
Ne_ch = 1; % number of expansions during charge

% Number of hot and cold stores IN SERIES
Ncld = 2; % number of cold stores. Not implemented for >1
Nhot = 1; % number of hot stores. Not implemented for >2

% Number of recuperators
Nrcp = 0 ; % Can be 0,1,2. If 0 may need two hot stores. If 2 may require a recompression. 
switch Nrcp
    case 0
        % Hot storage tanks
        fHname  = 'MineralOil'; % fluid name
        TH_dis0 = T0;           % initial temperature of discharged hot fluid, K
        MH_dis0 = 1e6;          % initial mass of discharged hot fluid, kg
        TH_chg0 = 250 + 273.15; % initial temperature of charged hot fluid, K
        MH_chg0 = 0.00*MH_dis0; % initial mass of charged hot fluid, kg
        TH_int  = 100 + 273.15 ;% Intermediate temperature between two hot stores
        % Cold storage tanks
        fCname  = 'INCOMP::MEG2[0.56]'; % fluid name
        TC_dis0 = 100 + 273.15; % initial temperature of discharged cold fluid, K
        MC_dis0 = 1e6;          % initial mass of discharged cold fluid, kg
        TC_chg0 = T0-5;         % initial temperature of charged cold fluid, K
        MC_chg0 = 0.00*MC_dis0; % initial mass of charged cold fluid, kg
        TC_int  = 50 + 273.15 ; % Intermediate temperature between two cold stores
    case 1
        % Hot storage tanks
        fHname  = 'SolarSalt';  % fluid name
        TH_dis0 = T0 + 273.15;  % initial temperature of discharged hot fluid, K
        MH_dis0 = 1e6;          % initial mass of discharged hot fluid, kg
        TH_chg0 = 550 + 273.15; % initial temperature of charged hot fluid, K
        MH_chg0 = 0.00*MH_dis0; % initial mass of charged hot fluid, kg
        % Cold storage tanks
        fCname  = 'INCOMP::MEG2[0.56]'; % fluid name
        TC_dis0 = T0 + 0;           % initial temperature of discharged cold fluid, K
        MC_dis0 = 1e6;          % initial mass of discharged cold fluid, kg
        TC_chg0 = T0-5;        % initial temperature of charged cold fluid, K
        MC_chg0 = 0.00*MC_dis0; % initial mass of charged cold fluid, kg
    case 2
end

% The Load structure stores information about the duration, type of cycle
% (charge, storage or discharge) and mass flow rate of each time period.
Load.mode = 0;
switch Load.mode
    case 0 % PTES
        Load.time = [10;4;10].*3600;          % time spent in each load period, s
        Load.type = ["chg";"str";"dis"]; % type of load period
        Load.mdot = [10;0;10];              % working fluid mass flow rate, kg/s
        Load.num  = numel(Load.time);
    case 1 % Heat pump
        Load.time = 10.*3600;                  % time spent in each load period, s
        Load.type = "chg";                     % type of load period
        Load.mdot = 10;                        % working fluid mass flow rate, kg/s
        Load.num  = numel(Load.time);
    case 2 % Heat engine (no cold tanks)
        error('not implemented')
        Load.time = [0,10].*3600;                  % time spent in each load period, s
        Load.type = ["sol","dis"];                 % type of load period
        Load.mdot = [0,10];                        % working fluid mass flow rate, kg/s
        Load.num  = numel(Load.time);
end

% Set working fluids, storage fluids, and heat rejection streams. 'WF' or
% 'SF' indicates working fluid or storage fluid. 'CP' or 'TAB' indicate
% CoolProp or Tabular reading modes. 'backend' is used by CoolProp (either
% 'HEOS', 'TTSE' or 'BICUBIC&HEOS'). 'num' indicates number of preallocated
% elements in state arrays.
% Working fluid
gas = fluid_class('CarbonDioxide','WF','CP','TTSE',Load.num,30);

% Set double tanks
switch Ncld
    case 1
        fluidC(1:Ne_ch) = fluid_class(fCname,'SF','TAB',NaN,Load.num,2); % Storage fluid
        CT  = double_tank_class(fluidC,TC_dis0,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0,Load.num+1); %cold double tank
    case 2
        fluidC(1:Ne_ch)  = fluid_class(fCname,'SF','TAB',NaN,Load.num,2);
        fluidC2(1:Ne_ch) = fluid_class(fCname,'SF','TAB',NaN,Load.num,2);
        CT  = double_tank_class(fluidC,TC_int,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0,Load.num+1); %cold double tank
        CT2 = double_tank_class(fluidC2,TC_dis0,p0,MC_dis0,TC_int,p0,MC_chg0,T0,Load.num+1); %cold double tank
    case 3
        error('Not implemented')
end
% Hot tanks
switch Nhot
    case 1
        fluidH(1:Nc_ch) = fluid_class(fHname,'SF','TAB',NaN,Load.num,2*Nhot);
        HT  = double_tank_class(fluidH,TH_dis0,p0,MH_dis0,TH_chg0,p0,MH_chg0,T0,Load.num+1); %hot double tank
    case 2
        fluidH(1:Nc_ch)  = fluid_class(fHname,'SF','TAB',NaN,Load.num,2);
        fluidH2(1:Nc_ch) = fluid_class(fHname,'SF','TAB',NaN,Load.num,2);
        HT  = double_tank_class(fluidH,TH_int,p0,MH_dis0,TH_chg0,p0,MH_chg0,T0,Load.num+1); %hot double tank
        HT2 = double_tank_class(fluidH2,TH_dis0,p0,MH_dis0,TH_int,p0,MH_chg0,T0,Load.num+1); %hot double tank
    case 3
        error('Not implemented')
end


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
%copyfile(['./PTES_scripts/',mfilename,'.m'],'./Outputs/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%