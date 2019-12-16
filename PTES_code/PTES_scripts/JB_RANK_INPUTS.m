% Set atmospheric conditions and cycle parameters
T0      = 30 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
pmax    = 25e5;         % top pressure, Pa
PRch    = 6.0;          % charge pressure ratio
PRr     = 1.5;          % discharge pressure ratio: PRdis = PRch*PRr
PRr_min = 0.1;          % minimum PRr for optimisation
PRr_max = 3.0;          % maximum PRr for optimisation
setTmax = 1;            % set Tmax? (this option substitutes PRch)
Tmax    = 570 + 273.15; % maximum temp at compressor outlet, K

% Set Rankine-specific parameters
Ran_ptop  = 100e5;
Ran_Tbot0 = T0+15; %when discharging against the environment
Ran_TbotC = 273.15+20; %when discharging against the cold stores

% Set component parameters
eta   = 0.90;  % polytropic efficiency
eff   = 0.97;  % heat exchanger effectiveness
ploss = 0.01;  % pressure loss in HEXs

% Number of intercooled/interheated compressions/expansions
Nc_ch = 1; % number of compressions during charge
Ne_ch = 2; % number of expansions during charge
nH    = max([2,Nc_ch]); % number of hot fluid streams
nC    = Ne_ch;          % number of cold fluid streams

% Number of hot and cold stores IN SERIES
Ncld = 1; % number of cold stores. Not implemented for >2
Nhot = 1; % number of hot stores. Not implemented for >2

% Set parameters of Load structure
switch Load.mode
    case 0 % PTES
        Load.time = [10;10;4;6;10].*3600;               % time spent in each load period, s
        Load.type = ["chg";"chg";"str";"str";"dis"];    % type of load period
        Load.mdot = [6;4;0;0;10];                       % working fluid mass flow rate, kg/s
        
    case 1 % Heat pump
        Load.time = 10.*3600;                  % time spent in each load period, s
        Load.type = "chg";                     % type of load period
        Load.mdot = 10;                        % working fluid mass flow rate, kg/s
        
    case 2 % Heat engine (no cold tanks)
        Load.time = [0,10].*3600;                  % time spent in each load period, s
        Load.type = ["sol","dis"];                 % type of load period
        Load.mdot = [0,10];                        % working fluid mass flow rate, kg/s
        
    case 3 % JB charge, Rankine discharge
        Load.time = [10;4;10;10].*3600;        % time spent in each load period, s
        Load.type = ["chg";"str";"ran";"ran"];    % type of load period
        Load.mdot = [10;0;1;1];              % working fluid mass flow rate, kg/s
        Load.options.useCold = [0,0,1,0]; %Use cold stores during Rankine discharge?
%         Load.time = [10;4;15].*3600;        % time spent in each load period, s
%         Load.type = ["chg";"str";"ran"];    % type of load period
%         Load.mdot = [10;0;1];              % working fluid mass flow rate, kg/s
%         Load.options.useCold = [0,0,0]; %Use cold stores during Rankine discharge?
end
Load.num  = numel(Load.time);
Load.ind  = 1:Load.num;

% Hot storage tanks
fHname  = 'SolarSalt';  % fluid name
TH_dis0 = 230 + 273.15; % initial temperature of discharged hot fluid, K
MH_dis0 = 1e6;          % initial mass of discharged hot fluid, kg
TH_chg0 = 550 + 273.15; % initial temperature of charged hot fluid, K
MH_chg0 = 0.00*MH_dis0; % initial mass of charged hot fluid, kg

% Cold storage tanks
fCname  = 'INCOMP::MEG2[0.56]'; % fluid name
TC_dis0 = T0;           % initial temperature of discharged cold fluid, K
MC_dis0 = 1e6;          % initial mass of discharged cold fluid, kg
TC_chg0 = T0-50;        % initial temperature of charged cold fluid, K
MC_chg0 = 0.00*MC_dis0; % initial mass of charged cold fluid, kg

% Set fluids. 'WF' or 'SF' indicates working fluid or storage fluid. 'CP'
% or 'TAB' indicate CoolProp or Tabular reading modes. 'backend' is used by
% CoolProp (either 'HEOS', 'TTSE' or 'BICUBIC&HEOS'). 'HEOS' is the most
% accurate method but is very slow. 'BICUBIC&HEOS' is recommended over
% 'TTSE' for speed and accuracy. 'num' indicates number of preallocated
% elements in state arrays.
gas = fluid_class('Nitrogen','WF','CP','BICUBIC&HEOS',Load.num,30);
if Load.mode==3
    % 'TTSE' interpolation is NOT recommended for steam when reading values
    % close to the saturation curve. Use 'HEOS' or 'BICUBIC&HEOS'
    steam = fluid_class('Water','WF','CP','HEOS',Load.num,30);
end

% Set heat exchangers (this is temporary, it will be done properly using a
% class constructor
HX_CONDEN.model = 'eff';
HX_CONDEN.eff = eff;
HX_CONDEN.ploss = 0;
HX_CONDEN.stage_type = 'hex';
HX_CONDEN.NX = 100;
HX_REHEAT.model = 'eff';
HX_REHEAT.eff = eff;
HX_REHEAT.ploss = ploss;
HX_REHEAT.stage_type = 'hex';
HX_REHEAT.NX = 100;
HX_BOILER.model = 'eff';
HX_BOILER.eff = eff;
HX_BOILER.ploss = ploss;
HX_BOILER.stage_type = 'hex';
HX_BOILER.NX = 100;
HX_ACC.model = 'eff';
HX_ACC.eff = eff;
HX_ACC.ploss = ploss;
HX_ACC.stage_type = 'regen';
HX_ACC.NX = 100;

% Save copy of input file in "Outputs" folder
copyfile(['./PTES_scripts/',mfilename,'.m'],'./Outputs/')