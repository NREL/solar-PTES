% Set atmospheric conditions and cycle parameters
T0      = 30 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
Tmax    = 570 + 273.15; % maximum temperature, K
pmax_LA = 200e5;         % top pressure of LA cycle, Pa
T0_inc  = 5.0 ; % Increment above ambient temperature that gas is cooled to

% Set compressor/expander parameters
eta   = 0.90;  % polytropic efficiency

% Number of intercooled/interheated compressions/expansions
Nc_ch = 2; % number of compressions during charge
Ne_ch = 1; % number of expansions during charge
nH    = 2; % number of hot fluid streams
nC    = 1; % number of cold fluid streams

% Number of hot and cold stores IN SERIES
Ncld = 1; % number of cold stores. Not implemented for >2
Nhot = 1; % number of hot stores. Not implemented for >2

% Set parameters of Load structure
% This is the load scenario the plant is designed for
Design_Load      = Load ;
Design_Load.time = [10;10].*3600;    % time spent in each load period, s
Design_Load.type = ["chgCC";"disCC"];    % type of load period
Design_Load.mdot = [10;10];  % working fluid mass flow rate, kg/s
Load = Design_Load ;
Design_Load.num  = numel(Design_Load.time);
Design_Load.ind  = 1:Design_Load.num;
Load.num  = numel(Load.time);
Load.ind  = 1:Load.num;

% Hot storage tanks
fHname  = "SolarSalt";  % fluid name
TH_dis0 = 300 + 273.15; % initial temperature of discharged hot fluid, K
MH_dis0 = 1e9;          % initial mass of discharged hot fluid, kg
TH_chg0 = 570 + 273.15; % initial temperature of charged hot fluid, K
MH_chg0 = 0.00*MH_dis0; % initial mass of charged hot fluid, kg

% Medium storage tanks
fMname  = "MineralOil"; % fluid name
TM_dis0 = T0;           % initial temperature of discharged medium fluid, K
MM_dis0 = 1e9;          % initial mass of discharged medium fluid, kg
TM_chg0 = 300 + 273.15; % initial temperature of charged medium fluid, K
MM_chg0 = 0.00*MM_dis0; % initial mass of charged medium fluid, kg

% Cold storage tanks
fCname  = 'Isopentane'; % fluid name
TC_dis0 = T0-150;       % initial temperature of discharged cold fluid, K
MC_dis0 = 1e9;          % initial mass of discharged cold fluid, kg
TC_chg0 = T0;           % initial temperature of charged cold fluid, K
MC_chg0 = 0.00*MC_dis0; % initial mass of charged cold fluid, kg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COST MODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CCMPmode = [12,13,15,15] ; % Charging compressor cost mode
CEXPmode = [41,42,43,44] ; % Charging expander cost mode
DCMPmode = [12,13,15,15] ; % Discharging compressor cost mode
DEXPmode = [41,42,43,44] ; % Discharging expander cost mode

PMPmode = [60,61,62]; % Pump cost mode
FANmode = [71,70,71,72]; % Fan cost mode

hotHXmode = [1,2,3,4,5,10,11]; % Heat exchanger - hot cost mode
cldHXmode = [1,2,3,4,5,10,11]; % Heat exchanger - cold cost mode
rcpHXmode = [1,2,3,4,5,10,11]; % Heat exchanger - recuperator cost mode
rejHXmode = [30,31,32,33]; % Heat exchanger - rejection cost mode

GENmode = [2,1,3] ; % Motor-generator

HTmode.tankmode  = [5,1,2,3,4,6] ; % Cost mode for hot tank container cost
HTmode.fld_cost  = [0.8,0.5,1.3] ; % Hot tank fluid cost, $/kg
HTmode.ins_cost  = [30,5,50] ; % Insulation material, %/kg
HTmode.ins_k     = 0.08 ; % Thermal conductivity of insulation
HTmode.ins_rho   = 150 ; % Density of insulation
HTmode.tau       = 200 ; % Number of days before all heat leaks out of tank
HTmode.AR        = 1.0 ; % Aspect ratio (L/D) of tank
HTmode.over_fac  = 1.1 ; % How much larger is inner tank volume than the fluid volume

CTmode.tankmode  = [5,1,2,3,4,6] ; % Cost mode for cold tank container cost
CTmode.fld_cost  = [0.56,0.1,1] ; % Cold tank fluid cost, $/kg
CTmode.ins_cost  = [30,5,50] ; % Insulation material, %/kg
CTmode.ins_k     = 0.08 ; % Thermal conductivity of insulation
CTmode.ins_rho   = 150 ; % Density of insulation
CTmode.tau       = 200 ; % Number of days before all heat leaks out of tank
CTmode.AR        = 1.0 ; % Aspect ratio (L/D) of tank
CTmode.over_fac  = 1.1 ; % How much larger is inner tank volume than the fluid volume

ATmode.tankmode  = 0 ;
ATmode.fld_cost  = 0 ; % Cold tank fluid cost, $/kg
ATmode.ins_cost  = 0 ; % Insulation material, %/kg

% Set fluids. 'WF' or 'SF' indicates working fluid or storage fluid. 'CP'
% or 'TAB' indicate CoolProp or Tabular reading modes. 'backend' is used by
% CoolProp (either 'HEOS', 'TTSE' or 'BICUBIC&HEOS'). 'HEOS' is the most
% accurate method but is very slow. 'BICUBIC&HEOS' is recommended over
% 'TTSE' for speed and accuracy. 'num' indicates number of preallocated
% elements in state arrays.
gas1 = fluid_class('Nitrogen','WF','CP','BICUBIC&HEOS',Load.num,30); %LAES side (air)
gas2 = fluid_class('Neon','WF','CP','BICUBIC&HEOS',Load.num,30);     %PTES side

% Save copy of input file in "Outputs" folder
copyfile(['./PTES_scripts/',mfilename,'.m'],'./Outputs/')
