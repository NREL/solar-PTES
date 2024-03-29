% Set atmospheric conditions and cycle parameters
T0      = 40 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
pmax    = 25e5;         % top pressure, Pa
PRch    = 2.5;          % charge pressure ratio
PRr     = 1.18;          % discharge pressure ratio: PRdis = PRch*PRr
PRr_min = 0.1;          % minimum PRr for optimisation
PRr_max = 3.0;          % maximum PRr for optimisation
setTmax = 1;            % set Tmax? (this option substitutes PRch)
Tmax    = 560 + 273.15; % maximum temp at compressor outlet, K

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
Ne_ch = 1; % number of expansions during charge
nH    = max([2,Nc_ch]); % number of hot fluid streams
nC    = Ne_ch;          % number of cold fluid streams

% Number of hot and cold stores IN SERIES
Ncld = 1; % number of cold stores. Not implemented for >2
Nhot = 1; % number of hot stores. Not implemented for >2

% Design or off-design analysis?
Loffdesign = 0; % L for Logical. 0 just run design case. 1 run design case then off-design load cycle.

% Set parameters of Load structure
switch Load.mode
    case 0 % PTES
        fac = 1.0; % This can be used to more easily set the mass flow to obtain a desired power output
        stH = 8 ;
        
        % This is the load scenario the plant is designed for
        Design_Load      = Load ;
        Design_Load.time = [stH;stH].*3600;  % time spent in each load period, s
        Design_Load.type = ["chg";"dis"];    % type of load period
        Design_Load.mdot = [10*fac;10*fac];  % working fluid mass flow rate, kg/s
        
        if Loffdesign
            % This is the actual load profile that the plant meets
            Load.time = [stH;4;stH].*3600;      % time spent in each load period, s
            Load.type = ["chg";"str";"dis"];    % type of load period
            Load.mdot = [10*fac;0;10*fac];      % working fluid mass flow rate, kg/s
        else
            Load = Design_Load ;
        end
                
    case 1 % Heat pump
        Load.time = 10.*3600;                  % time spent in each load period, s
        Load.type = "chg";                     % type of load period
        Load.mdot = 10;                        % working fluid mass flow rate, kg/s
        
    case 2 % Heat engine (no cold tanks)
        Load.time = [0,10].*3600;                  % time spent in each load period, s
        Load.type = ["sol","dis"];                 % type of load period
        Load.mdot = [0,10];                        % working fluid mass flow rate, kg/s
        
    case 3 % JB charge, Rankine discharge
         fac = 1.0 ; % THis can be used to more easily set the mass flow to obtain a desired power output 
         Load.time = [10;4;10;10].*3600;        % time spent in each load period, s
         Load.type = ["chg";"str";"ran";"ran"];    % type of load period
         Load.mdot = [10*fac;0;1*fac;1*fac];              % working fluid mass flow rate, kg/s
         Load.options.useCold = [0,0,1,0]; %Use cold stores during Rankine discharge?
end
Load.num  = numel(Load.time);
Load.ind  = 1:Load.num;

% Hot storage tanks
fHname  = 'SolarSalt';  % fluid name
TH_dis0 = 350 + 273.15; % initial temperature of discharged hot fluid, K
MH_dis0 = 1e9;          % initial mass of discharged hot fluid, kg
TH_chg0 = 570 + 273.15; % initial temperature of charged hot fluid, K
MH_chg0 = 0.00*MH_dis0; % initial mass of charged hot fluid, kg

% Cold storage tanks
fCname  = 'INCOMP::MEG2[0.56]'; % fluid name
TC_dis0 = T0;           % initial temperature of discharged cold fluid, K
MC_dis0 = 1e9;          % initial mass of discharged cold fluid, kg
TC_chg0 = T0-50;        % initial temperature of charged cold fluid, K
MC_chg0 = 0.00*MC_dis0; % initial mass of charged cold fluid, kg

% Set fluids. 'WF' or 'SF' indicates working fluid or storage fluid. 'CP'
% or 'TAB' indicate CoolProp or Tabular reading modes. 'backend' is used by
% CoolProp (either 'HEOS', 'TTSE' or 'BICUBIC&HEOS'). 'HEOS' is the most
% accurate method but is very slow. 'BICUBIC&HEOS' is recommended over
% 'TTSE' for speed and accuracy. 'num' indicates number of preallocated
% elements in state arrays.
gas = fluid_class('Nitrogen','WF','CP','BICUBIC&HEOS',Load.num,30);

% Set up an ideal gas - should run faster
% gas_temp = fluid_class('Nitrogen','WF','CP','BICUBIC&HEOS',Load.num,30);
% dat.T0   = 300 ;
% dat.P0   = 1e5 ;
% dat.cp   = CP1('PT_INPUTS',dat.P0,dat.T0,'CPMASS',gas_temp.handle) ;
% dat.cv   = CP1('PT_INPUTS',dat.P0,dat.T0,'CVMASS',gas_temp.handle) ;
% dat.mu0  = 17.81e-6 ;
% dat.TVref = 300.55 ;
% dat.S     = 111 ;
% dat.k     = 0.026 ;
% 
% gas      = fluid_class('Nitrogen','WF','IDL',dat,Load.num,30);

if Load.mode==3
    % 'TTSE' interpolation is NOT recommended for steam when reading values
    % close to the saturation curve. Use 'HEOS' or 'BICUBIC&HEOS'
    steam = fluid_class('Water','WF','CP','HEOS',Load.num,30);
end

% Make heat exchangers
switch Load.mode
    case {0,1,2}
        % Call HX classes for ideal-gas PTES cycle
        HX(1) = hx_class('hot',  'hex',   'eff', eff, ploss,  11, 100, Load.num, Load.num) ; % Hot heat exchanger
        HX(2) = hx_class('cold', 'hex',   'eff', eff, ploss,  11, 100, Load.num, Load.num) ; % Cold heat exchanger
        HX(3) = hx_class('regen','regen', 'eff', eff, ploss,  11, 100, Load.num, Load.num) ; % Recuperator
        HX(4) = hx_class('rej',  'hex',   'eff', eff, ploss, 2, 100, Load.num, Load.num) ; % Heat rejection unit
    case 3
        % Call HX classes for ideal-gas PTES heat pump with Rankine cycle discharge
        HX(1) = hx_class('hot',  'hex',   'eff', eff, ploss,  11, 100, Load.num, Load.num) ; % Hot heat exchanger
        HX(2) = hx_class('cold', 'hex',   'eff', eff, ploss,  11, 100, Load.num, Load.num) ; % Cold heat exchanger
        HX(3) = hx_class('regen','regen', 'eff', eff, ploss,  1, 100, Load.num, Load.num) ; % Recuperator
        HX(4) = hx_class('rej',  'hex',   'eff', eff, ploss, 0, 100, Load.num, Load.num) ; % Heat rejection unit
        
        HX(5) = hx_class('cold', 'hex',   'eff', eff, 0.1/100, 0, 100, Load.num, Load.num) ; % Condenser
        HX(6) = hx_class('hot',  'hex',   'eff', eff, ploss,  0, 100, Load.num, Load.num) ; % Reheat
        HX(7) = hx_class('hot',  'hex',   'eff', eff, ploss,  0, 100, Load.num, Load.num) ; % Boiler
        HX(8) = hx_class('rej',  'regen', 'eff', eff, 0.1/100,0, 100, Load.num, Load.num) ; % Air-cooled condenser
end

% Options for specifying heat exchanger geometry
% This will probably be expanded over time
for i = 1 : length(HX)
   HX(i).LestA = true ;
   HX(i).D1    = 0.025 ;
end


% Save copy of input file in "Outputs" folder
copyfile(['./PTES_scripts/',mfilename,'.m'],'./Outputs/')