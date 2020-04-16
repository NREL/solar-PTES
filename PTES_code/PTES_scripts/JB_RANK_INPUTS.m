% Set atmospheric conditions and cycle parameters
T0      = 30 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
pmax    = 25e5;         % top pressure, Pa
PRch    = 3.0;          % charge pressure ratio
PRr     = 1.18;          % discharge pressure ratio: PRdis = PRch*PRr
PRr_min = 0.1;          % minimum PRr for optimisation
PRr_max = 3.0;          % maximum PRr for optimisation
setTmax = 1;            % set Tmax? (this option substitutes PRch)
Tmax    = 560 + 273.15; % maximum temp at compressor outlet, K

% Set Rankine-specific parameters
Ran_ptop    = 100e5;
Ran_pbotMIN = 0.05e5 ; % If condenser pressure decreases below this, the final stage chokes. Can't go to pressures below this, because who knows what happens.
Ran_Tbot0   = T0+15; %when discharging against the environment. This sets design condenser pressure.
Ran_TbotC   = 273.15+20; %when discharging against the cold stores

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

% Set parameters of Load structure
switch Load.mode
    case 0 % PTES
        fac = 100.0; % This can be used to more easily set the mass flow to obtain a desired power output
        stH = 8 ;
        ee  = 1;%0.6031 ;
        % This is the load scenario the plant is designed for
        Design_Load      = Load ;
        Design_Load.time = [stH/ee;stH].*3600;  % time spent in each load period, s
        Design_Load.type = ["chg";"dis"];    % type of load period
        Design_Load.mdot = [10*fac*ee;10*fac];  % working fluid mass flow rate, kg/s
        
        if Loffdesign
            % This is the actual load profile that the plant meets
            Load.time = [stH;stH].*3600;      % time spent in each load period, s
            Load.type = ["chg";"dis"];    % type of load period
            Load.mdot = [8.*fac;8.*fac];      % working fluid mass flow rate, kg/s
            T0_inc    = 0.0 ; % Increase in ambient temperature
        else
            Load = Design_Load ;
        end
                
    case 1 % Heat pump
        Design_Load      = Load ;
        stH              = 10 ;
        Design_Load.time = stH*3600;                  % time spent in each load period, s
        Design_Load.type = "chg";                     % type of load period
        Design_Load.mdot = 1000;                        % working fluid mass flow rate, kg/s
        
        if Loffdesign
            % This is the actual load profile that the plant meets
            Load.time = [stH].*3600;      % time spent in each load period, s
            Load.type = ["chg"];    % type of load period
            Load.mdot = [1000];      % working fluid mass flow rate, kg/s
            T0_inc    = 10.0 ; % Increase in ambient temperature
        else
            Load = Design_Load ;
        end
        
    case 2 % Heat engine (no cold tanks)
        
        Design_Load      = Load ;
        Design_Load.time = [0,10].*3600;                  % time spent in each load period, s
        Design_Load.type = ["sol";"dis"];                     % type of load period
        Design_Load.mdot = [0,10];                        % working fluid mass flow rate, kg/s
        
        if Loffdesign
            % This is the actual load profile that the plant meets
            Load.time = [0,10].*3600;      % time spent in each load period, s
            Load.type = ["sol";"dis"];    % type of load period
            Load.mdot = [0,10];      % working fluid mass flow rate, kg/s
            T0_inc    = 0;  % Increase in ambient temperature
        else
            Load = Design_Load ;
        end
        
    case 3 % JB charge, Rankine discharge
        fac = 1.0 ; % This can be used to more easily set the mass flow to obtain a desired power output 
        
        % This is the load scenario the plant is designed for
        Design_Load      = Load ;
        Design_Load.time = [10;4;10;10].*3600;          % time spent in each load period, s
        Design_Load.type = ["chg";"str";"ran";"ran"];   % type of load period
        Design_Load.mdot = [10*fac;0;1*fac;1*fac];      % working fluid mass flow rate, kg/s
        Design_Load.options.useCold = [0,0,0,0];        % Use cold stores during Rankine discharge? This should be set to 0 for design cases of retrofits.
        
        if Loffdesign
            % This is the actual load profile that the plant meets
            Load.time = [10;4;10;10].*3600;         % time spent in each load period, s
            Load.type = ["chg";"str";"ran";"ran"];  % type of load period
            Load.mdot = [10*fac;0;1*fac;1*fac];     % working fluid mass flow rate, kg/s
            Load.options.useCold = [0,0,1,0];        % Use cold stores during Rankine discharge?
            T0_inc    = 0.0 ; % Increase in ambient temperature
        else
            Load = Design_Load ;
        end
end
Design_Load.num  = numel(Design_Load.time);
Design_Load.ind  = 1:Design_Load.num;
Load.num  = numel(Load.time);
Load.ind  = 1:Load.num;

% Hot storage tanks
switch PBmode
    case 0
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
    case 1
        % Hot storage - estimated temperatures
        TH_dis0 = T0;   % initial temperature of discharged hot fluid, K
        TH_chg0 = Tmax; % initial temperature of charged hot fluid, K
        % Cold storage - estimated temperatures
        TC_dis0 = T0;           % initial temperature of discharged cold fluid, K
        TC_chg0 = -150 + 273.15;        % initial temperature of charged cold fluid, K
        
        % Also go through and rename chg to chgPB so that correct routines
        % are called from PTES.m
        for ii = 1 : Design_Load.num
           if strcmp(Design_Load.type(ii),"chg")
               Design_Load.type(ii) = "chgPB" ;
           elseif strcmp(Design_Load.type(ii),"dis")
               Design_Load.type(ii) = "disPB" ;
           end
        end
        for ii = 1 : Load.num
           if strcmp(Load.type(ii),"chg")
               Load.type(ii) = "chgPB" ;
           elseif strcmp(Load.type(ii),"dis")
               Load.type(ii) = "disPB" ;
           end
        end
    case 2
        error('Not implemented')
end

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
% dat.cp   = RP1('PT_INPUTS',dat.P0,dat.T0,'CPMASS',gas_temp) ;
% dat.cv   = RP1('PT_INPUTS',dat.P0,dat.T0,'CVMASS',gas_temp) ;
% dat.mu0  = 17.81e-6 ;
% dat.TVref = 300.55 ;
% dat.S     = 111 ;
% dat.k     = 0.026 ;
% 
% gas      = fluid_class('Nitrogen','WF','IDL',dat,Load.num,30);

if Load.mode==3
    % 'TTSE' interpolation is NOT recommended for steam when reading values
    % close to the saturation curve. Use 'HEOS' or 'BICUBIC&HEOS'
    steam = fluid_class('Water','WF','CP','BICUBIC&HEOS',Load.num,30);
end

% Save copy of input file in "Outputs" folder
copyfile(['./PTES_scripts/',mfilename,'.m'],'./Outputs/')
