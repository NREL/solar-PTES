% Set atmospheric conditions and cycle parameters
T0      = 40 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
pmax    = 25e5;         % top pressure, Pa
PRch    = 1.5;          % charge pressure ratio
PRr     = 1.15;          % discharge pressure ratio: PRdis = PRch*PRr
PRr_min = 0.1;          % minimum PRr for optimisation
PRr_max = 3.0;          % maximum PRr for optimisation
LPRr    = 0 ;           % Logical. Estimate optimal PRr after charging run.
setTmax = 1;            % set Tmax? (this option substitutes PRch)
Tmax    = 570 + 273.15; % maximum temp at compressor outlet, K

% Set Rankine-specific parameters
Ran_ptop    = 100e5;
Ran_pbotMIN = 0.02e5 ; % If condenser pressure decreases below this, the final stage chokes. Can't go to pressures below this, because who knows what happens.
Ran_Tbot0   = T0 + 5; %when discharging against the environment. This sets design condenser pressure.
Ran_TbotC   = 273.15+20; %when discharging against the cold stores

% Set compressor/expander parameters
eta   = 0.90;  % polytropic efficiency

% Number of intercooled/interheated compressions/expansions
Nc_ch = 1; % number of compressions during charge
Ne_ch = 3; % number of expansions during charge
nH    = max([2,Nc_ch]); % number of hot fluid streams
nC    = Ne_ch;          % number of cold fluid streams

% Number of hot and cold stores IN SERIES
Ncld = 1; % number of cold stores. Not implemented for >2
Nhot = 1; % number of hot stores. Not implemented for >2

% Set parameters of Load structure
switch Load.mode
    case 0 % PTES
        fac = 10*100/1.1618/1.11; % This can be used to more easily set the mass flow to obtain a desired power output
        stH = 10 ;
        % This is the load scenario the plant is designed for
        Design_Load      = Load ;
        Design_Load.time = [stH/1;stH].*3600;  % time spent in each load period, s
        Design_Load.type = ["chg";"dis"];    % type of load period
        Design_Load.mdot = [fac*1;fac];  % working fluid mass flow rate, kg/s
        T0_inc    = 3.0 ; % Increment above ambient temperature that gas is cooled to
        
        if Loffdesign
            % This is the actual load profile that the plant meets
            if ~Lreadload
                Load.time = [stH;stH].*3600;      % time spent in each load period, s
                Load.type = ["chg";"dis"];    % type of load period
                %Load.mdot = mdotIN;      % working fluid mass flow rate, kg/s
                %T0_off    = T0IN;
                Load.mdot = [100.*fac;100.*fac];      % working fluid mass flow rate, kg/s
                T0_off    = [T0;T0] ;
            else
                fload     = './Data/load2.csv';
                fdat      = readmatrix(fload,'Range','A:B') ;
                T0_off    = fdat(:,1) ;
                Load.mdot = fdat(:,2) .* Design_Load.mdot(1) ;
                Load.type = readmatrix(fload,'Range','C','OutputType','string') ;
                Load.time = ones(numel(Load.mdot),1) * 3600.;
            end

        else
            Load = Design_Load ;
        end
                
    case 1 % Heat pump
        Design_Load      = Load ;
        stH              = 10 ;
        Design_Load.time = stH*3600;                  % time spent in each load period, s
        Design_Load.type = "chg";                     % type of load period
        Design_Load.mdot = 1000;                        % working fluid mass flow rate, kg/s
        T0_inc           = 5.0 ; % Increase in ambient temperature
        
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
        Design_Load.time = [0;10].*3600;                  % time spent in each load period, s
        Design_Load.type = ["sol";"dis"];                     % type of load period
        Design_Load.mdot = [0;10];                        % working fluid mass flow rate, kg/s
        T0_inc           = 5;  % Increase in ambient temperature
        
        if Loffdesign
            % This is the actual load profile that the plant meets
            Load.time = [0;10].*3600;      % time spent in each load period, s
            Load.type = ["sol";"dis"];    % type of load period
            Load.mdot = [0;10];      % working fluid mass flow rate, kg/s
            T0_inc    = 0;  % Increase in ambient temperature
        else
            Load = Design_Load ;
        end
        
    case 3 % JB charge, Rankine discharge
        fac = 100.0/1.3109 ; % This can be used to more easily set the mass flow to obtain a desired power output 
        
        % This is the load scenario the plant is designed for
        Design_Load      = Load ;
        Design_Load.time = [10/0.864;4;10;10].*3600;          % time spent in each load period, s
        Design_Load.type = ["chg";"str";"ran";"ran"];   % type of load period
        Design_Load.mdot = [10*fac;0;1.0*fac;1.0*fac];  % working fluid mass flow rate, kg/s
        Design_Load.options.useCold = [0;0;1;0];        % Use cold stores during Rankine discharge? This should be set to 0 for design cases of retrofits.
        
        T0_inc    = 5.0 ; % Increase in ambient temperature
        
        if Loffdesign
            % This is the actual load profile that the plant meets
            fac = 100 ;
            Load.time = [10;4;10;10].*3600;         % time spent in each load period, s
            Load.type = ["chg";"str";"ran";"ran"];  % type of load period
            %Load.mdot = [10*fac;0;1*fac;1*fac];     % working fluid mass flow rate, kg/s
            Load.options.useCold = [0;0;1;0];        % Use cold stores during Rankine discharge?
            %T0_off    = [T0-0;T0-0;T0-0;T0-0] ;
            Load.mdot = mdotIN;      % working fluid mass flow rate, kg/s
            T0_off    = T0IN;

        else
            Load = Design_Load ;
        end
        
    case 7 % Electric heater charge, Rankine discharge
        fac = 1.0 ; % This can be used to more easily set the mass flow to obtain a desired power output 
        
        % This is the load scenario the plant is designed for
        Design_Load      = Load ;
        Design_Load.time = [10;4;10].*3600;          % time spent in each load period, s
        Design_Load.type = ["sol";"str";"ran"];      % type of load period
        Design_Load.mdot = [10*fac;0;1*fac];         % working fluid mass flow rate, kg/s
        Design_Load.options.useCold = [0;0;0];       % Use cold stores during Rankine discharge? This should be set to 0 for design cases of retrofits.
        
        if Loffdesign
            % This is the actual load profile that the plant meets
            Load.time = [10;4;10].*3600;         % time spent in each load period, s
            Load.type = ["sol";"str";"ran"];     % type of load period
            Load.mdot = [10*fac;0;1*fac];        % working fluid mass flow rate, kg/s
            Load.options.useCold = [0;0;0];      % Use cold stores during Rankine discharge?
            T0_inc    = 0.0 ; % Increase in ambient temperature
        else
            Load = Design_Load ;
        end
end
Design_Load.num  = numel(Design_Load.time);
Design_Load.ind  = (1:Design_Load.num)';
Load.num  = numel(Load.time);
Load.ind  = (1:Load.num)';

% Hot storage tanks
switch PBmode
    case 0
        fHname  = 'SolarSalt';  % fluid name
        TH_dis0 = 300 + 273.15; % initial temperature of discharged hot fluid, K
        MH_dis0 = 1e9;          % initial mass of discharged hot fluid, kg
        TH_chg0 = 570 + 273.15; % initial temperature of charged hot fluid, K
        MH_chg0 = 0.00*MH_dis0; % initial mass of charged hot fluid, kg
        
        switch Load.mode
            case {0,1,2}
                % Cold storage tanks
                fCname  = 'Methanol'; % fluid name
                TC_dis0 = T0;           % initial temperature of discharged cold fluid, K
                MC_dis0 = 20e6;          % initial mass of discharged cold fluid, kg
                TC_chg0 = T0-50;        % initial temperature of charged cold fluid, K
                MC_chg0 = 0.00*MC_dis0; % initial mass of charged cold fluid, kg
                
            case 3
                % Cold storage tanks
                fCname  = 'Water';      % cold fluid name
                TC_dis0 = T0;           % initial temperature of discharged cold fluid, K
                MC_dis0 = 1e9;          % initial mass of discharged cold fluid, kg
                TC_chg0 = 273.15+1;     % initial temperature of charged cold fluid, K
                MC_chg0 = 0.00*MC_dis0; % initial mass of charged cold fluid, kg
                HTFname = 'INCOMP::MEG2[0.56]'; % secondary heat transfer fluid name
                HTF     = fluid_class(HTFname,'SF','TAB',NaN,Load.num,30); % Storage fluid
        end
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
    otherwise
        error('Not implemented')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COST MODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CCMPmode = [10,12,15,16] ; % Charging compressor cost mode. Removed mode 13.
CEXPmode = [41,43,44] ; % Charging expander cost mode. Removed mode 42.
DCMPmode = [10,12,15] ; % Discharging compressor cost mode. Removed mode 13.
DEXPmode = [41,43,44] ; % Discharging expander cost mode

PMPmode = [60,61,62]; % Pump cost mode
FANmode = [72,70,71]; % Fan cost mode. Removed mode 72.

hotHXmode = [1,3,4,5,10,11]; % Heat exchanger - hot cost mode. Removed mode 2.
cldHXmode = [1,3,4,5,10,11]; % Heat exchanger - cold cost mode. Removed mode 2.
rcpHXmode = [1,3,4,5,10,11]; % Heat exchanger - recuperator cost mode. Removed mode 2.
rejHXmode = [33,32]; % Heat exchanger - rejection cost mode. Removed mode 31. [33,32,30]

GENmode = [2,1,3] ; % Motor-generator

HTmode.tankmode  = [5,1,2,3,4,6] ; % Cost mode for hot tank container cost
HTmode.fld_cost  = [0.8,0.5,1.3] ; % Hot tank fluid cost, $/kg solarsalt : [0.8,0.5,1.3]. Mineral oil: [1.6,0.5,2.5]. Therminol: [3.1,0.5,2.5]. Chloride Salt: [0.64,0.4,1.1]
HTmode.ins_cost  = [30,5,50] ; % Insulation material, %/kg
HTmode.ins_k     = 0.08 ; % Thermal conductivity of insulation
HTmode.ins_rho   = 150 ; % Density of insulation
HTmode.tau       = 200 ; % Number of days before all heat leaks out of tank
HTmode.AR        = 1.0 ; % Aspect ratio (L/D) of tank
HTmode.over_fac  = 1.1 ; % How much larger is inner tank volume than the fluid volume

CTmode.tankmode  = [5,1,2,3,4,6] ; % Cost mode for cold tank container cost
CTmode.fld_cost  = [0.01,0.05,0.1] ; % Cold tank fluid cost, $/kg. Water: [0.01,0.05,0.1]. Methanol [0.3,0.1,1]
CTmode.ins_cost  = [30,5,50] ; % Insulation material, %/kg
CTmode.ins_k     = 0.08 ; % Thermal conductivity of insulation
CTmode.ins_rho   = 150 ; % Density of insulation
CTmode.tau       = 200 ; % Number of days before all heat leaks out of tank
CTmode.AR        = 1.0 ; % Aspect ratio (L/D) of tank
CTmode.over_fac  = 1.1 ; % How much larger is inner tank volume than the fluid volume

ATmode.tankmode = 0 ;
ATmode.fld_cost  = 0 ; % Cold tank fluid cost, $/kg
ATmode.ins_cost  = 0 ; % Insulation material, %/kg

steamC = [500,750,1000]; % Price of a steam turbine $/kW

% Set fluids. 'WF' or 'SF' indicates working fluid or storage fluid. 'CP'
% or 'TAB' indicate CoolProp or Tabular reading modes. 'backend' is used by
% CoolProp (either 'HEOS', 'TTSE' or 'BICUBIC&HEOS'). 'HEOS' is the most
% accurate method but is very slow. 'BICUBIC&HEOS' is recommended over
% 'TTSE' for speed and accuracy. 'num' indicates number of preallocated
% elements in state arrays.
gas = fluid_class('Nitrogen','WF','CP','BICUBIC&HEOS',Load.num,30);
%{
% Set up an ideal gas - should run faster

% gas_temp = fluid_class('Nitrogen','WF','CP','BICUBIC&HEOS',Load.num,30);
% dat.T0   = 300 ;
% dat.P0   = 1e5 ;
% dat.cp   = RPN('PT_INPUTS',dat.P0,dat.T0,'CPMASS',gas_temp) ;
% dat.cv   = RPN('PT_INPUTS',dat.P0,dat.T0,'CVMASS',gas_temp) ;
% dat.mu0  = 17.81e-6 ;
% dat.TVref = 300.55 ;
% dat.S     = 111 ;
% dat.k     = 0.026 ;
% 
% gas      = fluid_class('Nitrogen','WF','IDL',dat,Load.num,30);

%}
if any(Load.mode==[3,7])
    % 'TTSE' interpolation is NOT recommended for steam when reading values
    % close to the saturation curve. Use 'HEOS' or 'BICUBIC&HEOS'
    steam = fluid_class('Water','WF','CP','BICUBIC&HEOS',Load.num,30);
end

% Save copy of input file in "Outputs" folder
copyfile(['./PTES_scripts/',mfilename,'.m'],'./Outputs/')
