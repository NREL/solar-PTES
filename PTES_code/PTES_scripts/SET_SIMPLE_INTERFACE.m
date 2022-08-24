% Simple interface
% The aim of this script is to provide a mimimum number of inputs to be
% able to run the solar-PTES code. This provides a simple introduction to
% some of the main parameters and is therefore good for beginners. It is
% also a suitable tool to enable integration with other coding frameworks -
% e.g. if the solar-PTES code is to be called by another program. 

% The code below has been set up to model a molten-salt Joule-Brayton PTES
% cycle for a single design-point calculation. 

% A well-informed user can modify or create a new version of this file to
% provide a simple interface to access other features of the code, if
% desirable. 

% INPUTS CORRESPONDING TO 'SYSTEM DESIGN' TAB
heater_mult = 2 ;            % How much faster heat pump operates compared to heat engine
stH         = 10 ;            % Hours of storage, h (really discharging duration)
Tmax        = 570 + 273.15;   % Compressor outlet temperature (which is > Hot storage hot temperature, K)
TH_dis0     = 370.0+273.15;   % Hot storage cold temperautre, K
pow_out     = 100 ;           % Cycle thermodynamic power, MW-e

% INPUTS CORRESPONDING TO 'LOCATION AND RESOURCE' TAB
T0         = 25 + 273.15;    % ambient temp, K
p0         = 1e5;            % ambient pressure, Pa

% Some other important inputs that do not get entered into SAM
eff        = 0.95;           % heat exchanger effectiveness
ploss      = 0.01;           % pressure loss in HEXs
eta        = 0.90;           % polytropic efficiency
fluid_name = 'Nitrogen';     % Power cycle working fluid
fHname     = 'SolarSalt';    % hot storage fluid name
fCname     = 'Methanol';     % cold storage fluid name
TC_dis0    = T0;             % initial temperature of discharged cold fluid, K. This corresponds to the cold tank hot temperature. At the moment, the cycle requires this to be ambient temperature.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















% The code below checks that the program will model a single
% design-point run of a molten-salt Joule-Brayton PTES cycle. 
% DO NOT EDIT.

if (SIiter == 0)
    % Guess mass flow rate
    massflow = pow_out * 1e6 / (1000 * ((554-250)-(25--67)));
else
    % Calculate mass flow rate from power
    SIconv = 100*abs(W_out_dis/t_dis/1e6 - pow_out) / pow_out ;
    massflow = massflow * pow_out / (W_out_dis/t_dis/1e6) ;
end


Design_Load.mdot = [heater_mult*massflow;massflow];
Load.mdot = Design_Load.mdot ;
Load0.mdot = Design_Load.mdot ;

Design_Load.time = [stH/heater_mult;stH].*3600;  % time spent in each load period, s
Load.time = Design_Load.time;
Load0.time = Design_Load.time ;


SIiter = SIiter + 1;


if x ~= 0
   error('SIMPLE INTERFACE ERROR: Set run_mode = 0 in main.m') 
end

if Load.mode  ~= 0
   error('SIMPLE INTERFACE ERROR: Set Load.mode  = 0 in INPUTS.m') 
end

if Loffdesign ~= 0 
   error('SIMPLE INTERFACE ERROR: Set Loffdesign = 0  in INPUTS.m') 
end

if Lreadload  ~= 0 
   error('SIMPLE INTERFACE ERROR: Set Lreadload  = 0  in INPUTS.m') 
end

if PBmode     ~= 0 
   error('SIMPLE INTERFACE ERROR: Set PBmode     = 0  in INPUTS.m') 
end

if multi_run   ~= 0
   error('SIMPLE INTERFACE ERROR: Set multi_run   = 0 in INPUTS.m') 
end

if Lmulti_mdot ~= 0
   error('SIMPLE INTERFACE ERROR: Set Lmulti_mdot = 0 in INPUTS.m') 
end

       