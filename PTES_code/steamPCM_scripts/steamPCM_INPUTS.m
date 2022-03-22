%**** INPUT PARAMETERS ****%
% PIPE SIZING
dp          =       0.02 ;          % Pipe diameter

% STORAGE SIZING
Pmax        =       1e6 ;            % Maximum power input, W    
tN          =       6 ;              % Nominal duration of charge, hours

% STEAM CHARGING CONDITIONS
PsatC       =       55e5;%20e5 ;          % Saturation pressure
xsC         =       0.6 ;          % Initial dryness fraction of steam

% STEAM DISCHARGING CONDITIONS
PsatD       =       35e5;%12e5 ;          % Saturation pressure
xsD         =       0.0 ;           % Initial dryness fraction of steam

% PCM
xpC         =       1.0 ;
xpD         =       0.0 ;
cpl         =       1.216e3 ;         % Specific heat capacity of liquid PCM, J/kg.K
cps         =       1.216e3 ;         % Specific heat capacity of solid PCM, J/kg.K
rhopl       =       1.92e3 ;           % Density of liquid PCM
rhops       =       1.92e3 ;         % Density of solid PCM
Lpcm        =       245e3 ;           % Latent heat of PCM, J/kg
kpcm        =       5 ;             % Thermal conductivity of PCM

XPend_chg   =       0.85 ;           % When PCM melted fraction exceeds this, end the charging cycle
XPend_dis   =       0.2 ;           % When PCM melted fraction is less than this, end the discharging cycle

% GRID AND TIME STEPS
CFL         =       10 ;           % Courant-Freidrich-Lewy number
Nx          =       100 ;          % Number of gridsteps
Nload       =       2 ;             % Number of loads (charge and discharge are separate loads)
Load        =       ["c";"d";"c"];  % 'c' - charge. 'd' - discharge, duh. 
PCMiterations =     2 ;             % Number of times to iterate PCM equation. 2 gives good results. 1 is satisfactory.

% Read input load data
Lreadload   =       true ;         % Is load data being read from a file? If not, use the data below
fload       =       '.\steamPCM_scripts\data\hourly_data_short_270C.csv' ;
SM          =       2.0 ;             % Solar multiple
dsg_Tin     =       230+273.15;%180 + 273.15 ;  % Temperature into the DSG solar field
dsg_Pin     =       55e5;%20e5 ;          % Pressure of DSG solar field

mdotC_fac   =       1.0 ;
mdotD_fac   =       1.0 ;

% ECONOMICS
Cpcm        =       750.0 ;       % PCM cost, $/m3
Csteel      =       4800;        % Steel cost, $/m3