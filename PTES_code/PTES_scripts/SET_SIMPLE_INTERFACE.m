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

eff        = 0.95;           % heat exchanger effectiveness
ploss      = 0.01;           % pressure loss in HEXs
T0         = 25 + 273.15;    % ambient temp, K
p0         = 1e5;            % ambient pressure, Pa
Tmax       = 570 + 273.15;   % maximum temp at compressor outlet, K
eta        = 0.90;           % polytropic efficiency
stH        = 10 ;            % Discharging duration, h
fluid_name = 'Nitrogen';     % Power cycle working fluid
fHname     = 'SolarSalt';    % hot storage fluid name
TH_dis0    = 250.0+273.15;   % initial temperature of discharged hot fluid, K
fCname     = 'Methanol';     % cold storage fluid name
TC_dis0    = T0;             % initial temperature of discharged cold fluid, K
                
