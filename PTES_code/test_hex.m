%%% To run this script, place it inside the PTES_code folder, alongside the
%%% PTES.m file, and run the first section of the PTES.m file to gain
%%% access to the different folders and load CoolProp.

% Clear variables
clear

% Set fluids
steam  = fluid_class('Water','WF','CP','TTSE',1,5);
fluidH = fluid_class('SolarSalt','SF','TAB',NaN,1,5);

% Define heat exchanger geometry (shell-and-tube)
% 1 refers to the tube side, 2 refers to the shell side
HX.shape     = 'circular';
HX.sigma     = 1e8;        % Maximum allowable stress, Pa
HX.L         = 3.0;        % Tube length, m
HX.D1        = 0.5e-2;     % Tube diameter, m
HX.t1        = 0.1*HX.D1;  % Tube thickness, m
HX.AfT       = 1.0;        % Total flow area, m2
HX.AfR       = 1.00;       % Ratio of flow areas, Af2/Af1, -

% Code settings
HX.NX  = 100;               % Number of sections (grid)

iL = 1;
i  = 1;
steam.state(iL,i).p = 100*1e5;
steam.state(iL,i).T = 400;
steam.state(iL,i).mdot = 10;
[steam] = update(steam,[iL,i],1);
fluidH.state(iL,1).p = 1e5;
fluidH.state(iL,1).T = 830;
[fluidH] = update(fluidH,[iL,1],1);
[steam,fluidH,i,~] = hex_TQA(steam,[iL,i],fluidH,[iL,1],HX,'hex',2,1.20);