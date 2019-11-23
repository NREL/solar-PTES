%%% To run this script, place it inside the PTES_code folder, alongside the
%%% PTES.m file

%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clear;

% Enter debugging mode if an error occurs
dbstop if error
%dbclear all

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% INPUTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose scenario
% 1 = Helium and Helium
% 2 = SolarSalt and Water
% 3 = CO2 and Water
scenario = 1;

% Set indices
iL = 1; i1 = 1; i2 = 1;

% Declare fluids and specify initial conditions
switch scenario
    case 1
        % Helium (hot, high pressure)
        F1 = fluid_class('Helium','WF','CP','TTSE',1,5);
        F1.state(iL,i1).p = 100e5;
        F1.state(iL,i1).T = 600;
        F1.state(iL,i1).mdot = 10;
        
        % Helium (cold, low pressure)
        F2 = fluid_class('Helium','WF','CP','TTSE',1,5);
        F2.state(iL,i2).p = 20e5;
        F2.state(iL,i2).T = 300;
        F2.state(iL,i2).mdot = 10;
        
        % Set hex_mode
        hex_mode = 0; % Both mass flow rates specified
        var = 0;
        
    case 2        
        % Solar Salt
        F1 = fluid_class('SolarSalt','SF','TAB',NaN,1,5);
        F1.state(iL,i1).p = 1e5;
        F1.state(iL,i1).T = 830;
        
        % Water
        F2 = fluid_class('Water','WF','CP','TTSE',1,5);
        F2.state(iL,i2).p = 100*1e5;
        F2.state(iL,i2).T = 400;
        F2.state(iL,i2).mdot = 10;
        
        % Set hex_mode
        hex_mode = 2; % Compute mass flow rate of hot fluid with var=mH*CpH/(mC*CpC)
        var = 1.20;     
        
    case 3        
        % CO2
        F1 = fluid_class('CarbonDioxide','WF','CP','TTSE',1,5); % working fluid
        F1.state(1,i1).p = 85e5;
        F1.state(1,i1).T = 380;
        F1.state(1,i1).mdot = 0.75;
        
        % Water
        F2 = fluid_class('Water','SF','CP','TTSE',1,5); % storage fluid       
        F2.state(1,i2).p = 5e5;
        F2.state(1,i2).T = 300;
        F2.state(1,i2).mdot = 1;
        
        % Set hex_mode
        hex_mode = 0; % Both mass flow rates specified
        var = 0;
        
    otherwise
        error('not implemented')
end
% Update fluid states
[F1] = update(F1,[iL,i1],1);
[F2] = update(F2,[iL,i2],1);

% Specify HEX geometry
method = 'automatic';
switch method
    case 'manual'
        % Define heat exchanger geometry (shell-and-tube)
        % 1 refers to the tube side, 2 refers to the shell side
        HX.shape     = 'circular';
        HX.sigma     = 1e8;        % Maximum allowable stress, Pa
        HX.t1_min    = 0.1e-3;     % Minimum tube thickness, m
        HX.t1_D1_min = 0.05;       % Minimum tube thickness-to-diameter ratio
        HX.L         = 3.0;        % Tube length, m
        HX.D1        = 5e-3;       % Tube diameter, m
        HX.AfT       = 0.25;       % Total flow area, m2
        HX.AfR       = 1.00;       % Ratio of flow areas, Af2/Af1, -
        
    case 'automatic'
        % Obtain geometric parameters based on performance objectives,
        % using analytical solutions. It should be expected that objectives
        % will be met accurately (only) when using fluids with small
        % variations of thermophysical properties.
        eff   = 0.97;
        ploss = 0.01;
        D     = 1e-2;
        [HX]  = set_hex_geom(F1,[iL,i1],F2,[iL,i2],eff,ploss,D,hex_mode,var);        
end

% Hex code settings
HX.NX = 100;               % Number of sections (grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% COMPUTE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;
if n>1
    mdot1 = F1.state(iL,i1).mdot*linspace(0.5,2.0,n);
    mdot2 = F2.state(iL,i1).mdot*linspace(0.5,2.0,n);
else
    mdot1 = F1.state(iL,i1).mdot;
    mdot2 = F2.state(iL,i1).mdot;
end
NU_eff  = zeros(size(mdot1));
AN_eff  = zeros(size(mdot1));
NU_DppH = zeros(size(mdot1));
AN_DppH = zeros(size(mdot1));
NU_DppC = zeros(size(mdot1));
AN_DppC = zeros(size(mdot1));
for im = 1:n
    F1.state(iL,i1).mdot = mdot1(im);
    F2.state(iL,i1).mdot = mdot2(im);
    % Run HEX code
    [F1,F2,~,~,HX] = hex_TQA(F1,[iL,i1],F2,[iL,i2],HX,'hex',hex_mode,var);
    
    % Compare numerical results (NU) with analytical results (AN)
    NU_eff(im)  = HX.QS(HX.NX+1)/HX.QMAX;
    NU_DppH(im) = (HX.H.pin-HX.H.p(1))./HX.H.pin;
    NU_DppC(im) = (HX.C.pin-HX.C.p(HX.NX+1))./HX.C.pin;
    [AN_eff(im),AN_DppH(im),AN_DppC(im)] = hex_analytic(F1,[iL,i1],F2,[iL,i2],HX);
    fprintf(1,'\n      Numerical  Analytical\n')
    fprintf(1,'Eff  = %8.3f   %9.3f\n',NU_eff(im),AN_eff(im))
    fprintf(1,'DppH = %8.5f   %9.5f\n',NU_DppH(im),AN_DppH(im))
    fprintf(1,'DppC = %8.5f   %9.5f\n',NU_DppC(im),AN_DppC(im))
end
if n>1
    figure(25)
    plot(mdot1,NU_eff,'o',mdot1,AN_eff)
    xlabel('Mass flow rate, kg/s')
    ylabel('Effectiveness')
    legend('numerical','analytical')
    figure(26)
    semilogy(mdot1,NU_DppH,'o',mdot1,AN_DppH,mdot1,NU_DppC,'o',mdot1,AN_DppC)
    xlabel('Mass flow rate, kg/s')
    ylabel('$ \Delta p / p $')
    legend('FH numerical','FH analytical','FC numerical','FC analytical','Location','Best')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots
plot_hex(HX,20,'C');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = get_hex_area(fluidH, indH, fluidC, indC, HX, xvars)

% Import variable parameters from "xvars" and into "HX" structure
HX.D1  = xvars(1);
HX.L   = xvars(2);
HX.AfT = xvars(3);
HX.AfR = xvars(4);

% Import fluid.state
stateH = fluidH.state(indH(1),indH(2));
stateC = fluidC.state(indC(1),indC(2));

% Declare the two fluid streams
H = stream; H.pin = stateH.p;
C = stream; C.pin = stateC.p;

% Compute derived geometric parameters and mass fluxes
[~, ~, HX] = shell_and_tube_geom(C, H, HX);
result = HX.A2;

end

function [c,ceq] = hex_constraints(F1, ind1, F2, ind2, HX, stage_type, hex_mode, var, obj_eff, obj_ploss, xvars)

% Import variable parameters from "xvars" and into "HX" structure
HX.D1  = xvars(1);
HX.L   = xvars(2);
HX.AfT = xvars(3);
HX.AfR = xvars(4);

[~,~,~,~,~,out] = hex_TQA(F1,ind1,F2,ind2,HX,stage_type,hex_mode,var);

c = out - [1-obj_eff, obj_ploss, obj_ploss];

ceq = 0;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%