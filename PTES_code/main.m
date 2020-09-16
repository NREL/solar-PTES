% Determine Operating System
clc
clear all

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/','./Other/','./LIB/','./Optimization')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\','.\Other\','.\LIB\','.\Optimization')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop

%run_mode = 0: single run. Use INPUTS file to change the thermodynamic cycle
%run_mode = 1: Multi run. Set parameters in INPUTS file
%run_mode = 2: Optimizing PTES using NSGA2 (Two Objectives). Set Parameters in SET_OPTIMIZE_RUN
%run_mode = 3: Optimizing PTES using MOPSO (Three Objectives). Set Parameters in SET_OPTIMIZE_RUN
run_mode = 3;
fname='PTES_optimize';

switch run_mode
    case 0
        feval(fname, 0);
    case 1
        feval(fname, 1);
    case 2
        Optimize_NSGA2;
    case 3
        Optimize_MOPSO;
end
