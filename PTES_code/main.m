%clc
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

% Load CoolProp library (low-level interface)
load_coolprop

%run_mode = 0: Run PTES cases (may be single run, multiple run, or off-design cases). Use INPUTS file to change the thermodynamic cycle
%run_mode = 1: Optimizing PTES using MOPSO (Three Objectives). Set Parameters in SET_OPTIMIZE_RUN
run_mode = 0;
fname='PTES_optimize';

switch run_mode
    case 0
        PTES_optimize(0, 0, 0);
    case 1
        Optimize_MOPSO;
end
