% Determine Operating System
clc
clear all

%run_mode = 0: single run. Use INPUTS file to change the thermodynamic cycle
%run_mode = 1: Multi run. Set parameters in INPUTS file
%run_mode = 2: Optimizing PTES using NSGA2. Set Parameters in SET_OPTIMIZE_RUN
%run_mode = 3: Optimizing PTES using MOPSO. Set Parameters in SET_OPTIMIZE_RUN
run_mode = 1;

fname='PTES_optimize2';

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