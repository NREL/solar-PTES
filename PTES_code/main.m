% Determine Operating System
clc
clear all

%run_mode = 0: single run. Use INPUTS file to change the thermodynamic cycle
%run_mode = 1: Multi run. Set parameters in INPUTS file
%run_mode = 2: Optimizing PTES using NSGA2 (Two Objectives). Set Parameters in SET_OPTIMIZE_RUN
%run_mode = 3: Optimizing PTES using MOPSO (Three Objectives). Set Parameters in SET_OPTIMIZE_RUN
 run_mode = 2;

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

%Inputs for optimzers are in resprctive optimization files
%Current number of objectives is set to 3
%To change it two objectives do the following
%Go to Optimize_NSGA2, change M to 2. No need for any changes in
%Optimize_MOPSO.
%Go to CONSTRAINT_AND_OUTPUT file and ERROR file to change 'fit' accordingly.