% Determine Operating System
clc
clear all

%run_mode = 0: single run. Use INPUTS file to change the thermodynamic cycle
%run_mode = 1: Multi run. Set parameters in INPUTS file
%run_mode = 2: Optimizing PTES using NSGA2 (Two Objectives). Set Parameters in SET_OPTIMIZE_RUN
%run_mode = 3: Optimizing PTES using MOPSO (Three Objectives). Set Parameters in SET_OPTIMIZE_RUN
 run_mode = 1;
 
 %If optimization routines are called set the number of objectives here
 objectives = 3; %Number of objectives = 2 or 3.
 csvwrite('./Data/objectives',objectives);
 
 
% Mode 0: Joule-Brayton PTES cycle
% Mode 1: Joule-Brayton heat pump (charge only)
% Mode 2: Joule-Brayton heat engine (discharge only)
% Mode 3: Joule-Brayton heat pump with a steam-Rankine heat engine
% Mode 4: sCO2 PTES cycle
% Mode 5: sCO2 recompression cycle (discharge only)
% Mode 6: sCO2 'time-shifted' recompression cycle
% Mode 7: Steam-Rankine heat engine (discharge only)
%
% Mode 20: PTES-LAES. Combined cycle energy storage (CCES)
 
 mode  = 0 ;
 csvwrite('./Data/mode',mode);

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
%Current objective functions are:
%objective 1 = 1 - roundtrip efficiency
%Objective 2 = LCOS
%Objective 3 = Capital Cost
%To change the objectives do the following
%Go to CONSTRAINT_AND_OUTPUT file and ERROR file to change f1 f2 and f3 accordingly.
%Current input variables for mode 0: TH_dis0, PRr and eff
%Current input variables for mode 3: TH_dis0, Ne_ch and eff
%Go to SET_OPTIMIZE_RUN to change the input variables