%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTES
% This a thermodynamic model of a pumped thermal energy storage cycle. It
% is based on the Joule-Brayton cycle and employs liquid storage media.
% Author: Pau Farres-Antunez (pf298@cam.ac.uk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop

% SET INPUTS
PTES_INPUTS

% Open the output files and print the headers
PTES_MANAGE_FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic % start timer
for ix = 1:1
    %%% RUN CYCLE LOOP %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for icrv = 1:Ncrv
        for ipnt = 1:Npnt
            % Select Write Mode (if WM=1, then write cycle for plotting)
            if all([icrv,ipnt] == 1), WM=1; else, WM=0; end
            
            % Set multi_run variables
            if multi_run==1, PTES_SET_MULTI_RUN; end            
            
            % Reinitialise arrays (gas, fluids and tanks) to zero and do other
            % preliminary tasks
            PTES_INITIALISE
            
            for iL = 1:Load.num
                switch Load.type(iL)
                    case 'chg'
                        PTES_CHARGE
                        
                    case 'str'
                        PTES_TANKS_STORAGE
                        
                    case 'dis'
                        PTES_DISCHARGE
                        
                    case 'ran'
                        RANK_DISCHARGE
                        
                    case 'sol'
                        PTES_SOLAR_TANKS
                end
            end
            
            if optimise % obtain optimal PRr
                error('not implemented')
                f = @(PRr) ptes_discharge_function(gas, fluidH, fluidC, HT, CT, environ,...
                    T0, T1, pbot, PRr, PRch, Load.mdot(iL), Nc_dis, Ne_dis,...
                    eta, eff, ploss, Load.time(iL), mode);
                
                [PRr,ineff,xv,yv,iter] = golden_search(f,PRr_min,PRr_max,0.005,'Min',100);
            end
            
            % Compute energy balance
            PTES_ENERGY_BALANCE
            
            if multi_run
                PTES_PRINT_MULTI_RUN
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
toc %stop timer

%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if make_plots
    switch Load.mode
        case {0,3} % PTES
            %PTES_PLOT_HEXS
            if optimise
                PTES_PLOT_GOLDEN_SEARCH
            end
        case 2 % Heat engine only
            if optimise
                PTES_PLOT_GOLDEN_SEARCH
            end
    end
    PTES_PLOT_CYCLE
    PTES_PLOT_LOSSES
    if multi_run
        PTES_PLOT_MULTI_RUN %#ok<*UNRCH>
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% FINISH PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close files, save plots and release CoolProp AbstractStates
PTES_FINISH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%