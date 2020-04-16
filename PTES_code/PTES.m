%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTES
% This code employs thermodynamic and economic models to predict the
% performance and cost of different pumped thermal energy storage cycles.
%
% Authors: Pau Farres-Antunez and Joshua Dominic McTigue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace and global variables
clear;
clear global;

% Enter debugging mode if an error occurs
dbstop if error
%dbclear all

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/','./Other/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\','.\Other\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop

% SET INPUTS
INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic % start timer

for ix = 1:1
    %%% RUN CYCLE LOOP %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for icrv = 1:Ncrv
        for ipnt = 1:Npnt
            % Select Write Mode (if WM=1, then write cycle for plotting)
            if all([icrv,ipnt] == 1), WM=1; else, WM=0; end
            
            % Set multi_run variables
            if multi_run==1, SET_MULTI_RUN; end            
            
            % Reinitialise arrays (gas, fluids and tanks) to zero and do
            % other preliminary tasks
            INITIALISE
            
            for iix = 1:(Loffdesign+1)
                fprintf(['\n',line,txt(iix,:),line,'\n'])
            
                for iL = 1:Load.num
                    switch Load.type(iL)
                        case 'chg'
                            JB_CHARGE
                            
                        case 'dis'
                            JB_DISCHARGE
                            
                        case 'ran'
                            RANK_DISCHARGE
                            
                        case 'chgCO2'
                            sCO2_CHARGE
                            
                        case 'disCO2'
                            sCO2_DISCHARGE
                            
                        case 'rcmpCO2'
                            sCO2_RECOMP
                            
                        case 'chgTSCO2'
                            TSCO2_CHARGE
                            
                        case 'disTSCO2'
                            TSCO2_DISCHARGE
                            
                        case 'str'
                            TANKS_STORAGE
                            
                        case 'sol'
                            SOLAR_TANKS
                    end
                end
                
                if Loffdesign && iix==1
                    SET_DESIGN
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
            ENERGY_BALANCE
            
            % Evaluate the system cost
            PTES_ECONOMICS
            
            % Save results from multi_run call
            if multi_run, PRINT_MULTI_RUN; end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
toc %stop timer


%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if make_plots    
    PLOT_CYCLE
    PLOT_LOSSES
    PLOT_COSTS
    if multi_run
        PLOT_MULTI_RUN %#ok<*UNRCH>
    end
end
if make_hex_plots
    PLOT_HEXS
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% FINISH PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close files, save plots and release CoolProp AbstractStates
FINISH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%