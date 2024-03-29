%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTES
% This a thermodynamic model of a pumped thermal energy storage cycle. It
% is based on the Joule-Brayton cycle and employs liquid storage media.
% Author: Pau Farres-Antunez (pf298@cam.ac.uk)
%
% Sept 16 2019. This is an extension to the PTES model to see what
% modifications are required to develop a cycle that uses supercritical-CO2
% as the working fluid.
% Contributions: Josh McTigue (JoshuaDominic.McTigue@nrel.gov)
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
sCO2_PTES_INPUTS

% Open the output files and print the headers
PTES_MANAGE_FILES

%%% RUN CYCLE LOOP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for icrv = 1:Ncrv
    for ipnt = 1:Npnt        
        % Select Write Mode (if WM=1, then write cycle for plotting)
        if all([icrv,ipnt] == 1), WM=1; else, WM=0; end
        
        % Set multi_run variables
        if multi_run==1, PTES_SET_MULTI_RUN; end
        
        tic % start timer
        
        % Reinitialise arrays (gas, fluids and tanks) to zero and do other
        % preliminary tasks
        if Load.mode == 4
            SCO2_PTES_INITIALISE
        else
            PTES_INITIALISE
        end
        
        for iL = 1:Load.num
            switch Load.type(iL)
                case 'chg'
                    PTES_CHARGE
                    
                case 'str'
                    PTES_TANKS_STORAGE
                    
                case 'dis'
                    PTES_DISCHARGE
                    
                case 'sol'
                    PTES_SOLAR_TANKS
                    
                case 'chgCO2'
                    sCO2_PTES_CHARGE
                    
                case 'disCO2'
                    sCO2_PTES_DISCHARGE
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
        
        toc %stop timer
        
        if multi_run
            PTES_PRINT_MULTI_RUN
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if make_plots
    switch Load.mode
        case {0,4} % PTES
            PTES_WRITE_CHARGE
            PTES_WRITE_DISCHARGE
            %PTES_PLOT_HEXS
            if optimise
                PTES_PLOT_GOLDEN_SEARCH
            end
        case 1 % Heat pump only
            PTES_WRITE_CHARGE
        case 2 % Heat engine only
            PTES_WRITE_DISCHARGE
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
% Have to plot recuperators here ...
%{
switch Load.mode
    case 0
        if Nrcp == 0
            plot_hex(gas,[1,2],fluidH,[1,1],100,30);
            if Nhot==2, plot_hex(gas,[1,3],fluidH2,[1,1],100,32); end
            
            plot_hex(gas,[1,6],fluidC,[1,1],100,40);
            if Ncld==2, plot_hex(gas,[1,7],fluidC2,[1,1],100,42); end
        elseif Nrcp == 2
            plot_hex(gas,[1,7],gas,[1,3],100,20);
            plot_hex(gas,[3,5],gas,[3,1],100,22); % Recuperator
        elseif Nrcp == 2
            plot_hex(gas,[1,9],gas,[1,3],100,20);
            plot_hex(gas,[1,8],gas,[1,4],100,21);
            plot_hex(gas,[3,8],gas,[3,1],100,22); % Recuperator
            plot_hex(gas,[3,6],gas,[3,2],100,23); % Recuperator
        end

        if Lrcmp
            % Plot hot storage - currently set up for sCO2-recomp cycle
            plot_hex(gas,[1,2],fluidH,[1,1],100,30);
            plot_hex(gas,[3,9],fluidH,[3,1],100,31);

            % Plot coldhot storage - currently set up for sCO2-recomp cycle
            plot_hex(gas,[1,7],fluidC,[1,1],100,40);
            plot_hex(gas,[3,4],fluidC,[3,1],100,41);
        end
    case 2
end
%}

%%% FINISH PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close files, save plots and release CoolProp AbstractStates
PTES_FINISH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%