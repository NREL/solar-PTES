%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAES
% This a thermodynamic model of a liquid air energy storage cycle.
% Author: Pau Farres-Antunez (pf298@cam.ac.uk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clear;

% Enter debugging mode if an error occurs
%dbstop if error
%dbclear all

% Set paths
addpath('./Classes/')
addpath('./Generic/')
addpath('./Functions/')
addpath('./LAES_scripts/')

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop

% READ INPUTS FROM EXTERNAL FILE
LAES_INPUTS

% Open the output files and print the headers
LAES_MANAGE_FILES

% Declare arrays of gas states and fluid streams
LAES_DECLARE_ARRAYS

%%% RUN CYCLE LOOP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for icrv = 1:Ncrv
    for ipnt = 1:Npnt        
        % Select Write Mode (if WM=1, then write cycle for plotting)
        if all([icrv,ipnt] == 1), WM=1; else, WM=0; end
        
        % Set multi_run variables
        if multi_run==1, LAES_SET_MULTI_RUN; end
        
        tic %start timer
        
        % Reinitialise arrays (gas, fluids and tanks) to zero and do other
        % preliminary tasks
        LAES_INITIALISE
        
        % Run charge cycle, compute state of storage tanks and run
        % discharge cycle
        LAES_CHARGE
        LAES_TANKS_STORAGE % no storage loss considered at the moment
        LAES_DISCHARGE
        
        LAES_WRITE_CHARGE
        LAES_WRITE_DISCHARGE
        LAES_PLOT_CYCLE
        
        % Compute energy balance
        LAES_ENERGY_BALANCE
        
        toc %stop timer
        
        if multi_run
            fprintf(ID1,'%15.5g', PR , eta, eff, ploss, chi_PTES*100,...
                HT.B(2).T, CT.B(2).T, WR_dis, rhoE, rhoP_ch,...
                WL_turbo, WL_hexs, WL_reject, WL_mix_liq, WL_tanks);
            fprintf(ID1,'\n');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if make_plots
    switch mode
        case 0 % PTES
            LAES_WRITE_CHARGE
            LAES_WRITE_DISCHARGE
            LAES_PLOT_LOSSES
            %PLOT_HEXS
            if multi_run
                LAES_PLOT_MULTI_RUN %#ok<*UNRCH>
            end
        otherwise
            error('***not implemented***')
    end
    LAES_PLOT_CYCLE
    LAES_PLOT_LOSSES
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% FINISH PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close files and save copy of multi_run files
fprintf(1,'\n...END...\n');
fclose('all');
if ~multi_run
    diary off
end

% Save plots
if savefig == 1
    save_fig(1,'./Outputs/T-s',0,0,0)
    %save_fig(4,'./Outputs/Golden_search',0,0,0)
    save_fig(8,'./Outputs/Losses',0,0,0)
end

% Releasing CoolProp AbstractState
ierr = 0; buffer_size = 10;
herr= char((1:1:buffer_size));
for i0=1:1000
    calllib('coolprop','AbstractState_free',i0, ierr,herr,buffer_size)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%