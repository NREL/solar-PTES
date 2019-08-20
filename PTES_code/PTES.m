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
%dbstop if error
%dbclear all

% Set paths
addpath('./Classes/')
addpath('./Generic/')
addpath('./Functions/')
addpath('./PTES_scripts/')

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop

% SET INPUTS
PTES_INPUTS

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
        
        tic %start timer
        
        % Reinitialise arrays (gas, fluids and tanks) to zero and do other
        % preliminary tasks
        PTES_INITIALISE
        
        switch mode
            case 0 % PTES
                % Run charge cycle, compute state of storage tanks and run
                % discharge cycle
                PTES_CHARGE
                PTES_TANKS_STORAGE % no storage loss considered at the moment
                
                if optimise % obtain optimal PRr
                    f = @(PRr) ptes_discharge_function(gas, fluidH, fluidC, HT, CT, environ,...
                        T0, T1, pbot, PRr, PRch, mdot, Nc_dis, Ne_dis,...
                        eta, eff, ploss, t_ch, mode);
                    
                    [PRr,ineff,xv,yv,iter] = golden_search(f,PRr_min,PRr_max,0.005,'Min',100);
                end
                PTES_DISCHARGE
                
            case 1 % Heat pump only
                PTES_CHARGE
                
            case 2 % Heat engine only
                PTES_SOLAR_TANKS
                if optimise % obtain optimal PRr
                    f = @(PRr) ptes_discharge_function(gas, fluidH, fluidC, HT, CT, environ,...
                        T0, T1, pbot, PRr, PRch, mdot, Nc_dis, Ne_dis,...
                        eta, eff, ploss, t_ch, mode);
                    
                    [PRr,ineff,xv,yv,iter] = golden_search(f,PRr_min,PRr_max,0.005,'Min',100);
                end
                PTES_DISCHARGE
        end
        
        % Compute energy balance
        PTES_ENERGY_BALANCE
        
        toc %stop timer
        
        if multi_run
            switch mode
                case 0
                    chi = chi_PTES;
                    COP=0;
                    EFF=0;
                case 1
                    chi = chi_hot;
                    EFF=0;
                case 2
                    chi = chi_tot;
                    COP=0;
                    rhoP_ch=0;
            end
            fprintf(ID1,'%15.5g', PRch , eta, eff, ploss, chi,COP,EFF,...
                HT.B(2).T, CT.B(2).T, WR_dis, rhoE, rhoP_ch,...
                WL_comp, WL_exp, WL_hexs, WL_reject, WL_mix_liq, WL_tanks);
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


%%% FINISH PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close files and save copy of multi_run files
fprintf(1,'\n...END...\n');
fclose('all');
if ~multi_run
    diary off
    % Save plots
    if save_figs == 1
        save_fig(1,'./Outputs/T-s',0,0,0)
        %save_fig(3,'./Outputs/Golden_search',0,0,0)
        save_fig(8,'./Outputs/Losses',0,0,0)
    end
end

% Releasing CoolProp AbstractState
ierr = 0; buffer_size = 10;
herr= char((1:1:buffer_size));
for i0=1:1000
    calllib('coolprop','AbstractState_free',i0, ierr,herr,buffer_size)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%