function [fit, err, extra]=PTES_optimize(x,opt_par,Nobjs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTES
% This code employs thermodynamic and economic models to predict the
% performance and cost of different pumped thermal energy storage cycles.
%
% Authors: Pau Farres-Antunez and Joshua Dominic McTigue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dbstop if error
%dbclear all

% Set properties for plots
set_graphics

% SET INPUTS
INPUTS

%Check for optimization calling
if Nobjs > 1 %First variable in optimization is always a Temperature
   SET_OPTIMIZE_RUN
end

tic % start timer
try

    %%% RUN CYCLE LOOP %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for icrv = 1:Ncrv
        for ipnt = 1:Npnt
            % Select Write Mode (if WM=1, then write cycle for plotting)
            if all([icrv,ipnt] == 1), WM=1; else, WM=0; end
            
            % Set multi_run variables
            if multi_run==1, SET_MULTI_RUN; end
            
            SIiter = 0;
            SIconv = 100;
            while ((SIiter < 10 && simple_interface==1 && SIconv > 1 ) || (simple_interface==0 && SIiter < 1))
            
            % If using the simple interface, set these variables
            if simple_interface==1
                SET_SIMPLE_INTERFACE
            else
                SIiter = 1;
            end
            
            
            % Reinitialise arrays (gas, fluids and tanks) to zero and do
            % other preliminary tasks
            INITIALISE
            
            for iix = 1:(Loffdesign+1)
                %fprintf(['\n',line,txt(iix,:),line,'\n'])
                iL=1;
                while iL <= Load.num
                    switch Load.type(iL)
                        case 'chg'
                            JB_CHARGE
                            iL=iL+1;
                        case 'chgPB'
                            JB_CHARGE_PB
                            iL=iL+1;
                        case 'dis'
                            JB_DISCHARGE%_alt_Qrej
                            iL=iL+1;
                        case 'disPB'
                            JB_DISCHARGE_PB
                            
                        case 'ran'
                            RANK_DISCHARGE
                            iL=iL+1;
                        case 'chgCO2'
                            sCO2_CHARGE
                            iL=iL+1;
                        case 'disCO2'
                            sCO2_DISCHARGE
                            iL=iL+1;
                        case 'rcmpCO2'
                            sCO2_RECOMP
                            iL=iL+1;
                        case 'chgTSCO2'
                            TSCO2_CHARGE
                            iL=iL+1;
                        case 'disTSCO2'
                            TSCO2_DISCHARGE
                            iL=iL+1;
                        case 'str'
                            TANKS_STORAGE
                            iL=iL+1;
                        case 'sol'
                            SOLAR_TANKS
                            iL=iL+1;
                        case 'chgCC'
                            CHARGE_CCES
                            iL=iL+1;
                        case 'disCC'
                            DISCHARGE_CCES
                            iL=iL+1;
                    end
                end
                
                if Loffdesign && iix==1
                    SET_DESIGN
                end
                
            end
            
            % Compute energy balance
            %ENERGY_BALANCE
            ENERGY_BALANCE_v2
            
            end %end while loop for simple_interface
            
            % Evaluate the system cost
            PTES_ECONOMICS
            
            % Save results from multi_run call
            if multi_run, PRINT_MULTI_RUN; end
            
            % If using the simple interface, set these variables
            if simple_interface==1, PRINT_SIMPLE_INTERFACE; end
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if Nobjs > 1
    CONSTRAINT_AND_OUTPUT
end

catch ME
    if Nobjs > 1
        OPT_ERROR
    else
        if simple_interface==1
            dbclear all
        end
        rethrow(ME);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FINISH PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close files, save plots and release CoolProp AbstractStates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
