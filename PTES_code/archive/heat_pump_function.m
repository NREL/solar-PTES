%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heat_pump_function
% This a thermodynamic model of a pumped thermal energy storage cycle. It
% is based on the Joule-Brayton cycle and employs liquid storage media.
% Author: Pau Farres-Antunez (pf298@cam.ac.uk)
%
% This function is based on PTES.m (see above).
% It has been turned into a function so that it can be called and iterated
% upon, with different input variables. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Have to ensure diff is first output variable so can use function in
% golden search routine
function [diff, out1, out2, out3] = heat_pump_function(MShotT, MScoldT, cmpTIN, T0actual, plts)
%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_graphics
% SET INPUTS
PTES_INPUTS
TH_dis0    = MScoldT ; % Rewrite this variable (the molten salt cold temp)
Tmax       = cmpTIN ;
make_plots = plts;
T0         = T0actual ;

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
        PTES_INITIALISE
        
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
        if make_plots
            PTES_ENERGY_BALANCE
        end
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
% Close files, save plots and release CoolProp AbstractStates
PTES_FINISH

% Variables to return
diff = abs(fluidH.state(1,2).T - MShotT) ;
out1 = fluidH ;
out2 = fluidC ;
out3 = gas ;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%