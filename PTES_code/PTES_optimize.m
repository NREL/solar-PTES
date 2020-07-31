function [fit err extra]=PTES_optimize(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTES
% This code employs thermodynamic and economic models to predict the
% performance and cost of different pumped thermal energy storage cycles.
%
% Authors: Pau Farres-Antunez and Joshua Dominic McTigue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enter debugging mode if an error occurs
dbstop if error
%dbclear all

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/','./Other/','./LIB/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\','.\Other\','.\LIB\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop

%Optimization inputs
TH_dis0 = x(1);
Ne_ch   = round (x(2));
eff = x(3);

% SET INPUTS
INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic % start timer
try
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
                            JB_DISCHARGE
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
            
            % Evaluate the system cost
            PTES_ECONOMICS
            
            % Save results from multi_run call
            if multi_run, PRINT_MULTI_RUN; end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
toc %stop timer
err= zeros(1,1);
%if fluidH.state(2,3).T <TH_dis0
 %  f1=0.8+(0.9-0.8)*rand(1);
  % f2=0.8+(0.9-0.8)*rand(1);
  % extra=1000000000;
    
%else
    f1=1-chi_PTES_para;
    f2=Cdata.lcosM;
    extra=Cdata.cap_costM;
    
%end
fit=[f1 f2];
catch
    f1= 0.8+(0.9-0.8)*rand(1);
    f2= 0.8+(0.9-0.8)*rand(1);
    extra=1000000000;
    err= zeros(1,1);
    fit=[f1 f2];
    
end
end

