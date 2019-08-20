% Clear Workspace variables
clear;

% Set paths
addpath('./Classes/')
addpath('./Generic/')
addpath('./Functions/')

% Set properties for plots
set_graphics

% Load CoolProp library
load_coolprop


%%%% INPUTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose scenario
scenario = 1; % 1 = CO2 and water, 2 = CO2 and mineral oil

% Specify HEX performance
eff   = 0.97;
ploss = 0.01;

% Set fluid objects and specify initial conditions
switch scenario
    case 1
%         % Water
%         F1 = fluid_class('Water','SF','CP',30); % storage fluid       
%         F1.state(1,1).p = 1e5;
%         F1.state(1,1).T = 273.15 + 95;
%         F1.state(1,1).mdot = 0.9;
%         
%         % CO2
%         F2 = fluid_class('CarbonDioxide','WF','CP',30); % working fluid
%         F2.state(1,1).p = 79.2e5;
%         F2.state(1,1).T = 273.15 + 16.4;
%         F2.state(1,1).mdot = 1;
        
        % Water
        F1 = fluid_class('Water','SF','CP',30); % storage fluid       
        F1.state(1,1).p = 5e5;
        F1.state(1,1).T = 273.15 + 31.9;
        F1.state(1,1).mdot = 1;
        
        % CO2
        F2 = fluid_class('CarbonDioxide','WF','CP',30); % working fluid
        F2.state(1,1).p = 84.1e5;
        F2.state(1,1).T = 273.15 + 107.2;
        F2.state(1,1).mdot = 0.7319;
        
    case 2
        % Mineral oil
        F1 = fluid_class('MineralOil','SF','TAB',30); % storage fluid        
        F1.state(1,1).p = 1e5;
        F1.state(1,1).T = 273.15 + 355;
        F1.state(1,1).mdot = 0.65;
        [F1] = update(F1,[1,1],1);
        
        % CO2
        F2 = fluid_class('CarbonDioxide','WF','CP',30); % working fluid
        F2.state(1,1).p = 79.2e5;
        F2.state(1,1).T = 273.15 + 16.4;
        F2.state(1,1).mdot = 1;
        [F2] = update(F2,[1,1],1);
        
    case 3
        % Neon
        F1 = fluid_class('Neon','WF','CP',30); % working fluid        
        F1.state(1,1).p = 1e5;
        F1.state(1,1).T = 70;
        F1.state(1,1).mdot = 2.0;
        [F1] = update(F1,[1,1],1);
        
        % CO2
        F2 = fluid_class('Air','WF','CP',30); % working fluid
        F2.state(1,1).p = 100e5;
        F2.state(1,1).T = 273.15 + 15;
        F2.state(1,1).mdot = 1;
        [F2] = update(F2,[1,1],1);
    otherwise
        error('not implemented')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% RUN CHARGE PROCESS
% Update initial states
[F1] = update(F1,[1,1],1);
[F2] = update(F2,[1,1],1);
% Run HEX routine
[F1,F2,~,~] = hex_TQ_cond(F1,[1,1],F2,[1,1],eff,0,ploss,'hex',4,0);
% Plot T-Q diagram
%plot_hex(F2,[1,1],F1,[1,1],100,1)
plot_hex(F1,[1,1],F2,[1,1],100,1)


% RUN DISCHARGE PROCESS
% Set initial conditions
F1.state(2,1) = F1.state(1,2);
F2.state(2,1) = F2.state(1,2);
% Run HEX routine
[F2,F1,~,~] = hex_TQ_cond(F2,[2,1],F1,[2,1],eff,0,ploss,'hex',4,0);
% Plot T-Q diagram
plot_hex(F2,[2,1],F1,[2,1],100,2)


% Save figures
switch scenario
    case 1
        save_fig(1,'./Results/TQ_H2O_CO2_chg',0,0,0)
        save_fig(2,'./Results/TQ_H2O_CO2_dis',0,0,0)
    case 2
        save_fig(1,'./Results/TQ_Oil_CO2_chg',0,0,0)
        save_fig(2,'./Results/TQ_Oil_CO2_dis',0,0,0)
    case 3
        save_fig(1,'./Results/TQ_Air_Neon_chg',0,0,0)
        save_fig(2,'./Results/TQ_Air_Neon_dis',0,0,0)
end