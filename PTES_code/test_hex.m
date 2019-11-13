%%% To run this script, place it inside the PTES_code folder, alongside the
%%% PTES.m file

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
        addpath('./Classes/','./Generic/','./Functions/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% INPUTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose scenario
% 1 = Helium and Helium
% 2 = SolarSalt and Water
% 3 = CO2 and Water
scenario = 1;

% Set indices
iL = 1; i1 = 1; i2 = 1;

% Declare fluids and specify initial conditions
switch scenario
    case 1
        % Helium (hot, high pressure)
        F1 = fluid_class('Argon','WF','CP','TTSE',1,5);
        F1.state(iL,i1).p = 100e5;
        F1.state(iL,i1).T = 600;
        F1.state(iL,i1).mdot = 10;
        
        % Helium (cold, low pressure)
        F2 = fluid_class('Argon','WF','CP','TTSE',1,5);
        F2.state(iL,i2).p = 20e5;
        F2.state(iL,i2).T = 300;
        F2.state(iL,i2).mdot = 10;
        
        % Set hex_mode
        hex_mode = 0; % Both mass flow rates specified
        var = 0;
        
    case 2        
        % Solar Salt
        F1 = fluid_class('SolarSalt','SF','TAB',NaN,1,5);
        F1.state(iL,i1).p = 1e5;
        F1.state(iL,i1).T = 830;
        
        % Water
        F2 = fluid_class('Water','WF','CP','TTSE',1,5);
        F2.state(iL,i2).p = 100*1e5;
        F2.state(iL,i2).T = 400;
        F2.state(iL,i2).mdot = 10;
        
        % Set hex_mode
        hex_mode = 2; % Compute mass flow rate of hot fluid with var=mH*CpH/(mC*CpC)
        var = 1.20;     
        
    case 3        
        % CO2
        F1 = fluid_class('CarbonDioxide','WF','CP','TTSE',1,5); % working fluid
        F1.state(1,i1).p = 85e5;
        F1.state(1,i1).T = 380;
        F1.state(1,i1).mdot = 0.75;
        
        % Water
        F2 = fluid_class('Water','SF','CP','TTSE',1,5); % storage fluid       
        F2.state(1,i2).p = 5e5;
        F2.state(1,i2).T = 300;
        F2.state(1,i2).mdot = 1;
        
        % Set hex_mode
        hex_mode = 0; % Both mass flow rates specified
        var = 0;
        
    otherwise
        error('not implemented')
end

% Specify HEX geometry
method = 'manual';
switch method
    case 'manual'
        % Define heat exchanger geometry (shell-and-tube)
        % 1 refers to the tube side, 2 refers to the shell side
        HX.shape     = 'circular';
        HX.sigma     = 1e8;        % Maximum allowable stress, Pa
        HX.L         = 3.0;        % Tube length, m
        HX.D1        = 0.5e-2;     % Tube diameter, m
        HX.t1        = 0.1*HX.D1;  % Tube thickness, m
        HX.AfT       = 0.25;        % Total flow area, m2
        HX.AfR       = 1.00;       % Ratio of flow areas, Af2/Af1, -
        
    case 'automatic'
        % Obtain geometric parameters based on performance objectives,
        % using analytical solutions. It should be expected that objectives
        % will be met accurately (only) when using fluids with small
        % variations of thermophysical properties.
        eff   = 0.97;
        ploss = 0.01;
        warning('under development')
        [HX]  = set_hex_geom(F1,[iL,i1],F2,[iL,i2],eff,ploss,5e-3,0.5e-3,1e8,hex_mode,var);
end

% Hex code settings
HX.NX = 100;               % Number of sections (grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% COMPUTE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update fluid states
[F1] = update(F1,[iL,i1],1);
[F2] = update(F2,[iL,i2],1);

mdot = F1.state(iL,i1).mdot*linspace(0.6,1.4,20);
NU_eff = zeros(size(mdot));
AN_eff = zeros(size(mdot));
for im = 1:length(mdot)
    F1.state(iL,i1).mdot = mdot(im);
    % Run HEX code
    [F1,F2,~,~,HX] = hex_TQA(F1,[iL,i1],F2,[iL,i2],HX,'hex',hex_mode,var);
    
    % Compare numerical results (NU) with analytical results (AN)
    NU_eff(im)  = HX.QS(HX.NX+1)/HX.QMAX;
    NU_DppH = (HX.H.pin-HX.H.p(1))./HX.H.pin;
    NU_DppC = (HX.C.pin-HX.C.p(HX.NX+1))./HX.C.pin;
    [AN_eff(im),AN_DppH,AN_DppC] = hex_analytic(F1,[iL,i1],F2,[iL,i2],HX);
    fprintf(1,'\n      Numerical  Analytical\n')
    fprintf(1,'Eff  = %8.3f   %9.3f\n',NU_eff(im),AN_eff(im))
    fprintf(1,'DppH = %8.5f   %9.5f\n',NU_DppH,AN_DppH)
    fprintf(1,'DppC = %8.5f   %9.5f\n',NU_DppC,AN_DppC)
end
figure(25)
plot(mdot,NU_eff,mdot,AN_eff)
xlabel('mass flow rate')
ylabel('effectiveness')
legend('numerical','analytical')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots
plot_hex_TQA(HX,20,'C');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Save figures
% switch scenario
%     case 1
%         save_fig(1,'./Results/TQ_H2O_CO2_chg',0,0,0)
%         save_fig(2,'./Results/TQ_H2O_CO2_dis',0,0,0)
% end