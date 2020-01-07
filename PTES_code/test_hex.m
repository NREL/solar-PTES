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
        addpath('./Classes/','./Generic/','./Functions/','./Other')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','./Other')
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
% 4 = Steam and MEG
% 5 = sCO2 and sCO2
scenario = 5;

% Set indices
iL = 1; i1 = 1; i2 = 1;

% Declare fluids and specify initial conditions
switch scenario
    case 1
        % Helium (hot, high pressure)
        F1 = fluid_class('Helium','WF','CP','TTSE',1,5);
        F1.state(iL,i1).p = 100e5;
        F1.state(iL,i1).T = 600;
        F1.state(iL,i1).mdot = 10;
        
        % Helium (cold, low pressure)
        F2 = fluid_class('Helium','WF','CP','TTSE',1,5);
        F2.state(iL,i2).p = 20e5;
        F2.state(iL,i2).T = 300;
        F2.state(iL,i2).mdot = 10;
        
        % Set hex_mode
        hex_mode = 0; % Both mass flow rates specified
        par = 0;
        
    case 2        
        % Solar Salt
        F1 = fluid_class('SolarSalt','SF','TAB',NaN,1,5);
        F1.state(iL,i1).p = 1e5;
        F1.state(iL,i1).T = 800;
        F1.state(iL,i2).mdot = 45;
        
        % Water
        F2 = fluid_class('Water','WF','CP','HEOS',1,5);
        F2.state(iL,i2).p = 100*1e5;
        F2.state(iL,i2).T = 400;
        F2.state(iL,i2).mdot = 10;
        
        % Set hex_mode
        hex_mode = 0; % Compute mass flow rate of hot fluid with var=mH*CpH/(mC*CpC)
        par = 1.10;
        %hex_mode = 4;
        %par = 300+273;
        
    case 3        
        % CO2
        F1 = fluid_class('CarbonDioxide','WF','CP','TTSE',1,5); % working fluid
        F1.state(1,i1).p = 85e5;
        F1.state(1,i1).T = 380;
        F1.state(1,i1).mdot = 0.75;
        
        % Water
        F2 = fluid_class('Water','SF','CP','HEOS',1,5); % storage fluid       
        F2.state(1,i2).p = 5e5;
        F2.state(1,i2).T = 300;
        F2.state(1,i2).mdot = 1;
        
        % Set hex_mode
        hex_mode = 0; % Both mass flow rates specified
        par = 0;
        
    case 4        
        % MEG
        F1 = fluid_class('INCOMP::MEG2[0.56]','SF','TAB',NaN,1,5);
        F1.state(iL,i1).p = 1e5;
        F1.state(iL,i1).T = 245;
        F1.state(iL,i2).mdot = 200;
        
        % Water
        F2 = fluid_class('Water','WF','CP','HEOS',1,5);
        F2.state(iL,i2).p = 0.1*1e5;
        F2.state(iL,i2).T = 322;
        F2.state(iL,i2).mdot = 10;
        
        % Set hex_mode
        hex_mode = 5;
        par = 315;
        %hex_mode = 0;
        %par = 0;
        
    case 5        
        % CO2
        F1 = fluid_class('CarbonDioxide','WF','CP','TTSE',1,5); % working fluid
        F1.state(1,i1).p = 30.86e5;
        F1.state(1,i1).T = 776.9 + 273.15;
        F1.state(1,i1).mdot = 58.0;
        
        % Water
        F2 = fluid_class('CarbonDioxide','WF','CP','TTSE',1,5); % storage fluid       
        F2.state(1,i2).p = 297.62e5;
        F2.state(1,i2).T = 81.9 + 273.15;
        F2.state(1,i2).mdot = 58.0;
        
        % Set hex_mode
        hex_mode = 0; % Both mass flow rates specified
        par = 0;
        
    otherwise
        error('not implemented')
end
% Update fluid states
[F1] = update(F1,[iL,i1],1);
[F2] = update(F2,[iL,i2],1);

% Specify HX settings
NX = 100; % Number of sections (grid)
model = 'geom'; % either 'eff', 'UA' or 'geom'
stage_type = 'hex';
method = 'automatic';
HX = hx_class('hot', stage_type, model, 1.00, 0.00,  4, NX, 2, 2);

switch HX.model
    case 'eff'
        HX.eff   = 1.00;
        HX.ploss = 0.01;
        
    case 'UA'
        HX.UA    = 1e6;
        HX.ploss = 0.01;        
        
    case 'geom'
        % Specify HEX geometry
        switch method
            case 'manual'
                % Define heat exchanger geometry (shell-and-tube)
                % 1 refers to the tube side, 2 refers to the shell side
                HX.shape     = 'circular';
                HX.sigma     = 1e8;        % Maximum allowable stress, Pa
                HX.t1_min    = 0.1e-3;     % Minimum tube thickness, m
                HX.t1_D1_min = 0.05;       % Minimum tube thickness-to-diameter ratio
                HX.L         = 3.0;        % Tube length, m
                HX.D1        = 5e-3;       % Tube diameter, m
                HX.AfT       = 0.25;       % Total flow area, m2
                HX.AfR       = 1.00;       % Ratio of flow areas, Af2/Af1, -
                
            case 'automatic'
                % Obtain geometric parameters based on performance objectives,
                % using analytical solutions.
                NTU   = 4.5;
                ploss = 0.026;
                D     = 0.95e-3;
                [HX]  = set_hex_geom(HX,iL,F1,i1,F2,i2,hex_mode,par,NTU,ploss,D);
        end
        
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% COMPUTE AND MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[HX,~,~,~,~] = hex_func(HX,iL,F1,i1,F2,i2,hex_mode,par);

% Make plots
plot_hex(HX,1,20,'C');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare specifications from set_hex_geom with numerical results
if all([strcmp(HX.model,'geom'),strcmp(method,'automatic')])
    
    fprintf(1,'\n      Specification  Result\n')
    fprintf(1,'NTU_min = %8.3f   %9.3f\n',NTU,HX.NTU)
    fprintf(1,'DppH    = %8.5f   %9.5f\n',ploss,HX.DppH)
    fprintf(1,'DppC    = %8.5f   %9.5f\n',ploss,HX.DppC)
    
end

%%% COMPARE GEOMETRICAL MODEL WITH ANALYTICAL SOLUTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 30;
if n>1
    mdot1 = F1.state(iL,i1).mdot*linspace(0.25,1.0,n);
    mdot2 = F2.state(iL,i1).mdot*linspace(0.25,1.0,n);
    %mdot2 = F2.state(iL,i1).mdot*ones(size(mdot1));
else
    mdot1 = F1.state(iL,i1).mdot;
    mdot2 = F2.state(iL,i1).mdot;
end

AN_TH1  = zeros(size(mdot1));
AN_TC2  = zeros(size(mdot1));
AN_pH1  = zeros(size(mdot1));
AN_pC2  = zeros(size(mdot1));
AN_DppH = zeros(size(mdot1));
AN_DppC = zeros(size(mdot1));
NU_TH1  = zeros(size(mdot1));
NU_TC2  = zeros(size(mdot1));
NU_pH1  = zeros(size(mdot1));
NU_pC2  = zeros(size(mdot1));
NU_DppH = zeros(size(mdot1));
NU_DppC = zeros(size(mdot1));
for im = 1:n
    F1.state(iL,i1).mdot = mdot1(im);
    F2.state(iL,i1).mdot = mdot2(im);
    % Run HEX code
    [HX,F1,~,F2,~] = hex_func(HX,iL,F1,i1,F2,i2,hex_mode,par);
    
    % Compare numerical results (NU) with analytical results (AN)
    NU_TH1(im)  = HX.H.T(1);
    NU_TC2(im)  = HX.C.T(NX+1);
    NU_pH1(im)  = HX.H.p(1);
    NU_pC2(im)  = HX.C.p(HX.NX+1);
    NU_DppH(im) = (HX.H.pin-HX.H.p(1))/HX.H.pin;
    NU_DppC(im) = (HX.C.pin-HX.C.p(HX.NX+1))/HX.C.pin;
    [AN_TH1(im),AN_TC2(im),AN_pH1(im),AN_pC2(im)] = hex_analytic(HX,iL,F1,i1,F2,i2);
    AN_DppH(im) = (HX.H.pin-AN_pH1(im))/HX.H.pin;
    AN_DppC(im) = (HX.C.pin-AN_pC2(im))/HX.C.pin;
%     fprintf(1,'\n      Numerical  Analytical\n')
%     fprintf(1,'TH1 [C] = %8.1f   %8.1f\n',NU_TH1(im)-273.15,AN_TH1(im)-273.15)
%     fprintf(1,'TC2 [C] = %8.1f   %8.1f\n',NU_TC2(im)-273.15,AN_TC2(im)-273.15)
%     fprintf(1,'DppH    = %8.5f   %8.5f\n',NU_DppH(im),AN_DppH(im))
%     fprintf(1,'DppC    = %8.5f   %8.5f\n',NU_DppC(im),AN_DppC(im))
    
    % Make plots
    %plot_hex(HX,1,20,'C');
    %pause(2);
end
% if n>1
%     figure(25)
%     plot(mdot1,NU_TH1,'o',mdot1,AN_TH1)
%     xlabel('Mass flow rate, kg/s')
%     ylabel('Outlet temperature, K')
%     legend('numerical','analytical')
%     figure(26)
%     plot(mdot1,NU_TC2,'o',mdot1,AN_TC2)
%     xlabel('Mass flow rate, kg/s')
%     ylabel('Outlet temperature, K')
%     legend('numerical','analytical')
%     figure(27)
%     semilogy(mdot1,NU_DppH,'o',mdot1,AN_pH1,mdot1,NU_DppC,'o',mdot1,AN_pC2)
%     xlabel('Mass flow rate, kg/s')
%     ylabel('$ \Delta p / p $')
%     legend('FH numerical','FH analytical','FC numerical','FC analytical','Location','Best')
% end

% Load data from Hoopes2016
data = load('Validation.csv');

% Compare numerical results with data from Hoopes2016
figure(30)
%plot(mdot1,NU_TH1-273.15,data(:,1),data(:,7),'s')
plot(mdot1,NU_TH1-273.15,data(:,1),data(:,3),'s')
xlabel('Mass flow rate [kg/s]')
ylabel('Outlet temperature (hot) [K]')
legend('numerical','Hoopes2016','Location','Best')
figure(31)
%plot(mdot1,NU_TC2-273.15,data(:,1),data(:,8),'s')
plot(mdot1,NU_TC2-273.15,data(:,1),data(:,4),'s')
xlabel('Mass flow rate [kg/s]')
ylabel('Outlet temperature (cold) [K]')
legend('numerical','Hoopes2016','Location','Best')
figure(32)
%plot(mdot1,NU_pH1/1e5,data(:,1),data(:,9),'s')
plot(mdot1,NU_pH1/1e5,data(:,1),data(:,5),'s')
xlabel('Mass flow rate [kg/s]')
ylabel('Outlet pressure (hot) [bar]')
legend('numerical','Hoopes2016','Location','Best')
figure(33)
%plot(mdot1,NU_pC2/1e5,data(:,1),data(:,10),'s')
plot(mdot1,NU_pC2/1e5,data(:,1),data(:,6),'s')
xlabel('Mass flow rate [kg/s]')
ylabel('Outlet pressure (cold) [bar]')
legend('numerical','Hoopes2016','Location','Best')
% figure(34)
% plot(mdot1,NU_DppH,data(:,1),(HX.H.pin-data(:,5)*1e5)/HX.H.pin,'s')
% xlabel('Mass flow rate, kg/s')
% ylabel('$ \Delta p / p $ (hot)')
% legend('numerical','Hoopes2016','Location','Best')
% figure(35)
% plot(mdot1,NU_DppC,data(:,1),(HX.C.pin-data(:,6)*1e5)/HX.C.pin,'s')
% xlabel('Mass flow rate, kg/s')
% ylabel('$ \Delta p / p $ (hot)')
% legend('numerical','Hoopes2016','Location','Best')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%