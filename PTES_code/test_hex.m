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
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/','./Other/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\','.\Other\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% INPUTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose scenario
% 1 = Gas regenerator
% 2 = SolarSalt and Water
% 3 = CO2 and Water
% 4 = Steam and Water
% 5 = sCO2 and sCO2
% 6 = Heat rejection unit (Nigroten and Nitrogen)
scenario = 4;

% Save figures?
save_figures = 0;

% Set indices
iL = 1; i1 = 1; i2 = 1;

% Declare fluids and specify initial conditions
switch scenario
    case 1
        % Helium (hot, high pressure)
        F1 = fluid_class('Nitrogen','WF','CP','TTSE',1,5);
        F1.state(iL,i1).p = 100e5;
        F1.state(iL,i1).T = 600;
        F1.state(iL,i1).mdot = 10;
        
        % Helium (cold, low pressure)
        F2 = fluid_class('Nitrogen','WF','CP','TTSE',1,5);
        F2.state(iL,i2).p = 20e5;
        F2.state(iL,i2).T = 300;
        F2.state(iL,i2).mdot = 10;
        
        % Set hex_mode
        hex_mode = 0;
        par = 0;
        
    case 2        
        % Solar Salt
        F1 = fluid_class('SolarSalt','SF','TAB',NaN,1,5);
        F1.state(iL,i1).p = 1e5;
        F1.state(iL,i1).T = 565+273.15;
        F1.state(iL,i2).mdot = 600;
        
        % Water
        F2 = fluid_class('Water','WF','CP','TTSE',1,5);
        F2.state(iL,i2).p = 100*1e5;
        F2.state(iL,i2).T = 180+273;
        F2.state(iL,i2).mdot = 100;
        
        % Set hex_mode and stage_type
        hex_mode = 0;
        par = 1.00;
        %hex_mode = 4;
        %par = 300+273;
        %hex_mode = 3;
        %par = 450+273;
        
    case 3        
        % CO2
        F1 = fluid_class('CarbonDioxide','WF','CP','TTSE',1,5);
        F1.state(1,i1).p = 85e5;
        F1.state(1,i1).T = 380;
        F1.state(1,i1).mdot = 0.75;
        
        % Water
        F2 = fluid_class('Water','SF','CP','TTSE',1,5);
        F2.state(1,i2).p = 5e5;
        F2.state(1,i2).T = 300;
        F2.state(1,i2).mdot = 1;
        
        % Set hex_mode
        hex_mode = 0;
        par = 0;
        
    case 4
        % Steam
        F1 = fluid_class('Water','WF','CP','TTSE',1,5);
        F1.state(iL,i2).p = 0.1*1e5;
        F1.state(iL,i2).T = 350;
        F1.state(iL,i2).mdot = 10;
        
        % Water
        F2 = fluid_class('Water','SF','TAB',NaN,1,5);
        F2.state(iL,i1).p = 1e5;
        F2.state(iL,i1).T = 273.15+1;
        F2.state(iL,i2).mdot = 200;
        %{
        % MEG
        F2 = fluid_class('INCOMP::MEG2[0.56]','SF','TAB',NaN,1,5);
        F2.state(iL,i1).p = 1e5;
        F2.state(iL,i1).T = 245;
        F2.state(iL,i2).mdot = 200;
        %}
        
        % Set hex_mode
        hex_mode = 5;
        par = 300;
        %hex_mode = 0;
        %par = 0;
        
    case 5        
        % CO2
        F1 = fluid_class('CarbonDioxide','WF','CP','TTSE',1,5);
        F1.state(1,i1).p = 30.86e5;
        F1.state(1,i1).T = 776.9 + 273.15;
        F1.state(1,i1).mdot = 58.0;
        
        % CO2
        F2 = fluid_class('CarbonDioxide','WF','CP','TTSE',1,5); % storage fluid       
        F2.state(1,i2).p = 297.62e5;
        F2.state(1,i2).T = 81.9 + 273.15;
        F2.state(1,i2).mdot = 58.0;
        
        % Set hex_mode
        hex_mode = 0;
        par = 0;
        
    case 6
        % Nitrogen (WF)
        F1 = fluid_class('Nitrogen','WF','CP','TTSE',1,5);
        F1.state(1,i1).p = 10e5;
        F1.state(1,i1).T = 100 + 273.15;
        F1.state(1,i1).mdot = 1e3;
        
        % Nitrogen (ENV)
        F2 = fluid_class('Nitrogen','ENV','CP','TTSE',1,5); % storage fluid       
        F2.state(1,i2).p = 1e5;
        F2.state(1,i2).T = 20 + 273.15;
        F2.state(1,i2).mdot = 0;
        
        % Set hex_mode
        hex_mode = 5;
        par = F2.state(1,i2).T + 5;
        
    otherwise
        error('not implemented')
end
% Update fluid states
[F1] = update(F1,[iL,i1],1);
[F2] = update(F2,[iL,i2],1);

% Specify HX settings
NX   = 100; % Number of sections (grid)
name = 'hot';
stage_type = 'hex';
model = 'geom';

switch model
    case 'eff'
        eff   = 0.97;
        ploss = 0.01;
        par1  = eff;
        par2  = ploss;
        par3  = [];
        par4  = [];
        
    case 'UA'
        UA    = 1e6;
        ploss = 0.01;
        par1  = UA;
        par2  = ploss;
        par3  = [];
        par4  = [];
        
    case 'geom'
        % Specify HEX geometry based on performance objectives
        switch scenario
            case 1
                % For comparisson with analytical results
                eff   = 0.95;
                ploss = 0.01;
                D1    = 2e-2;
            case 5
                % For comparisson with scaling model (Hoopes2016)
                %NTU   = 4.45;
                %ploss = 0.03;
                %D     = 2e-2;
                % For comparisson with numerical model (Hoopes2016)
                eff   = 0.95;
                ploss = 0.025;
                D1    = 1.00e-3;
            otherwise
                eff   = 0.97;
                ploss = 0.01;
                D1    = 0.01;                
        end
        
        par1  = eff;
        par2  = ploss;
        par3  = D1;
        par4  = 'circular';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CONSTRUCT HEAT EXCHANGER %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HX = hx_class(name, stage_type, 4, NX, 1, 1, model, par1, par2, par3, par4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% DESIGN PERFORMANCE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run heat exchanger model under design conditions
[HX,~,~,~,~] = hex_func(HX,iL,F1,i1,F2,i2,hex_mode,par);

%% SUMMARY
% Compare specifications from hex_set_geom with numerical results
if strcmp(HX.model,'geom')    
    fprintf(1,'\n      Specification  Result\n')
    fprintf(1,'Eff     = %8.3f   %9.3f\n',eff,1-min(HX.H(1).T-HX.C(1).T)/(HX.H(1).T(end)-HX.C(1).T(1)))
    fprintf(1,'DppH    = %8.5f   %9.5f\n',ploss,HX.DppH)
    fprintf(1,'DppC    = %8.5f   %9.5f\n',ploss,HX.DppC)
    
    print_hexs(HX,1,'Summary:\n');
end

% Make plots
plot_hex(HX,1,10,'C');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%% OFF-DESIGN PERFORMANCE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create mass flow rate arrays
switch scenario
    case 1
        n = 20;
        mdot1 = F1.state(iL,i1).mdot*linspace(0.25,1.0,n);
        mdot2 = F2.state(iL,i1).mdot*linspace(0.25,1.0,n);
        %mdot2 = F2.state(iL,i1).mdot*ones(size(mdot1));
    case 5
        % Load data from Hoopes2016
        data  = load('Validation.csv');
        mdot1 = data(:,1)';
        mdot2 = mdot1;
    otherwise
        mdot1 = F1.state(iL,i1).mdot;
        mdot2 = F2.state(iL,i1).mdot;
end
n = length(mdot1);

% Use easier nomenclature for inlet conditions
TC1 = HX.C.T(1);
TH2 = HX.H.T(HX.NX+1);
pC1 = HX.C.pin;
pH2 = HX.H.pin;

% COMPUTE OFF-DESIGN PERFORMANCE
% Allocate arrays (numerical results)
num_TH1  = zeros(size(mdot1));
num_TC2  = zeros(size(mdot1));
num_pH1  = zeros(size(mdot1));
num_pC2  = zeros(size(mdot1));
num_DppH = zeros(size(mdot1));
num_DppC = zeros(size(mdot1));
% Allocate arrays (analytical results)
an_TH1  = zeros(size(mdot1));
an_TC2  = zeros(size(mdot1));
an_pH1  = zeros(size(mdot1));
an_pC2  = zeros(size(mdot1));
an_DppH = zeros(size(mdot1));
an_DppC = zeros(size(mdot1));
for im = 1:n
    % Select mass flow rate
    F1.state(iL,i1).mdot = mdot1(im);
    F2.state(iL,i1).mdot = mdot2(im);
    
    % Run HEX code
    [HX,F1,~,F2,~] = hex_func(HX,iL,F1,i1,F2,i2,hex_mode,par);
    
    % Extract data for plotting
    num_TH1(im)  = HX.H.T(1);
    num_TC2(im)  = HX.C.T(NX+1);
    num_pH1(im)  = HX.H.p(1);
    num_pC2(im)  = HX.C.p(HX.NX+1);
    num_DppH(im) = (HX.H.pin-HX.H.p(1))/HX.H.pin;
    num_DppC(im) = (HX.C.pin-HX.C.p(HX.NX+1))/HX.C.pin;
    
    switch scenario
        case 1
            % Generate analytical results
            [an_TH1(im),an_TC2(im),an_pH1(im),an_pC2(im)] = hex_analytic(HX,iL,F1,i1,F2,i2);
            an_DppH(im) = (HX.H.pin-an_pH1(im))/HX.H.pin;
            an_DppC(im) = (HX.C.pin-an_pC2(im))/HX.C.pin;
            
            % Print comparison
            fprintf(1,'\n      Numerical  Analytical\n')
            fprintf(1,'TH1 [C] = %8.1f   %8.1f\n',num_TH1(im)-273.15,an_TH1(im)-273.15)
            fprintf(1,'TC2 [C] = %8.1f   %8.1f\n',num_TC2(im)-273.15,an_TC2(im)-273.15)
            fprintf(1,'DppH    = %8.5f   %8.5f\n',num_DppH(im),an_DppH(im))
            fprintf(1,'DppC    = %8.5f   %8.5f\n',num_DppC(im),an_DppC(im))
    end
    
    % Make plots
    %plot_hex(HX,1,20,'C');
    %pause(2);
end

% COMPARE WITH ANALYTICAL SOLUTIONS
switch scenario
    case 1
        % Plot figures
        figure(25)
        plot(mdot1,num_TH1,'ro',mdot1,an_TH1,'r')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet temperature (hot) [K]')
        legend('numerical','analytical')
        figure(26)
        plot(mdot1,num_TC2,'bo',mdot1,an_TC2,'b')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet temperature (cold) [K]')
        legend('numerical','analytical')
        figure(27)
        semilogy(mdot1,num_DppH,'ro',mdot1,an_DppH,'r',mdot1,num_DppC,'bo',mdot1,an_DppC,'b')
        xlabel('Mass flow rate [kg/s]')
        ylabel('$ \Delta p / p $')
        legend('numerical (hot)','analytical (hot)','numerical (cold)','analytical (cold)','Location','Best')
        
        % Compute errors
        DTmax   = TH2 - TC1;
        errTH1  = max(abs(num_TH1 - an_TH1)/DTmax*100);
        errTC2  = max(abs(num_TC2 - an_TC2)/DTmax*100);
        errpH1  = max(abs(num_pH1 - an_pH1)./num_pH1*100);
        errpC2  = max(abs(num_pC2 - an_pC2)./num_pC2*100);
        errDppH = max(abs(num_DppH - an_DppH)./num_DppH*100);
        errDppC = max(abs(num_DppC - an_DppC)./num_DppC*100);
end


% COMPARE WITH DATA OR MODELS FROM LITERATURE
switch scenario
    case 5 % Compare numerical results with data from Hoopes2016
        % Compute errors
        DTmax   = TH2 - TC1;
        errTH1  = abs(num_TH1-273.15 - data(:,3)')/DTmax*100;
        errTC2  = abs(num_TC2-273.15 - data(:,4)')/DTmax*100;
        errDppH = abs(num_pH1 - data(:,5)'*1e5)./pH2*100;
        errDppC = abs(num_pC2 - data(:,6)'*1e5)./pC1*100;
        
        % Plot figures
        figure(30)
        yyaxis left
        plot(mdot1,num_TH1-273.15,data(:,1),data(:,3),'s');
        ylim([85 105])
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet temperature (hot) [K]')        
        yyaxis right
        plot(mdot1,errTH1,'d')
        ylim([0 0.2])
        ylabel('Error [$\%$]')
        legend('Numerical','Hoopes2016','Error','Location','NorthWest')        
        
        figure(31)
        yyaxis left
        plot(mdot1,num_TC2-273.15,data(:,1),data(:,4),'s')
        ylim([615 635])
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet temperature (cold) [K]')
        yyaxis right
        plot(mdot1,errTC2,'d')
        ylim([0 0.10])
        ylabel('Error [$\%$]')
        legend('Numerical','Hoopes2016','Error','Location','Best')
        
        figure(32)
        yyaxis left
        plot(mdot1,num_pH1/1e5,data(:,1),data(:,5),'s')
        ylim([29 31])
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet pressure (hot) [bar]')
        yyaxis right
        plot(mdot1,errDppH,'d')
        ylim([0 0.10])
        ylabel('Error [$\%$]')
        legend('Numerical','Hoopes2016','Error','Location','Best')
        
        figure(33)
        yyaxis left
        plot(mdot1,num_pC2/1e5,data(:,1),data(:,6),'s')
        ylim([290 300])
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet pressure (cold) [bar]')
        yyaxis right
        plot(mdot1,errDppC,'d')
        ylim([0 2.5])
        ylabel('Error [$\%$]')
        legend('Numerical','Hoopes2016','Error','Location','West')
        
        % Compare against Hoopes' scaling method
        %figure(30)
        %plot(mdot1,num_TH1-273.15,data(:,1),data(:,7),'s');
        %figure(31)
        %plot(mdot1,num_TC2-273.15,data(:,1),data(:,8),'s')
        %figure(32)
        %plot(mdot1,num_pH1/1e5,data(:,1),data(:,9),'s')
        %figure(33)
        %plot(mdot1,num_pC2/1e5,data(:,1),data(:,10),'s')
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% SAVE FIGURES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch save_figures
    case 1
        formats = {'fig','epsc','emf'};
        
        save_fig(10,'./Results/T_Q',formats)
        save_fig(11,'./Results/T_A',formats)
        save_fig(12,'./Results/p_A',formats)
        save_fig(13,'./Results/Re_A',formats)
        save_fig(14,'./Results/ht_Q',formats)
        
        switch scenario
            case 1
                save_fig(25,'./Results/TH1_analytic',formats)
                save_fig(26,'./Results/TC2_analytic',formats)
                save_fig(27,'./Results/Dpp_analytic',formats)
            case 5
                save_fig(30,'./Results/TH1',formats)
                save_fig(31,'./Results/TC2',formats)
                save_fig(32,'./Results/pH1',formats)
                save_fig(33,'./Results/pC2',formats)
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%}