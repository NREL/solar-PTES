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
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/','./Other/','./LIB/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\','.\Other\','.\LIB\')
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
% 4 = Steam condenser -- Steam and Water, counter-flow
% 5 = sCO2 and sCO2 (Hoopes2016)
% 6 = Heat rejection unit -- Nigroten and Nitrogen, Xflow
% 7 = sCO2 and sCO2 (Marchionni2019 - CFD)
% 8 = sCO2 and sCO2 (Marchionni2019 - experiment)
% 9 = Helium and Helium (Figley2013)
%10 = Water and air (Nellis example 8.1-1), Xflow
%11 = Steam condenser -- Steam and air, Xflow
scenario = 11;

% Save figures?
save_figures = 0;

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
        F1 = fluid_class('CarbonDioxide','WF','CP','BICUBIC&HEOS',1,5);
        F1.state(1,i1).p = 85e5;
        F1.state(1,i1).T = 380;
        F1.state(1,i1).mdot = 0.75;
        
        % Water
        F2 = fluid_class('Water','SF','CP','BICUBIC&HEOS',1,5);
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
        F1.state(iL,i2).mdot = 100;
        
        % Water
        F2 = fluid_class('Water','SF','TAB',NaN,1,5);
        F2.state(iL,i1).p = 1e5;
        F2.state(iL,i1).T = 273.15+1;
        F2.state(iL,i2).mdot = 2000;
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
        
    case 5 % Hoopes2016
        % CO2
        F1 = fluid_class('CarbonDioxide','WF','CP','BICUBIC&HEOS',1,5);
        F1.state(1,i1).p = 30.86e5;
        F1.state(1,i1).T = 776.9 + 273.15;
        F1.state(1,i1).mdot = 58.0;
        
        % CO2
        F2 = fluid_class('CarbonDioxide','WF','CP','BICUBIC&HEOS',1,5);
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
        F1.state(1,i1).T = 30 + 273.15;
        F1.state(1,i1).mdot = 1e3;
        
        % Nitrogen (ENV)
        F2 = fluid_class('Nitrogen','ENV','CP','TTSE',1,5); % storage fluid       
        F2.state(1,i2).p = 1e5;
        F2.state(1,i2).T = 15 + 273.15;
        F2.state(1,i2).mdot = 0;
        
        % Set hex_mode
        %hex_mode = 5;
        %par = F2.state(1,i2).T + 5;
        hex_mode = 1;
        par = 0.5;
        
    case 7 % Marchionni2019 - CFD validation
        % CO2
        F1 = fluid_class('CarbonDioxide','WF','CP','BICUBIC&HEOS',1,5);
        F1.state(1,i1).p = 75e5;
        F1.state(1,i1).T = 400 + 273.15;
        F1.state(1,i1).mdot = 1.6;
        
        % CO2
        F2 = fluid_class('CarbonDioxide','WF','CP','BICUBIC&HEOS',1,5);
        F2.state(1,i2).p = 150e5;
        F2.state(1,i2).T = 100 + 273.15;
        F2.state(1,i2).mdot = 1.6;
        
        % Set hex_mode
        hex_mode = 0;
        par = 0;
        
    case 8 % Marchionni2019 - experimental data
        % CO2
        F1 = fluid_class('CarbonDioxide','WF','CP','BICUBIC&HEOS',1,5);
        F1.state(1,i1).p = 75e5;
        F1.state(1,i1).T = 344.3 + 273.15;
        F1.state(1,i1).mdot = 2.62;
        
        % CO2
        F2 = fluid_class('CarbonDioxide','WF','CP','BICUBIC&HEOS',1,5);
        F2.state(1,i2).p = 125e5;
        F2.state(1,i2).T = 72.9 + 273.15;
        F2.state(1,i2).mdot = 2.62;
        
        % Set hex_mode
        hex_mode = 0;
        par = 0;
        
    case 9 % Figley2013
        % Helium
        F1 = fluid_class('Helium','WF','CP','BICUBIC&HEOS',1,5);
        F1.state(1,i1).p = 30e5;
        F1.state(1,i1).T = 1173;
        F1.state(1,i1).mdot = 80/3600;
        
        % Helium
        F2 = fluid_class('Helium','WF','CP','BICUBIC&HEOS',1,5);
        F2.state(1,i2).p = 30e5;
        F2.state(1,i2).T = 813;
        F2.state(1,i2).mdot = 80/3600;
        
        % Set hex_mode
        hex_mode = 0;
        par = 0;
        
    case 10 % Nellis&Klein2009
        % Water
        F1 = fluid_class('Water','WF','CP','BICUBIC&HEOS',1,5);
        F1.state(1,i1).p = 1e5;
        F1.state(1,i1).T = 60 + 273.15;
        F1.state(1,i1).mdot = 0.03;
        
        % Nitrogen
        F2 = fluid_class('Nitrogen','ENV','CP','BICUBIC&HEOS',1,5);
        F2.state(1,i2).p = 1e5;
        F2.state(1,i2).T = 20 + 273.15;
        rho = RPN('PT_INPUTS',F2.state(1,i2).p,F2.state(1,i2).T,'D',F2);
        F2.state(1,i2).mdot = 0.06*rho;
        
        % Set hex_mode
        hex_mode = 0;
        par = 0;
        
    case 11 % Steam condenser (cross-flow air cooler)
        % Steam
        F1 = fluid_class('Water','WF','CP','HEOS',1,5);
        F1.state(1,i1).p = 0.1e5;
        F1.state(1,i1).Q = 0.8;
        F1.state(1,i1).mdot = 100;
        
        % Nitrogen
        F2 = fluid_class('Nitrogen','ENV','CP','BICUBIC&HEOS',1,5);
        F2.state(1,i2).p = 1e5;
        F2.state(1,i2).T = 10 + 273.15;
        F2.state(1,i2).mdot = 1e3;
        
        % Update fluid states
        [F1] = update(F1,[iL,i1],3);
        [F2] = update(F2,[iL,i2],1);
        
        % Set hex_mode
        %hex_mode = 1;
        %par = 0.5;
        hex_mode = 5;
        par = F1.state(1,i1).T-2;
        
    otherwise
        error('not implemented')
end
switch scenario
    case 11
    otherwise
        % Update fluid states
        [F1] = update(F1,[iL,i1],1);
        [F2] = update(F2,[iL,i2],1);
end

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
                shape = 'circular';
                
            case 5
                % For comparisson with scaling model (Hoopes2016)
                %NTU   = 4.45;
                %ploss = 0.02625;
                %D     = 2e-2;
                % For comparisson with numerical model (Hoopes2016)
                eff   = 0.9697;  %design value
                ploss = 0.02625;   %design value
                D1    = 0.611*2.00e-3; %design value
                %shape = 'circular';
                shape = 'PCHE';
                
            case 6
                eff   = 0.80;
                ploss = 1e-3;
                D1    = [];
                shape = 'cross-flow';
                
            case 7
                % For comparisson with Marchionni2019 - CFD
                eff   = 0.622;  %design value
                ploss = 0.001;   %design value
                D1    = 1.22e-3; %design value
                shape = 'PCHE';
                
            case 10
                % Nellis&Klein (these performance values are just guesses
                eff   = 0.521;
                ploss = 6e-5;
                D1    = [];
                shape = 'cross-flow';
                
            case 11
                % Nellis&Klein (these performance values are just guesses
                eff   = 0.80;
                ploss = 0.001;
                D1    = [];
                shape = 'cross-flow';
                
            otherwise
                eff   = 0.95;
                ploss = 0.01;
                D1    = 0.01;
                shape = 'circular';
        end
        
        par1  = eff;
        par2  = ploss;
        par3  = D1;
        par4  = shape;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CONSTRUCT HEAT EXCHANGER %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HX = hx_class(name, stage_type, 4, NX, 1, 1, model, par1, par2, par3, par4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify geometry manually
switch scenario        
    case 8 % Marchionni2019 - experimental
        HX.D1  = 1.22e-3;
        HX.D2  = 1.22e-3;
        HX.Af1 = 1.57e-6*54*42/2;
        HX.Af2 = 1.57e-6*54*42/2;
        HX.L   = 1.012;
        HX.A1  = 4*HX.Af1*HX.L/HX.D1;
        HX.A2  = 4*HX.Af2*HX.L/HX.D2;
        HX.shape = 'PCHE';
        HX.Lgeom_set = 1;
    
    case 9 % Figley2013
        HX.D1  = 1.222e-3;
        HX.D2  = 1.22e-3;
        HX.Af1 = 1.571e-6*10*12;
        HX.Af2 = 1.571e-6*10*12;
        HX.L   = 0.2472;
        HX.A1  = 4*HX.Af1*HX.L/HX.D1;
        HX.A2  = 4*HX.Af2*HX.L/HX.D2;
        HX.shape = 'PCHE';
        HX.Lgeom_set = 1;
    
    case 10 % Nellis&Klein2009 - cross-flow
        W1 = 0.2;
        H1 = 0.26;
        L1 = 0.06;
        D1 = 3.63e-3;
        Dout = 10.2e-3;
        t2 = 0.9e-3;
        D2 = Dout - 2*t2;
        Ntubes = 20;
        sigma1 = 0.534;
        sigma2 = 0.099;
        pfins  = 1/315;
        Nfins  = W1/pfins;
        Afins  = 2*Nfins*(H1*L1 - pi/4*Dout^2*Ntubes);
        Anotf  = pi*Dout*W1*Ntubes;
        A1tot  = Afins + Anotf;
        Afin_A1= Afins/A1tot;
        HX.D1  = D1;
        HX.D2  = D2;
        HX.Af1 = W1*H1*sigma1;
        HX.Af2 = (HX.D2)^2*pi/4;
        HX.L1  = L1;
        HX.L2  = W1*Ntubes;
        HX.A1  = 4*HX.Af1*HX.L1/HX.D1;
        HX.A2  = 4*HX.Af2*HX.L2/HX.D2;
        HX.Afin_A1 = Afin_A1;
        HX.shape = 'cross-flow';
        HX.Lgeom_set = 1;
end

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

switch scenario
    case 10
        ht_w1 = 3496; %W/(m2.K)
        ht_w2 = interp1(HX.H.T-273.15,HX.H.ht,40,'machima');
        ht_a1 = 43.7; %W/(m2.K)
        ht_a2 = interp1(HX.C.T-273.15,HX.C.ht,40,'machima');
        UA1   = 58.4;
        UA2   = HX.UA;
        Dp_a1 = 6.0; %Pa
        Dp_a2 = HX.DppC*1e5;
        Dp_w1 = NaN; %Pa
        Dp_w2 = HX.DppH*1e5;
        
        fprintf(1,'\n\n');
        fprintf(1,'                    %8s   %8s   %8s\n','Nellis','Result','Error')
        fprintf(1,'ht air   [W/m2/K] = %8.2f   %8.2f   %6.1f %%\n',ht_a1,ht_a2,abs(ht_a1-ht_a2)/ht_a1*100)
        fprintf(1,'ht water [W/m2/K] = %8.2f   %8.2f   %6.1f %%\n',ht_w1,ht_w2,abs(ht_w1-ht_w2)/ht_w2*100)
        fprintf(1,'UA total [W/K]    = %8.2f   %8.2f   %6.1f %%\n',UA1,UA2,abs(UA1-UA2)/UA1*100)
        fprintf(1,'Dp air   [Pa]     = %8.2f   %8.2f   %6.1f %%\n',Dp_a1,Dp_a2,abs(Dp_a1-Dp_a2)/Dp_a1*100)
        fprintf(1,'Dp water [Pa]     = %8.2f   %8.2f   %6.1f %%\n',Dp_w1,Dp_w2,abs(Dp_w1-Dp_w2)/Dp_w1*100)
end

% Make plots
plot_hex(HX,1,10,'C',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%% OFF-DESIGN PERFORMANCE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use easier nomenclature for inlet conditions
TC1 = HX.C.T(1);
TH2 = HX.H.T(HX.NX+1);
pC1 = HX.C.pin;
pH2 = HX.H.pin;

% Load and pre-process data from literature (when needed for comparisson)
% and create mass flow rate arrays
switch scenario
    case 1
        n = 20;
        mdot1 = F1.state(iL,i1).mdot*linspace(0.25,1.0,n)';
        mdot2 = F2.state(iL,i1).mdot*linspace(0.25,1.0,n)';
        %mdot2 = F2.state(iL,i1).mdot*ones(size(mdot1));
        
    case 5 % Hoopes2016
        data     = load('Validation_Hoopes2016.csv');
        
        mdot_dat = data(:,1);
        TH1_dat  = data(:,3) + 273.15;
        TC2_dat  = data(:,4) + 273.15;
        pH1_dat  = data(:,5) * 1e5;
        pC2_dat  = data(:,6) * 1e5;
        
        DT1_dat = TH1_dat - TC1;
        DT2_dat = TH2 - TC2_dat;
        DpH_dat = pH2 - pH1_dat;
        DpC_dat = pC1 - pC2_dat;
        
        mdot1   = mdot_dat;
        mdot2   = mdot_dat;
        
    case 9 % Figley2013
        data    = load('Validation_Figley2013.csv');
        
        mdot_dat = data(:,1)/3600;        
        %DpH_dat  = data(:,2)*1e3; %Correlation
        %DpC_dat  = data(:,4)*1e3; %Correlation
        DpH_dat  = data(:,3)*1e3; %CFD
        DpC_dat  = data(:,5)*1e3; %CFD
        ReH_dat  = data(:,7);
        ReC_dat  = data(:,6);
        htH_dat  = data(:,9);
        htC_dat  = data(:,8);
        U_dat    = data(:,11);
        
        mdot1 = mdot_dat;
        mdot2 = mdot_dat;
        
    otherwise
        mdot1 = F1.state(iL,i1).mdot;
        mdot2 = F2.state(iL,i1).mdot;
end
n = length(mdot1);

% COMPUTE OFF-DESIGN PERFORMANCE
% Allocate arrays (numerical results)
TH1_num = zeros(size(mdot1));
TC2_num = zeros(size(mdot1));
pH1_num = zeros(size(mdot1));
pC2_num = zeros(size(mdot1));
DpH_num = zeros(size(mdot1));
DpC_num = zeros(size(mdot1));
ReH_num = zeros(size(mdot1));
ReC_num = zeros(size(mdot1));
htH_num = zeros(size(mdot1));
htC_num = zeros(size(mdot1));
U_num   = zeros(size(mdot1));
% Allocate arrays (analytical results)
TH1_an = zeros(size(mdot1));
TC2_an = zeros(size(mdot1));
pH1_an = zeros(size(mdot1));
pC2_an = zeros(size(mdot1));
DpH_an = zeros(size(mdot1));
DpC_an = zeros(size(mdot1));
for im = 1:n
    % Select mass flow rate
    F1.state(iL,i1).mdot = mdot1(im);
    F2.state(iL,i1).mdot = mdot2(im);
    
    % Run HEX code
    [HX,F1,~,F2,~] = hex_func(HX,iL,F1,i1,F2,i2,hex_mode,par);
    
    % Extract data for plotting
    TH1_num(im) = HX.H.T(1);
    TC2_num(im) = HX.C.T(NX+1);
    pH1_num(im) = HX.H.p(1);
    pC2_num(im) = HX.C.p(HX.NX+1);
    DpH_num(im) = HX.H.pin-HX.H.p(1);
    DpC_num(im) = HX.C.pin-HX.C.p(HX.NX+1);
    ReH_num(im) = sum(0.5*(HX.H.Re(1:end-1)+HX.H.Re(2:end)).*HX.dL')./HX.L;
    ReC_num(im) = sum(0.5*(HX.C.Re(1:end-1)+HX.C.Re(2:end)).*HX.dL')./HX.L;
    htH_num(im) = sum(0.5*(HX.H.ht(1:end-1)+HX.H.ht(2:end)).*HX.dL')./HX.L;
    htC_num(im) = sum(0.5*(HX.C.ht(1:end-1)+HX.C.ht(2:end)).*HX.dL')./HX.L;
    U_num(im)   = HX.UA/HX.A1;
    
    switch scenario
        case 1
            % Generate analytical results
            [TH1_an(im),TC2_an(im),pH1_an(im),pC2_an(im)] = hex_analytic(HX,iL,F1,i1,F2,i2);
            DpH_an(im) = HX.H.pin-pH1_an(im);
            DpC_an(im) = HX.C.pin-pC2_an(im);
            
            % Print comparison
            fprintf(1,'\n      Numerical  Analytical\n')
            fprintf(1,'TH1 [C] = %8.1f   %8.1f\n',TH1_num(im)-273.15,TH1_an(im)-273.15)
            fprintf(1,'TC2 [C] = %8.1f   %8.1f\n',TC2_num(im)-273.15,TC2_an(im)-273.15)
            fprintf(1,'DppH    = %8.5f   %8.5f\n',DpH_num(im)/HX.H.pin,DpH_an(im)/HX.H.pin)
            fprintf(1,'DppC    = %8.5f   %8.5f\n',DpC_num(im)/HX.C.pin,DpC_an(im)/HX.C.pin)
    end
    
    % Make plots
    %plot_hex(HX,1,20,'C');
    %pause(2);
end
DT1_num = TH1_num - TC1;
DT2_num = TH2 - TC2_num;
DpH_num = pH2 - pH1_num;
DpC_num = pC1 - pC2_num;


% COMPARE WITH ANALYTICAL SOLUTIONS OR DATA FROM LITERATURE
switch scenario
    case 1
        % Plot figures
        figure(25)
        plot(mdot1,TH1_num,'ro',mdot1,TH1_an,'r')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet temperature (hot) [K]')
        legend('numerical','analytical')
        figure(26)
        plot(mdot1,TC2_num,'bo',mdot1,TC2_an,'b')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet temperature (cold) [K]')
        legend('numerical','analytical')
        figure(27)
        semilogy(mdot1,DpH_num,'ro',mdot1,DpH_an,'r',mdot1,DpC_num,'bo',mdot1,DpC_an,'b')
        xlabel('Mass flow rate [kg/s]')
        ylabel('$ \Delta p / p $')
        legend('numerical (hot)','analytical (hot)','numerical (cold)','analytical (cold)','Location','Best')
        
        %{
        % Compute errors
        DTmax  = TH2 - TC1;
        errTH1 = max(abs(num_TH1 - an_TH1)/DTmax*100);
        errTC2 = max(abs(num_TC2 - an_TC2)/DTmax*100);
        errpH1 = max(abs((num_pH1 - an_pH1)./num_pH1*100));
        errpC2 = max(abs((num_pC2 - an_pC2)./num_pC2*100));
        errDpH = max(abs((num_DpH - an_DpH)./num_DpH*100));
        errDpC = max(abs((num_DpC - an_DpC)./num_DpC*100));
        %}

    case 5 % Compare numerical results with data from Hoopes2016
        % Compute errors
        %%{
        
        errDT1 = abs(DT1_num - DT1_dat)./DT1_dat*100;
        errDT2 = abs(DT2_num - DT2_dat)./DT2_dat*100;
        errDpH = ((DpH_num - DpH_dat)./DpH_dat)*100;
        errDpC = ((DpC_num - DpC_dat)./DpC_dat)*100;
        %%}
        %{
        DTmax  = TH2 - TC1;
        errDT1 = abs(num_TH1-273.15 - data(:,3)')/DTmax*100;
        errDT2 = abs(num_TC2-273.15 - data(:,4)')/DTmax*100;
        errDpH = abs(num_pH1 - data(:,5)'*1e5)./pH2*100;
        errDpC = abs(num_pC2 - data(:,6)'*1e5)./pC1*100;
        %}
        
        % Plot figures
        figure(30)
        yyaxis left
        plot(mdot1,TH1_num,mdot_dat,TH1_dat,'s');
        %ylim([85 105])
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet temperature (hot) [K]')        
        yyaxis right
        plot(mdot1,errDT1,'d')
        %ylim([0 0.50])
        ylabel('Error [$\%$]')
        legend('Numerical','Hoopes2016','Error','Location','NorthWest')        
        
        figure(31)
        yyaxis left
        plot(mdot1,TC2_num,mdot_dat,TC2_dat,'s')
        %ylim([615 635])
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet temperature (cold) [K]')
        yyaxis right
        plot(mdot1,errDT2,'d')
        %ylim([0 0.50])
        ylabel('Error [$\%$]')
        legend('Numerical','Hoopes2016','Error','Location','Best')
        
        figure(32)
        yyaxis left
        plot(mdot1,pH1_num/1e5,mdot_dat,pH1_dat/1e5,'s')
        %ylim([29 31])
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet pressure (hot) [bar]')
        yyaxis right
        plot(mdot1,errDpH,'d')
        %ylim([0 0.10])
        ylabel('Error [$\%$]')
        legend('Numerical','Hoopes2016','Error','Location','Best')
        
        figure(33)
        yyaxis left
        plot(mdot1,pC2_num/1e5,mdot_dat,pC2_dat/1e5,'s')
        %ylim([290 300])
        xlabel('Mass flow rate [kg/s]')
        ylabel('Outlet pressure (cold) [bar]')
        yyaxis right
        plot(mdot1,errDpC,'d')
        %ylim([0 2.5])
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
        
    case 9 % Compare numerical results with data from Figley2013
        
        % Compute errors
        errDpH = ((DpH_num - DpH_dat)./DpH_dat)*100;
        errDpC = ((DpC_num - DpC_dat)./DpC_dat)*100;
        errReH = ((ReH_num - ReH_dat)./ReH_dat)*100;
        errReC = ((ReC_num - ReC_dat)./ReC_dat)*100;
        errhtH = ((htH_num - htH_dat)./htH_dat)*100;
        errhtC = ((htC_num - htC_dat)./htC_dat)*100;
        errU   = ((U_num - U_dat)./U_dat)*100;
        
        figure(40)
        yyaxis left
        plot(mdot1,DpH_num/1e3,mdot_dat,DpH_dat/1e3,'s')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Pressure drop (hot) [kPa]')
        yyaxis right
        plot(mdot1,errDpH,'d')
        ylabel('Error [$\%$]')
        legend('Numerical','Figley2013','Error','Location','Best')
        
        figure(41)
        yyaxis left
        plot(mdot1,DpC_num/1e3,mdot_dat,DpC_dat/1e3,'s')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Pressure drop (cold) [kPa]')
        yyaxis right
        plot(mdot1,errDpC,'d')
        ylabel('Error [$\%$]')
        legend('Numerical','Figley2013','Error','Location','Best')
        
        figure(42)
        yyaxis left
        plot(mdot1,ReH_num,mdot_dat,ReH_dat,'s')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Reynolds number (hot side)')
        yyaxis right
        plot(mdot1,errReH,'d')
        ylabel('Error [$\%$]')
        legend('Numerical','Figley2013','Error','Location','Best')
        
        figure(43)
        yyaxis left
        plot(mdot1,ReC_num,mdot_dat,ReC_dat,'s')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Reynolds number (cold side)')
        yyaxis right
        plot(mdot1,errReC,'d')
        ylabel('Error [$\%$]')
        legend('Numerical','Figley2013','Error','Location','Best')
        
        figure(44)
        yyaxis left
        plot(mdot1,htH_num,mdot_dat,htH_dat,'s')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Heat transfer coeff. (hot) [W/m$^2$/K]')
        yyaxis right
        plot(mdot1,errhtH,'d')
        ylabel('Error [$\%$]')
        legend('Numerical','Figley2013','Error','Location','Best')
        
        figure(45)
        yyaxis left
        plot(mdot1,htC_num,mdot_dat,htC_dat,'s')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Heat transfer coeff. (cold) [W/m$^2$/K]')
        yyaxis right
        plot(mdot1,errhtC,'d')
        ylabel('Error [$\%$]')
        legend('Numerical','Figley2013','Error','Location','Best')
        
        figure(46)
        yyaxis left
        plot(mdot1,U_num,mdot_dat,U_dat,'s')
        xlabel('Mass flow rate [kg/s]')
        ylabel('Heat transfer coefficient [W/m$^2$/K]')
        yyaxis right
        plot(mdot1,errU,'d')
        ylabel('Error [$\%$]')
        legend('Numerical','Figley2013','Error','Location','Best')
        
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