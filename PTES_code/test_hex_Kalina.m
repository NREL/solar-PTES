%%% To run this script, place it inside the PTES_code folder, alongside the
%%% PTES.m file

%%%%%%%%%%%%%%%%%%%%%
%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clear;

% Enter debugging mode if an error occurs
dbstop if error

% Add path to solarPTES folder
folderName = '~/Dropbox/Work/Matlab/solar-PTES/PTES_code';
addpath([folderName,'/Classes/'],[folderName,'/Generic/'],...
    [folderName,'/Functions/'],[folderName,'/LIB/']);

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop


%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%

% Set indices
iL = 1; i1 = 1; i2 = 1;

% Declare fluids and specify initial conditions
% Propane
F1 = fluid_class('AlkaneMix','WF','TAB',NaN,1,5);
F1.state(iL,i1).p = 20e5;
F1.state(iL,i1).T = 60.8+273.15;
F1.state(iL,i1).mdot = 120;

% Water
F2 = fluid_class('Water','WF','CP','HEOS',1,5);
F2.state(iL,i2).p = 1e5;
F2.state(iL,i2).T = 11.7+273.15;
F2.state(iL,i2).mdot = F1.state(iL,i1).mdot*0.63;

% Set hex_mode and stage_type
hex_mode = 0.0;
par = 1.00;

% Update fluid states
[F1] = update(F1,[iL,i1],1);
[F2] = update(F2,[iL,i2],1);

% Specify HX settings
NX   = 100; % Number of sections (grid)
name = 'hot';
stage_type = 'hex';
model = 'geom';

% Specify performance objectives and channel type
eff = (58.8-11.7)/(60.8-11.7); %=0.9593;
ploss = 0.06;
D1    = 0.002;
shape = 'circular';
par1  = eff;
par2  = ploss;
par3  = D1;
par4  = shape;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DESIGN HEAT EXCHANGER %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do it a number of times, for different values of ploss
np = 10;
%plossV = logspace(log10(0.01),log10(0.1),np);
plossV = linspace(0.01,0.1,np);
L   = zeros(1,np);
A   = zeros(1,np);
Af  = zeros(1,np);
L_H = zeros(1,np);
h1  = zeros(1,np);
h2  = zeros(1,np);
U   = zeros(1,np);
Cf1  = zeros(1,np);
Cf2  = zeros(1,np);
Re1 = zeros(1,np);
Re2 = zeros(1,np);
Dpp1 = zeros(1,np);
Dpp2 = zeros(1,np);
Wpump = zeros(1,np);
paras = zeros(1,np);
for ip=1:np
    ploss = plossV(ip);
    par2  = ploss;
    
    % Construct HX object
    HX = hx_class(name, stage_type, 4, NX, 1, 1, model, par1, par2, par3, par4);
    
    % Find geometry that satisfies performance objectives under design
    % conditions
    [HX,F1,~,F2,~] = hex_func(HX,iL,F1,i1,F2,i2,hex_mode,par);
    
    % Compute work to pump water (compared to 5 MW discharge power)
    Dp2    = HX.DppC*F2.state(iL,i2).p;
    mdot2  = F2.state(iL,i2).mdot;
    rho2   = F2.state(iL,i2).rho;
    Wpump(ip) = Dp2*mdot2/(rho2*0.70);
    paras(ip) = Wpump(ip)/5e6;
    
    % Extract outputs
    L(ip)   = HX.L1; 
    A(ip)   = HX.A1;
    Af(ip)  = HX.Af1;
    L_H(ip) = HX.L1/sqrt(HX.Af1);
    h1(ip)  = mean(HX.H.ht);
    h2(ip)  = mean(HX.C.ht);
    U(ip)   = HX.UA/HX.A1;
    Cf1(ip)  = mean(HX.H.Cf);
    Cf2(ip)  = mean(HX.C.Cf);
    Re1(ip) = mean(HX.H.Re);
    Re2(ip) = mean(HX.C.Re);
    Dpp1(ip)= HX.DppH;
    Dpp2(ip)= HX.DppC;
end


%% MAKE PLOTS
figure(1)
plot(plossV*100,h1,plossV*100,h2,plossV*100,U)
xlabel('$(\Delta p/p)_{\mathrm{max}}$ [$$\%$$]')
ylabel('Heat transfer coeff. [W/m$$^2$$/K]')
legend('$$h_\mathrm{mix}$$','$$h_\mathrm{H2O}$$','$$U$$','Location','Best')

figure(2)
plot(plossV*100,Re1,plossV*100,Re2)
xlabel('$(\Delta p/p)_{\mathrm{max}}$ [$$\%$$]')
ylabel('Reynolds number')
legend('$$\mathrm{Re_{mix}}$$','$$\mathrm{Re_{H2O}}$$','Location','Best')

figure(3)
yyaxis left
plot(plossV*100,Wpump/1e3)
xlabel('$(\Delta p/p)_{\mathrm{max}}$ [$$\%$$]')
ylabel('Pump power [kW]')
ylim([0 Inf])
yyaxis right
plot(plossV*100,paras*100,'--')
ylabel('Parasitic [$$\%$$]')
ylim([0 Inf])

figure(4)
plot(plossV*100,Dpp1*100,plossV*100,Dpp2*100)
xlabel('$(\Delta p/p)_{\mathrm{max}}$ [$$\%$$]')
ylabel('$(\Delta p/p)$ [$$\%$$]')
ylim(xlim())
legend('$$(\Delta p/p)_{\mathrm{mix}}$$','$$(\Delta p/p)_{\mathrm{H2O}}$$','Location','Best')

figure(5)
plot(plossV*100,Cf1,plossV*100,Cf2)
xlabel('$(\Delta p/p)_{\mathrm{max}}$ [$$\%$$]')
ylabel('$c_f$')
legend('$$c_{f,\mathrm{mix}}$$','$$c_{f,\mathrm{H2O}}$$','Location','Best')

% Print summary
fprintf(1,'\n')
fprintf(1,'        %8s  %9s\n','Objective','Result')
fprintf(1,'Eff   = %8.3f   %9.3f\n',eff,1-min(HX.H(1).T-HX.C(1).T)/(HX.H(1).T(end)-HX.C(1).T(1)))
fprintf(1,'DppH  = %8.5f   %9.5f\n',ploss,HX.DppH)
fprintf(1,'DppC  = %8.5f   %9.5f\n',ploss,HX.DppC)

%print_hexs(HX,1,'Summary:\n');

fprintf('\n- Geometry -\n')
fprintf('L     = %6.1f m\n',HX.L1)
fprintf('A     = %6.1f m2\n',HX.A1)
fprintf('Af    = %6.1f m2\n',HX.Af1)
fprintf('L/H   = %6.1f \n',HX.L1/sqrt(HX.Af1))

fprintf('\n- Performance -\n')
fprintf('h1     = %6.1f W/m2/K\n',mean(HX.H.ht))
fprintf('h2     = %6.1f W/m2/K\n',mean(HX.C.ht))
fprintf('Re1    = %6.1f \n',mean(HX.H.Re))
fprintf('Re2    = %6.1f \n\n',mean(HX.C.Re))

fprintf('U      = %6.1f W/m2/K\n',HX.UA/HX.A1)
fprintf('DppMAX = %6.3f \n',max([HX.DppH,HX.DppC]))
fprintf('Wpump  = %6.1f W\n',Wpump(end))
fprintf('paras  = %6.3f %%\n',paras(end)*100)

% Make plots
plot_hex(HX,1,10,'C',true);

