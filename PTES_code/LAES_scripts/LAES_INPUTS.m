%%% SET INPUT VARIABLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Atmospheric conditions and cycle parameters
T0      = 15 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
TC1_0   = T0-100;       % temperature of cold fluid (end of discharge), K
TC2_0   = T0-180;       % temperature of cold fluid (end of discharge), K
TH_0    = T0;           % temperature of hot fluid  (end of discharge), K
mdot    = 10;           % working fluid mass flow rate, kg/s
xAir    = 0.1;          % fraction of mdot for second stream (this should be optimised)
t_ch    = 10 * 3600;    % charge time, s
pLA     = 3e5;          % pressure of liquid air tanks, Pa
pmax    = 150e5;        % top pressure during charge, Pa
PRch    = pmax/p0;      % charge pressure ratio
PRr     = 1.0;          % PRdis = PRch*PRr
PRr_min = 0.9;
PRr_max = 1.1;
gasName     = 'Air'; % working fluid
fluidHname  = 'MineralOil'; % hot storage fluid
fluidC1name = 'Methanol'; % cold storage fluid
fluidC2name = 'Propane'; % cold storage fluid

% Set component parameters
eta    = 0.90;  % polytropic efficiency
eff    = 0.97;  % heat exchanger effectiveness
ploss  = 0.01;  % pressure loss in HEXs

% Set operation modes
mode       = 0; % cycle mode. 0:Basic, 1:+heat pump, 2:+regenerator
Nc_ch      = 3; % number of compressions during charge
Ne_dis = Nc_ch; % expansions during discharge

multi_run  = 0; % run cycle several times with different parameters?
optimise   = 0; % optimise cycle?
make_plots = 0; % make plots?
savefig    = 0; % save figures at the end?

% Variables to run cycle multiple times and plot curves. The variables must
% have been defined in the LAES_SET_MULTI_RUN script
if multi_run==1
    Vpnt = 'eta';   % variable along curve
    Npnt = 5;       % points on curve    
    pnt1 = 0.85;    % min value
    pnt2 = 0.95;    % max value
    Apnt = linspace(pnt1,pnt2,Npnt); % array
    Vcrv = 'eff';  % variable between curves
    Ncrv = 3;      % number of curves    
    crv1 = 0.95;   % min value
    crv2 = 0.99;   % max value
    Acrv = linspace(crv1,crv2,Ncrv); % array
    SET_MULTI_VAR;
else
    Npnt=1; Ncrv=1;
    % Start new logfile
    delete ./Outputs/log.txt
    diary  ./Outputs/log.txt
end

% Save copy of input file in "Outputs" folder
copyfile(['./LAES_scripts/',mfilename,'.m'],'./Outputs/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%