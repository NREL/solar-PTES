% Plot out characteristic curves of compressors and turbines
%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clear;

% Enter debugging mode if an error occurs
%dbstop if error
%dbclear all

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop

% Set atmospheric conditions and cycle parameters
T0      = 30 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
pr      = 2.0 ;         % Design pressure ratio

% Set component parameters
eta0   = 0.90;  % polytropic efficiency
dp0    = 0.06;  % Design pressure loss that occurs between comp and turb etc.

Load.time = [10; 10].*3600;               % time spent in each load period, s
Load.type = ["chg"; "chg"];    % type of load period
Load.mdot = [1; 0.7];                       % working fluid mass flow rate, kg/s
Load.num  = numel(Load.time);
Load.ind  = 1:Load.num;

%gas = fluid_class('CarbonDioxide','WF','CP','TTSE',Load.num,30);
gas = fluid_class('Nitrogen','WF','CP','BICUBIC&HEOS',Load.num,30);

CMP = compexp_class('comp', 'isen', 2, eta0, Load.num) ; % Compressors
EXP = compexp_class('exp', 'isen', 2, eta0, Load.num) ; % Expanders

% Gas inlet to compressor
gas.state(1,1).p    = p0; 
gas.state(1,1).T    = T0;
gas.state(1,1).mdot = Load.mdot(1);
[gas]               = update(gas,[1,1],1);


% Design points of compressor
CMP.pr(1)   = 1 ;
[CMP,gas,~] = compexp_func (CMP,1,gas,1,'Paim',p0*pr) ;

% Gas inlet to turbine
gas.state(1,3).p    = gas.state(1,2).p * (1.-dp0) ; 
gas.state(1,3).T    = 500 + 273.15;
gas.state(1,3).mdot = Load.mdot(1);
[gas]               = update(gas,[1,3],1);
 
% Design point of expander
EXP.pr(1)   = 1 ;
[EXP,gas,~] = compexp_func (EXP,1,gas,3,'Paim',p0/(1.-dp0)) ;

%*****
% Now let's try to figure out what happens when go to an off-design point
% Gas inlet to compressor
gas.state(2,1).p    = 0.7*p0; 
gas.state(2,1).T    = T0;
gas.state(2,1).mdot = Load.mdot(2);
[gas]               = update(gas,[2,1],1);

% Off-design points of compressor
[CMP,gas,~] = compexp_func (CMP,2,gas,1,'Paim',p0*pr) ;

% Gas inlet to turbine
% Assume pressure losses scale with mass flow squared
dp = dp0 * (Load.mdot(2) / Load.mdot(1))^2 ;
gas.state(2,3).p    = gas.state(2,2).p * (1.-dp) ; 
gas.state(2,3).T    = 500 + 273.15;
gas.state(2,3).mdot = Load.mdot(2);
[gas]               = update(gas,[2,3],1);
 
% Off-design point of expander
[EXP,gas,~] = compexp_func (EXP,2,gas,3,'Paim',p0/(1.-dp0)) ;


print_states(gas,1,1:4,Load);
print_states(gas,2,1:4,Load);

fprintf('Compressor design efficiency, %%:         %8.2f\n',100*CMP.eta(1));
fprintf('Compressor off-design efficiency, %%:     %8.2f\n',100*CMP.eta(2));
fprintf('Compressor design pressure ratio:        %8.2f\n',CMP.pr(1));
fprintf('Compressor off-design pressure ratio:    %8.2f\n\n',CMP.pr(2));

fprintf('Expander design efficiency, %%:         %8.2f\n',100*EXP.eta(1));
fprintf('Expander off-design efficiency, %%:     %8.2f\n',100*EXP.eta(2));
fprintf('Expander design pressure ratio:        %8.2f\n',EXP.pr(1));
fprintf('Expander off-design pressure ratio:    %8.2f\n\n',EXP.pr(2));

% Design
W0 = Load.mdot(1) * (EXP.w(1) + CMP.w(1))/1e3; % Work
H0 = Load.mdot(1) * (gas.state(1,3).h - gas.state(1,2).h)/1e3; % Heat
HEeff0 = 100 * W0 / H0 ; % Efficiency
% Off-design 
W1 = Load.mdot(2) * (EXP.w(2) + CMP.w(2))/1e3;    % work
H1 = Load.mdot(2) * (gas.state(2,3).h - gas.state(2,2).h)/1e3; % Heat
HEeff1 = 100 * W1 / H1 ; % Efficiency

fprintf('Design work, kW:       %8.2f\n',W0);
fprintf('Design heat, kW:       %8.2f\n',H0);
fprintf('Design efficiency, %%:%8.2f\n\n',HEeff0);
fprintf('Off-design work, kW:   %8.2f\n',W1);
fprintf('Off-design heat, kW:       %8.2f\n',H1);
fprintf('Off-design efficiency, %%:%8.2f\n\n',HEeff1)