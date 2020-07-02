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

tic

% Other inputs
% Heat rejection streams
% Set atmospheric conditions and cycle parameters
T0      = 30 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa

stH       = 8 ;
Load.time = [stH;stH;].*3600;      % time spent in each load period, s
Load.type = ["chg";"dis"];    % type of load period
Load.mdot = [10;10];      % working fluid mass flow rate, kg/s

Load.num  = numel(Load.time);
Load.ind  = 1:Load.num;

environ = environment_class(T0,p0,Load.num,10);

% Fluid
FSname = 'MineralOil';  % fluid name
fluidS = fluid_class(FSname,'SF','TAB',NaN,Load.num,30); % Storage fluid
Nhot = 1 ;
% Packed bed inputs
PACKED_BED_INPUTS_orig

i    = 1 ;
iCYC = 1 ;
iF   = 1 ;
Lcyc = false ;
den  = 0.0 ;
en_prev = 0.0 ;

for ii = 1 : Nhot
    [pbH(ii), pbH(ii).H(iCYC), pbH(ii).S(iCYC)] = PB_ENERGY(pbH(ii), fluidS) ; % Evaluate energy at the start of time
end
iCYC = iCYC + 1;

while ~Lcyc 
    
    [pbH, TsC, TfC, fluidS, iF, iCYC] = PB_RUN_HEX(pbH, Nhot, fluidS, [iCYC,iF], environ, iCYC, 'chg');
    %[pbH, TsD, TfD, iCYC] = PB_RUN(pbH, Nhot, fluidS, iCYC, 'dis');
    fprintf(1,'COMPLETED CYCLE %5i\n\n',i) ;
    
    den     = 100.0 * abs(pbH(1).H(1) - en_prev) / en_prev ;
    en_prev = pbH(1).H(1) ;
    
    if i == pbH.Ncyc || den < 0.1 
        Lcyc = true ;
    end
    
    i = i + 1 ;
    iCYC = 2 ;

end

error1 = 100*(((pbH.Hflux(2,1)-pbH.Hflux(2,2)))-pbH.H(2))/((pbH.Hflux(2,1)-pbH.Hflux(2,2))) ;
fprintf(1,'FIRST LAW ERROR %8.3f\n\n',error1) ;

plot(TsC,'DisplayName','TsC')
toc