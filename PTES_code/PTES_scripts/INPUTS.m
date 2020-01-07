%%% SET INPUT VARIABLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call the correct input file
Load.mode = 5;
switch Load.mode
    case {0,1,2,3} % Joule-Bratyon PTES / Joule-Brayton + Rankine
        JB_RANK_INPUTS
    case {4, 5, 6} % sCO2-PTES type cycles
        sCO2_INPUTS
end

% Code options
multi_run  = 0; % run cycle several times with different parameters?
optimise   = 0; % optimise cycle?
make_plots = 1; % make plots?
save_figs  = 0; % save figures at the end?
make_hex_plots = 0; % make plots of heat exchangers?

if (Nc_ch > 1 || Ne_ch > 1) && (Ncld > 1 || Nhot > 1)
    error('Have not implemented multiple compressions/expansions AND multiple storage tanks in series')
end

% Set double tanks
if Ncld == 1
    fluidC = fluid_class(fCname,'SF','TAB',NaN,Load.num,30); % Storage fluid
    CT  = double_tank_class(fluidC,TC_dis0,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0,Load.num+1); %cold double tank
else
    for ii = 1 : Ncld
        fluidC(ii)  = fluid_class(char(fCname(ii,:)),'SF','TAB',NaN,Load.num,30);
        CT(ii)      = double_tank_class(fluidC(ii),TC_dis0(ii),p0,MC_dis0(ii),TC_chg0(ii),p0,MC_chg0(ii),T0,Load.num+1); %cold double tank
    end
end

% Hot tanks
if Nhot == 1
    fluidH = fluid_class(fHname,'SF','TAB',NaN,Load.num,30); % Storage fluid
    HT  = double_tank_class(fluidH,TH_dis0,p0,MH_dis0,TH_chg0,p0,MH_chg0,T0,Load.num+1); %hot double tank
else
    for ii = 1 : Nhot
        fluidH(ii)  = fluid_class(char(fHname(ii,:)),'SF','TAB',NaN,Load.num,30);
        HT(ii)  = double_tank_class(fluidH(ii),TH_dis0(ii),p0,MH_dis0(ii),TH_chg0(ii),p0,MH_chg0(ii),T0,Load.num+1); %hot double tank
    end
end

% Set 'atmospheric' air tanks
air  = fluid_class('Air','ENV','CP','HEOS',Load.num,30);
huge = max(Load.mdot)*3600*1e6; % represents a very large mass
AT   = double_tank_class(air,T0,p0,huge,T0,p0,huge,T0,Load.num+1);

% Use new heat exchanger calls?
new_hex_calls = 1;

% Heat rejection streams
environ = environment_class(T0,p0,Load.num,10);

% Variables to run cycle multiple times and plot curves. The variables must
% have been defined in the SET_MULTI_RUN script
if multi_run==1
    Vpnt = 'Ran_TbotC';  % variable along curve
    Npnt = 8;            % points on curve
    pnt1 = 2+273.15;     % min value
    pnt2 = 40+273.15;    % max value
    Apnt = linspace(pnt1,pnt2,Npnt); % array
    Vcrv = 'Ne_ch';  % variable between curves
    Acrv = 1;%[1,2,3];
    Ncrv = numel(Acrv);
    %     Vcrv = 'eta';  % variable between curves
    %     Ncrv = 3;      % number of curves
    %     crv1 = 0.95;   % min value
    %     crv2 = 0.99;   % max value
    %     Acrv = linspace(crv1,crv2,Ncrv); % array
else
    Npnt=1; Ncrv=1;
    % Start new logfile
    delete ./Outputs/log.txt
    diary  ./Outputs/log.txt
end

% Save initial values of Load structure (to be reset during the INITIALISE
% subroutine, in case any values changed during running time, e.g.
% discharge time)
Load0 = Load;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%