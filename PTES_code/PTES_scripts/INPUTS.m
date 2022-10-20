%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMPLE INTERFACE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The 'simple interface' is to introduce new users to the code or to
% faciliate integration with other tools.
% The simple interface is specifically to model a molten-salt Joule-Brayton
% PTES cycle for a single design-point calculation.

% Set 'simple_interface' to 1 to use this mode.
% Then edit the variables in SET_SIMPLE_INTERFACE.m
% Ignore the rest of this file
simple_interface = 0;

% To use the detailed interface set 'simple_interface' to 0
% Then make edits to the file below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DETAILED INTERFACE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode 0: Joule-Brayton PTES cycle
% Mode 1: Joule-Brayton heat pump (charge only)
% Mode 2: Joule-Brayton heat engine (discharge only)
% Mode 3: Joule-Brayton heat pump with a steam-Rankine heat engine
% Mode 4: sCO2 PTES cycle
% Mode 5: sCO2 recompression cycle (discharge only)
% Mode 6: sCO2 'time-shifted' recompression cycle
% Mode 7: Steam-Rankine heat engine (discharge only)
%
% Mode 20: PTES-LAES. Combined cycle energy storage (CCES)

% Call the correct input file
Load.mode  = 0 ;
Loffdesign = 0 ; % 'L' for Logical. 0 just run design case. 1 run design case then off-design load cycle.
Lreadload  = 0 ;
PBmode     = 0 ; % Liquid stores = 0; Packed beds = 1; Heat exchangers between power cycle and a storage fluid, which then passes through packed beds = 2

switch Load.mode
    case {0,1,2,3,7} % Joule-Bratyon PTES / Joule-Brayton + Rankine
        JB_RANK_INPUTS
    case {4, 5, 6}   % sCO2-PTES type cycles
        sCO2_INPUTS
    case 20          % PTES-LAES combined cycle
        CCES_INPUTS
    otherwise
        error('not implemented')
end

% If Wdis_req is zero, then the power output will be based on the mass flow rates in JB_RANK_INPUTS.m
Wdis_req = 100e6 ; % Required average power output during discharge
HP_mult  = 1; % Heat pump multiplier - i.e. the power rating of the heat pump compared to Wdis_req
wf_mdot  = 1000 ; % Mass flow rate of working fluid during discharge. If Wdis_req > 0 then this is used as a first guess, otherwise, this is the mass flow that is actually used
dis_dur  = 10; % Discharge duration, h
str_dur  = 10; % Storage duration, h


% Set heat exchanger parameters
HX_model  = 'geom';
eff       = 0.95;  % heat exchanger effectiveness
ploss     = 0.01;  % pressure loss in HEXs
HX_D1     = 0.005; %hydraulic diameter
HX_shape  = 'circular';   %channel shape for counter-flow HEXs
HX_NX     = 100; % number of sections for HEX algorithm
% For cross-flow HEXs
effX      = 0.90; 
plossX    = 0.001;
HX_shapeX = 'cross-flow';

% Code options
multi_run   = 0; % run cycle several times with different parameters?
Lmulti_mdot = 0; % Read data from previous multirun to recalculate what the actual mass flow rates should be for a desired power
optimise    = 0; % optimise cycle?
make_plots  = 1; % make plots?
save_figs   = 0; % save figures at the end?
make_hex_plots = 0; % make plots of heat exchangers?

%if (Nc_ch > 1 || Ne_ch > 1) && (Ncld > 1 || Nhot > 1)
%    error('Have not implemented multiple compressions/expansions AND multiple storage tanks in series')
%end

switch PBmode
    case 0

        % Set double tanks
        if Ncld == 1
            fluidC = fluid_class(fCname,'SF','TAB',NaN,Load.num,30); % Storage fluid
            %fluidC = fluid_class(fCname,'SF','CP','HEOS',Load.num,30); % Storage fluid
            CT  = double_tank_class(fluidC,TC_dis0,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0,CTmode,Load.num+1); %cold double tank
        else
            for ii = 1 : Ncld
                fluidC(ii)  = fluid_class(char(fCname(ii,:)),'SF','TAB',NaN,Load.num,30); %#ok<*SAGROW>
                CT(ii)      = double_tank_class(fluidC(ii),TC_dis0(ii),p0,MC_dis0(ii),TC_chg0(ii),p0,MC_chg0(ii),T0,CTmode,Load.num+1); %cold double tank
            end
        end
        
        % Hot tanks
        if Nhot == 1
            fluidH = fluid_class(fHname,'SF','TAB',NaN,Load.num,30); % Storage fluid
            %fluidH = fluid_class(fHname,'SF','CP','HEOS',Load.num,30);
            HT  = double_tank_class(fluidH,TH_dis0,p0,MH_dis0,TH_chg0,p0,MH_chg0,T0,HTmode,Load.num+1); %hot double tank
        else
            for ii = 1 : Nhot
                fluidH(ii)  = fluid_class(char(fHname(ii,:)),'SF','TAB',NaN,Load.num,30);
                HT(ii)  = double_tank_class(fluidH(ii),TH_dis0(ii),p0,MH_dis0(ii),TH_chg0(ii),p0,MH_chg0(ii),T0,HTmode,Load.num+1); %hot double tank
            end
        end
        
        switch Load.mode
            case 20
                % Medium tanks
                fluidM = fluid_class(fHname,'SF','TAB',NaN,Load.num,30); % Storage fluid
                MT     = double_tank_class(fluidM,TM_dis0,p0,MM_dis0,TM_chg0,p0,MM_chg0,T0,HTmode,Load.num+1); %medium double tank
        end
        
    case 1
        PACKED_BED_INPUTS
        
    case 2
        error ('Not implemented')
end


% Set 'atmospheric' air tanks
%air  = fluid_class('Air','ENV','CP','HEOS',Load.num,30);
air  = fluid_class('Nitrogen','ENV','CP','BICUBIC&HEOS',Load.num,30);
AT   = double_tank_class(air,T0,p0,0,T0,p0,0,T0,ATmode,Load.num+1);

% Heat rejection streams
environ = environment_class(T0,p0,Load.num,10);

% Variables to run cycle multiple times and plot curves. The variables must
% have been defined in the SET_MULTI_RUN script
if multi_run==1
    % Set variable along curves
    Vpnt = 'HT_B_dT';  % variable along curve
    Npnt = 7;            % points on curve
    pnt1 = -20;    % min value
    pnt2 = 40;    % max value
    Apnt = linspace(pnt1,pnt2,Npnt); % array
   
    % Set variable between curves
    Vcrv = 'mdot_off';
    %Acrv = [250,350,450]+273.15;
    %Acrv = [0.6 0.8 1.0 1.2 1.4];
    Acrv = 1.0;
    Ncrv = numel(Acrv);
    
    if Lmulti_mdot
        multi_mdot = mdot_extract(Npnt,Ncrv,Wdis_req) ;
    end
    
    if ~exist('./Outputs/Multi_run','dir')
        mkdir('./Outputs/Multi_run')
    end
    
    % Delete previous files
    delete('./Outputs/Multi_run/*.mat')
    
    % Store information on the variables being changed along the multi-run
    % calls
    save('./Outputs/Multi_run/Multi_run_var.mat',...
        'Vpnt','Npnt','Apnt','Vcrv','Ncrv','Acrv','Loffdesign');
    
else
    Npnt=1; Ncrv=1;
    % Start new logfile
    delete ./Outputs/log.txt
    diary  ./Outputs/log.txt
end

% Save some initial values
setTmax0 = setTmax ;
pmax0    = pmax ;
Design_Load0 = Design_Load ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

txt = ['RUN DESIGN CASE FIRST:     \n';'RUN OFF-DESIGN LOAD CYCLES:\n'];
line = '---------------------------\n';