% This file should be split up further
% Set atmospheric conditions and cycle parameters
T0      = 50 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
pmax    = 90e5;         % top pressure, Pa
PRch    = 2.0;          % charge pressure ratio
PRr     = 0.85;          % discharge pressure ratio: PRdis = PRch*PRr
PRr_min = 0.1;          % minimum PRr for optimisation
PRr_max = 3.0;          % maximum PRr for optimisation
setTmax = 0;            % set Tmax? (this option substitutes PRch)
Tmax    = 570 + 273.15; % maximum temp at compressor outlet, K

% Set component parameters
eta   = 0.80;  % polytropic efficiency
eff   = 0.85;  % heat exchanger effectiveness
ploss = 0.02;  % pressure loss in HEXs

% Number of intercooled/interheated compressions/expansions
Nc_ch = 1; % number of compressions during charge
Ne_ch = 1; % number of expansions during charge
nH    = Nc_ch;        % number of hot fluid streams
nC    = Ne_ch;        % number of cold fluid streams

% Number of hot and cold stores IN SERIES
Ncld = 1; % number of cold stores. Not implemented for >2
Nhot = 1; % number of hot stores. Not implemented for >2

switch Load.mode
    case 4
        fac = 3.5211 ; % THis can be used to more easily set the mass flow to obtain a desired power output
        Load.time = [10;4;10].*3600;          % time spent in each load period, s
        Load.type = ["chgCO2";"str";"disCO2"]; % type of load period
        Load.mdot = [1000*fac;0;1000*fac];              % working fluid mass flow rate, kg/s
    case 5
        Load.time = [10;4;10].*3600;          % time spent in each load period, s
        Load.type = ["sol";"str";"rcmpCO2"]; % type of load period
        Load.mdot = [10;0;10];              % working fluid mass flow rate, kg/s
    case 6
        Load.time = [10;4;10].*3600;          % time spent in each load period, s
        Load.type = ["chgTSCO2";"str";"disTSCO2"]; % type of load period
        Load.mdot = [1000;0;1000];              % working fluid mass flow rate, kg/s
end
Load.num = numel(Load.time);
Load.ind = 1:Load.num;

Lcld    = false ;       % Make cold store as cold as possible?
Lrcmp   = false ;       % Is there a recompressor?
Trej    = 0.0 ;

% Number of recuperators
Nrcp = 0 ; % Can be 0,1,2.
if (Nrcp > 0) && (Nhot > 1 || Ncld > 1)
    error('Have not implemented recuperators and multiple storage tanks in series')
end

switch Nrcp
    % If the sCO2-PTES cycle is not recuperator, then efficiency is
    % increased by having several stores in series
    case 0
        % Hot storage tanks
        %fHname  = ["SolarSalt";"MineralOil"]; % fluid name
        %fHname  = ["MineralOil";"MineralOil"]; % fluid name
        fHname  = 'Water'; % fluid name
        MH_dis0(1:Nhot) = 1e9;          % initial mass of discharged hot fluid, kg
        MH_chg0(1:Nhot) = 0.00*MH_dis0; % initial mass of charged hot fluid, kg
        
        Td = T0+14; % Temperature of the coldest 'hot' tank when it is discharged
        Tc = 200 + 273.15; % Temperature of the hottest 'hot' tank when it is charged
        Tint = [125 + 273.15] ; % intermediate temperature between tanks. This array should be size = Nhot - 1
        TH_chg0 = zeros(1,Nhot);
        TH_dis0 = zeros(1,Nhot);
        % Allocate temperatures to the correct tanks
        if Nhot == 1
            TH_chg0(1) = Tc ;
            TH_dis0(1) = Td ;
        else
            for ii = 1 : Nhot
                if ii == 1
                    TH_chg0(ii) = Tc ;
                    TH_dis0(ii) = Tint(ii) ;
                elseif ii == Nhot
                    TH_chg0(ii) = Tint(ii-1) ;
                    TH_dis0(ii) = Td ;
                else
                    TH_chg0(ii) = Tint(ii-1) ;
                    TH_dis0(ii) = Tint(ii) ;
                end
            end
        end
        % Cold storage tanks
        %fCname  = ["INCOMP::MEG2[0.56]";"INCOMP::MEG2[0.56]"]; % fluid name
        fCname  = 'Water'; % fluid name
        MC_dis0(1:Ncld) = 1e9;          % initial mass of discharged cold fluid, kg
        MC_chg0(1:Ncld) = 0.00*MC_dis0; % initial mass of charged cold fluid, kg
        
        % Cold store temperatures
        Td = 115 + 273.15; % Temperature of the hottest 'cold' tank when it is discharged
        Tc = T0-5; % Temperature of the coldest 'cold' tank when it is charged
        Tint = [55 + 273.15] ; % intermediate temperature between tanks. This array should be size = Nhot - 1
        TC_chg0 = zeros(1,Ncld);
        TC_dis0 = zeros(1,Ncld);
        % Allocate temperatures to the correct tanks
        if Ncld == 1
            TC_chg0(1) = Tc ;
            TC_dis0(1) = Td ;
        else
            for ii = 1 : Ncld
                if ii == 1
                    TC_chg0(ii) = Tc ;
                    TC_dis0(ii) = Tint(ii) ;
                elseif ii == Ncld
                    TC_chg0(ii) = Tint(ii-1) ;
                    TC_dis0(ii) = Td ;
                else
                    TC_chg0(ii) = Tint(ii-1) ;
                    TC_dis0(ii) = Tint(ii) ;
                end
            end
        end
    case 1
        % Hot storage tanks
        fHname  = 'MineralOil';  % fluid name
        TH_dis0 = 100 + 273.15;  % initial temperature of discharged hot fluid, K
        MH_dis0 = 1e9;          % initial mass of discharged hot fluid, kg
        TH_chg0 = 200 + 273.15; % initial temperature of charged hot fluid, K
        MH_chg0 = 0.00*MH_dis0; % initial mass of charged hot fluid, kg
        % Cold storage tanks
        fCname  = 'INCOMP::MEG2[0.56]'; % fluid name
        TC_dis0 = T0 + 0;           % initial temperature of discharged cold fluid, K
        MC_dis0 = 1e9;          % initial mass of discharged cold fluid, kg
        TC_chg0 = T0-5;        % initial temperature of charged cold fluid, K
        MC_chg0 = 0.00*MC_dis0; % initial mass of charged cold fluid, kg
        
        % If there are two recuperators, also use a recompressor during discharge
    case 2
        Lrcmp    = true ;         % Is there a recompression
        switch Load.mode
            case 4
                % Hot storage tanks
                fHname  = 'SolarSalt';  % fluid name
                TH_dis0 = 410. + 273.15;  % initial temperature of discharged hot fluid, K
                MH_dis0 = 1e9;              % initial mass of discharged hot fluid, kg
                TH_chg0 = 570 + 273.15;     % initial temperature of charged hot fluid, K
                MH_chg0 = 0.0*1.e6;         % initial mass of charged hot fluid, kg
                % Cold storage tanks
                fCname  = 'INCOMP::MEG2[0.56]'; % fluid name
                TC_dis0 = T0 + 0;           % initial temperature of discharged cold fluid, K
                MC_dis0 = 1e9;          % initial mass of discharged cold fluid, kg
                TC_chg0 = T0-5;        % initial temperature of charged cold fluid, K
                MC_chg0 = 0.0*1.e6; % initial mass of charged cold fluid, kg
                % Choose a threshold temperature between the recuperators
                TthreshC = 38. + 273.15 ; % Charge - threshold is on low-pressure side
                TthreshD = 200. + 273.15 ; % Discharge - threshold is on high-pressure side
                Trej     = 3.5 ;  % Reject heat to this temperature above ambient (helps cold store behaviour)
            case 5
                % Hot storage tanks
                fHname  = 'SolarSalt';  % fluid name
                TH_dis0 = 410. + 273.15;    % initial temperature of discharged hot fluid, K
                MH_dis0 = 1e9;              % initial mass of discharged hot fluid, kg
                TH_chg0 = 565 + 273.15;     % initial temperature of charged hot fluid, K
                MH_chg0 = 0.0*1.e6;         % initial mass of charged hot fluid, kg
                % Cold storage tanks
                fCname  = 'INCOMP::MEG2[0.56]'; % fluid name
                TC_dis0 = T0 + 0;           % initial temperature of discharged cold fluid, K
                MC_dis0 = 1e9;              % initial mass of discharged cold fluid, kg
                TC_chg0 = T0-5;             % initial temperature of charged cold fluid, K
                MC_chg0 = 0.0*1.e6;         % initial mass of charged cold fluid, kg
                % Choose a threshold temperature between the recuperators
                TthreshD = 250. + 273.15 ; % Discharge - threshold is on high-pressure side
            case 6
                % Time-shifted recompression sCO2 power cycles
                % Hot storage tanks
                Nhot = 2 ;
                Ncld = 1 ;
                fHname  = ["SolarSalt";"MineralOil"]; % fluid name
                MH_dis0(1:Nhot) = 1e9;          % initial mass of discharged hot fluid, kg
                MH_chg0(1:Nhot) = 0.00*MH_dis0; % initial mass of charged hot fluid, kg
                
                % Temperatures of tanks (some of these will be reset later)
                TH_dis0(1) = 410. + 273.15;    % initial temperature of discharged hot fluid, K
                TH_chg0(1) = 565. + 273.15;
                
                TH_dis0(2) = T0;    % initial temperature of discharged hot fluid, K
                TH_chg0(2) = 200 + 273.15 ;
        
                % Cold storage tanks
                fCname  = 'INCOMP::MEG2[0.56]'; % fluid name
                MC_dis0(1:Ncld) = 1e9;          % initial mass of discharged cold fluid, kg
                MC_chg0(1:Ncld) = 0.00*MC_dis0; % initial mass of charged cold fluid, kg
        
                % Cold store temperatures
                TC_dis0 = 150 + 273.15; % Temperature of the hottest 'cold' tank when it is discharged
                TC_chg0 = T0-5; % Temperature of the coldest 'cold' tank when it is charged
                
                % Choose a threshold temperature between the recuperators
                TthreshD = 250. + 273.15 ; % Discharge - threshold is on high-pressure side
                Lcld = false ;
        end
                
end


% Set working fluids, storage fluids, and heat rejection streams. 'WF' or
% 'SF' indicates working fluid or storage fluid. 'CP' or 'TAB' indicate
% CoolProp or Tabular reading modes. 'backend' is used by CoolProp (either
% 'HEOS', 'TTSE' or 'BICUBIC&HEOS'). 'num' indicates number of preallocated
% elements in state arrays.
% Working fluid
gas = fluid_class('CarbonDioxide','WF','CP','TTSE',Load.num,30);

% Make heat exchangers


% Heat exchangers set up to match Ty's work
if Load.mode == 5
    HX(1) = hx_class('hot', 'hex', 'eff', 0.879, ploss, 25, 100, Load.num, Load.num) ; % Hot heat exchanger
    HX(2) = hx_class('cold', 'hex', 'eff', eff, ploss, 25, 100, Load.num, Load.num) ; % Hot heat exchanger
    HX(3) = hx_class('regen', 'regen', 'eff', 0.968, ploss, 25, 100, Load.num, Load.num) ; % Hot heat exchanger
    HX(4) = hx_class('regen', 'regen', 'eff', 0.937, ploss, 25, 100, Load.num, Load.num) ; % Hot heat exchanger
elseif Load.mode ==6
    HX(1) = hx_class('hot', 'hex', 'eff', eff, ploss, 0, 100, Load.num, Load.num) ; % Hot heat exchanger
    HX(2) = hx_class('hot', 'hex', 'eff', eff, ploss, 25, 100, Load.num, Load.num) ; % Hot heat exchanger
    HX(3) = hx_class('cold', 'hex', 'eff', eff, ploss, 25, 100, Load.num, Load.num) ; % Cold heat exchanger
    HX(4) = hx_class('regen', 'regen', 'eff', eff, ploss, 0, 100, Load.num, Load.num) ; % Recuperator exchanger
    HX(5) = hx_class('regen', 'regen', 'eff', eff, ploss, 0, 100, Load.num, Load.num) ; % Recuperator heat exchanger    
else
    iHX = 1 ; % Heat exchanger counter
    for ii = 1 : Nhot
        HX(iHX) = hx_class('hot', 'hex', 'eff', eff, ploss, 25, 100, Load.num, Load.num) ; % Hot heat exchanger
        iHX = iHX + 1 ;
    end
    for ii = 1 : Ncld
        HX(iHX) = hx_class('cold', 'hex', 'eff', eff, ploss, 25, 100, Load.num, Load.num) ; % Hot heat exchanger
        iHX = iHX + 1 ;
    end
    if (Nhot < 2) && (Ncld < 2) && (Nrcp > 0)
        for ii = 1 : Nrcp
            HX(iHX) = hx_class('regen', 'regen', 'eff', eff, ploss, 25, 100, Load.num, Load.num) ; % Hot heat exchanger
            iHX = iHX + 1 ;
        end
    end
end

% Options for specifying heat exchanger geometry
% This will probably be expanded over time
for i = 1 : length(HX)
   HX(i).LestA = true ;
   HX(i).D1    = 0.025 ;
end



% Save copy of input file in "Outputs" folder
copyfile(['./PTES_scripts/',mfilename,'.m'],'./Outputs/')