%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% Author: Josh McTigue, Pau Farres-Antunez
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clear;

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./steamPCM_scripts/','./Other/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\steamPCM_scripts\','.\Other\')
end
% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop
tic

% Inputs
etaC        =   0.8 ;   % Compression efficiency
etaT        =   0.0 ;   % Turbine efficiency (set to 0 if throttle)
etaF        =   0.6 ;   % Fan efficiency
etaP        =   0.6 ;   % Water pump efficiency
effE        =   0.8 ;   % Effectiveness of evaporator
effC        =   0.5 ;   % Effectiveness of condenser

dP_water1   = 0.02 ;    % Pressure drop of water in condenser
dP_water2   = 0.10 ;    % Pressure drop of water around house
dP_air      = 0.01 ;    % Pressure drop of air through fan
dP_HP       = 0.02 ;    % Pressure drop in each condenser/evaporator of the heat pump

degC        =   273.15 ;
T0          =   15 + degC ;     % Outside ambient temperature
P0          =   1.013e5 ;       % Ambient pressure
Ti          =   21 + degC ;     % Required inside temperature

Te          =   1 + degC ;      % Evaporation temperature of working fluid (should be above freezing)
Tmax        =   80 + degC ;     % Max temperature to compress to
Twat        =   30 + degC ;     % Return temperature of water from house
Tlowest     =   -5 + degC ;     % Expect ambient temperature to be greater than this for 98% of the year

HPname   = 'R410A' ;            % Working fluid in the heat pump
HPmdot   = 1.0 ;                % Mass flow rate of working fluid  

% Set up fluid classes
fld   = fluid_class(HPname,'WF','CP','BICUBIC&HEOS',2,30);  % Heat pump working fluid
water = fluid_class('water','WF','CP','BICUBIC&HEOS',2,30); % Water
air   = fluid_class('Nitrogen','ENV','CP','BICUBIC&HEOS',2,30); % Atmopheric air

% Set up initial air conditions
air.state(1,1).p = P0; 
air.state(1,1).T = T0;
[air] = update(air,[1,1],1);

% Set up initial water conditions
water.state(1,1).p = P0 * 1.5; 
water.state(1,1).T = Twat;
[water] = update(water,[1,1],1);

% Set up working fluid conditions at pump inlet
fld.state(1,1).T = Te ;
fld.state(1,1).Q = 1.0 ;
fld.state(1,1).p = RP1('QT_INPUTS',1.0,fld.state(1,1).T,'P',fld) ;
[fld] = update(fld,[1,1],3);

% Set up classes for compressors, expanders, and pumps
CMP = compexp_class('comp', 'poly', 1, etaC, 2) ; % Compressor
EXP = compexp_class('exp',  'poly', 1, etaT, 2) ; % Expanders
PMP = compexp_class('pump', 'isen', 1, etaP, 2) ; % Water pump
FAN = compexp_class('comp', 'isen', 1, etaF, 2) ; % Air fan

%%% HEAT PUMP %%%
iHP = 1 ; % Counter as go around the heat pump

% Compress fluid
[CMP,fld,iHP] = compexp_func (CMP,1,fld,iHP,'Taim',Tmax, 1) ;

% Condense fluid
iHP = iHP + 1 ;
fld.state(1,iHP).p = fld.state(1,iHP-1).p * (1. - dP_HP) ;
fld.state(1,iHP).Q = 0 ;
[fld] = update(fld,[1,iHP],3);

% Expand fluid
p_aim         = fld.state(1,1).p / (1. - dP_HP) ;
[EXP,fld,iHP] = compexp_func (EXP,1,fld,iHP,'Paim',p_aim, 1) ; 

% Evaporate fluid
iHP = iHP + 1 ;
fld.state(1,iHP).p = fld.state(1,iHP-1).p * (1. - dP_HP) ;
fld.state(1,iHP).Q = 1 ;
[fld] = update(fld,[1,iHP],3);

% SPECIFIC HEAT AND WORK TERMS FOR HEAT PUMP
win = -(CMP.w(1) + EXP.w(1)) ;                  % specific work input, J/kg
qe  = fld.state(1,1).h - fld.state(1,4).h ;     % Heat into evaporator, J/kg
qc  = fld.state(1,2).h - fld.state(1,3).h ;     % Heat out of condenser evaporator, J/kg
COP = qc / win ;                                % Coefficient of performance - OF THE CLOSED LOOP PUMP

% CONSIDER AIR STREAMS
air.state(1,2).p = air.state(1,1).p * (1. - dP_air); 
air.state(1,2).T = air.state(1,1).T * (1. - effE) + effE * fld.state(1,4).T ; % From effectiveness method
[air] = update(air,[1,2],1);

% Calculate air mass flow rate (relative to heat pump working fluid)
cp = 0.5 * (RP1('PT_INPUTS',air.state(1,1).p,air.state(1,1).T,'CPMASS',air) ...
          + RP1('PT_INPUTS',air.state(1,2).p,air.state(1,2).T,'CPMASS',air)) ;
      
air.state(1,1).mdot = qe / (cp * (air.state(1,1).T - air.state(1,2).T)) ;

% For good measure
airNTU = -log(1. - effE) ; % NTU based on system with boiling/condensation
airUA  = airNTU * air.state(1,1).mdot * cp;

% Have to run a fan through the pressure drop
p_aim       = air.state(1,1).p  ;
[FAN,air,~] = compexp_func (FAN,1,air,2,'Paim',p_aim, 1) ; 
win_fan     = -FAN.w(1) * air.state(1,1).mdot ; % specific work input to fan, per mass flow of heat pump fluid

% CONSIDER WATER STREAM
water.state(1,2).p = water.state(1,1).p * (1. - dP_water1); 
water.state(1,2).T = water.state(1,1).T * (1. - effC) + effC * fld.state(1,2).T ; % From effectiveness method
[water] = update(water,[1,2],1);

% Calculate water mass flow rate (relative to heat pump working fluid)
cp = 0.5 * (RP1('PT_INPUTS',water.state(1,1).p,water.state(1,1).T,'CPMASS',water) ...
          + RP1('PT_INPUTS',water.state(1,2).p,water.state(1,2).T,'CPMASS',water)) ;
      
water.state(1,1).mdot = qc / (cp * (water.state(1,2).T - water.state(1,1).T)) ;

% For good measure
waterNTU = -log(1. - effC) ; % NTU based on system with boiling/condensation
waterUA  = waterNTU * water.state(1,1).mdot * cp;

% Water then delivers heat to house and goes through another pressure drop
water.state(1,3).p = water.state(1,2).p * (1. - dP_water2); 
water.state(1,3).T = Twat ;
[water] = update(water,[1,3],1);

% Now have to pump water back up to design pressure
p_aim         = water.state(1,1).p  ;
PMP.mdot0     = water.state(1,1).mdot ;
[PMP,water,~] = compexp_func (PMP,1,water,3,'Paim',p_aim, 1) ; 
win_pmp       = -PMP.w(1) * water.state(1,1).mdot ; % specific work input to pump per mass flow of heat pump fluid

COP_para      = qc / (win + win_fan + win_pmp) ;

% PLOT T-Q DIAGRAM OF CONDENSER
Npnts = 100 ;
Hfld  = linspace(fld.state(1,3).h,fld.state(1,2).h,Npnts) ;
Pfld  = linspace(fld.state(1,3).p,fld.state(1,2).p,Npnts) ;
Tfld  = RP1('HmassP_INPUTS',Hfld,Pfld,'T',fld) ;

Hwater  = linspace(water.state(1,1).h,water.state(1,2).h,Npnts) ;
Pwater  = linspace(water.state(1,1).p,water.state(1,2).p,Npnts) ;
Twater  = RP1('HmassP_INPUTS',Hwater,Pwater,'T',water) ;

qcond = linspace(0,qc,Npnts) ;

TpinchC = min(Tfld - Twater) ; % Check pinch point
if TpinchC < 0
    warning('Non-physical heat transfer in condenser')
end

figure(1)
plot(qcond/qc,Twater-degC) ;
hold on
plot(qcond/qc,Tfld-degC) ;
hold off
xlabel('Heat transferred (normalized)')
ylabel('Temperature, $$^{\circ}$$C')
legend('Water','Heat transfer fluid','location','northwest') ;

% PLOT T-Q DIAGRAM OF EVAPORATOR
Hfld  = linspace(fld.state(1,4).h,fld.state(1,1).h,Npnts) ;
Pfld  = linspace(fld.state(1,1).p,fld.state(1,1).p,Npnts) ;
Tfld  = RP1('HmassP_INPUTS',Hfld,Pfld,'T',fld) ;

Hair  = linspace(air.state(1,2).h,air.state(1,1).h,Npnts) ;
Pair  = linspace(air.state(1,2).p,air.state(1,1).p,Npnts) ;
Tair  = RP1('HmassP_INPUTS',Hair,Pair,'T',air) ;

qevap = linspace(0,qe,Npnts) ;

TpinchE = min(Tair-Tfld) ; % Check pinch point
if TpinchE < 0
    warning('Non-physical heat transfer in evaporator')
end

figure(2)
plot(qevap/qe,Tair-degC) ;
hold on
plot(qevap/qe,Tfld-degC) ;
hold off
xlabel('Heat transferred (normalized)')
ylabel('Temperature, $$^{\circ}$$C')
legend('Air','Heat transfer fluid','location','northwest') ;


% Can I plot a T-s diagram of the bloody thing?
% Saturation curve
figure(3) 
Npnts = 100 ;
Pcrit = RP1(0,0,0,'Pcrit',fld);
Psat  = linspace(0.5e5,Pcrit*0.9999,Npnts) ;
Tfld  = RP1('PQ_INPUTS',Psat,zeros(Npnts,1),'T',fld) ;
Sfld  = RP1('PQ_INPUTS',Psat,zeros(Npnts,1),'S',fld) ;
plot(Sfld/1e3,Tfld-degC,'color',c_grey) ;
hold on

Tfld  = RP1('PQ_INPUTS',Psat,ones(Npnts,1),'T',fld) ;
Sfld  = RP1('PQ_INPUTS',Psat,ones(Npnts,1),'S',fld) ;
plot(Sfld/1e3,Tfld-degC,'color',c_grey) ;

% Compressor
Tfld = [fld.state(1,1).T,fld.state(1,2).T] ;
Sfld = [fld.state(1,1).s,fld.state(1,2).s] ;
plot(Sfld/1e3,Tfld-degC,'color',c_dark_blue) ;
Tmax = max(Tfld) ; Smax = max(Sfld) ;

% Condenser
Pfld = logspace(log10(fld.state(1,2).p),log10(fld.state(1,3).p),Npnts);
Hfld = linspace(fld.state(1,2).h,fld.state(1,3).h,Npnts);
Tfld = RP1('HmassP_INPUTS',Hfld,Pfld,'T',fld) ;
Sfld = RP1('HmassP_INPUTS',Hfld,Pfld,'S',fld) ;
plot(Sfld/1e3,Tfld-degC,'color',c_dark_blue) ;

% Expander
Tfld = [fld.state(1,3).T,fld.state(1,4).T] ;
Sfld = [fld.state(1,3).s,fld.state(1,4).s] ;
plot(Sfld/1e3,Tfld-degC,'color',c_dark_blue) ;
Tmin = min(Tfld) ; Smin = min(Sfld) ;

% Evaporator
Pfld = logspace(log10(fld.state(1,4).p),log10(fld.state(1,5).p),Npnts);
Hfld = linspace(fld.state(1,4).h,fld.state(1,5).h,Npnts);
Tfld = RP1('HmassP_INPUTS',Hfld,Pfld,'T',fld) ;
Sfld = RP1('HmassP_INPUTS',Hfld,Pfld,'S',fld) ;
plot(Sfld/1e3,Tfld-degC,'color',c_dark_blue) ;

hold off

xlabel('Entropy, kJ/kg.K')
ylabel('Temperature, $$^{\circ}$$C') ;

xlim ([0.9*Smin/1e3,1.05*Smax/1e3]) ;
ylim ([0.9*Tmin-degC,1.05*Tmax-degC]) ;


%%% MODEL HEAT LOSS FROM HOUSE %%%
Nfloor          = 2 ;   % Number of floors
Ndoor           = 2;    % Number of doors
Nwindow_front   = 4;    % Number of windows on front and back of hours
Nwindow_side    = 2;    % Number of windows on the sides of the house
pitch           = 30;   % Pitch of roof, degrees

area_tot        = 100 ; % Floor area, m2
area_win        = 1.5 ; % Each window area, m2
area_door       = 2 ;   % Each door area

draft           = 1 ;   % Number of air exchanges per hour

% Dimensions
area_floor = area_tot / Nfloor ;
L          = sqrt(area_floor) ;
W          = L ;
H          = 2.5 * Nfloor ;
V          = L * W * H ;

area_door  = area_door * Ndoor ;
area_win   = area_win * 2 * (Nwindow_front + Nwindow_side) ;
area_wall  = 4 * L * H - area_door - area_win ;
area_roof  = 2 * L * 0.5 * W / cos (pitch * pi /180) ;
area_roof  = area_roof + 2 * W * 0.5 * tan(pitch * pi / 180) ;

area_roof_south = L * 0.5 * W / cos (pitch * pi /180) ;

% Convection and radiation
hi = 2 ; % Convection inside, W/m2.K
ho = 10 ; % Convection outside, W/m2.K

sb = 5.67e-8 ;
radi = 4 * sb * Ti^3 ; % Radiation heat transfer coefficient inside
rado = 4 * sb * T0^3 ; % Radiation heat transfer coefficient outside

Ri = 1. / (hi + radi) ; % Internal thermal resistance, W/m2.K
Ro = 1. / (ho + rado) ; % External thermal resistance, W/m2.K

% These values from https://dothemath.ucsd.edu/2012/11/this-thermal-house/#:~:text=For%20many%20building%20materials%2C%20%CE%BA,degree%20Celsius%20presented%20across%20it.
Rwall  = 13.3 / 5.7 ; % Thermal resistance of wall, W/m2.K
Rfloor = 13.3 / 5.7; % Thermal resistance of floor, W/m2.K
Rroof  = 13.2 / 5.7; % Thermal resistance of roof, W/m2.K
Rwind  = 2.0 / 5.7;  % Thermal resistance of window, W/m2.K
Rdoor  = 3.0 / 5.7;  % Thermal resistance of door, W/m2.K

% Calculate the thermal admittance, W/K
Qwall  = area_wall / (Rwall + Ri + Ro) ;
Qfloor = area_floor / (Rfloor + Ri + Ro) ;
Qroof  = area_roof / (Rroof + Ri + Ro) ;
Qwind  = area_win / (Rwind + Ri + Ro) ;
Qdoor  = area_door / (Rdoor + Ri + Ro) ;

Qstruc = Qwall + Qfloor + Qroof + Qwind + Qdoor ;

% Now consider air exchanges
mflow = V * 1.225 / 3600 ; % mass exchange, kg/s
Qflow = mflow * 1.005e3 ; % Thermal admittance (mflow * heat capacity of air).

Qtot  = Qstruc + Qflow ;

Pmax  = Qtot * (Ti - Tlowest) ; % Maximum thermal power input required, W

% Mass flow rate of heat pump to deliver this power
fld.state(1,1).mdot   = Pmax / qc ;
water.state(1,1).mdot = water.state(1,1).mdot * fld.state(1,1).mdot ;
air.state(1,1).mdot   = air.state(1,1).mdot * fld.state(1,1).mdot ;

% Now find actual power inputs etc.
Qc  = fld.state(1,1).mdot * qc ;
Qe  = fld.state(1,1).mdot * qe ;
Win = fld.state(1,1).mdot * win ;

Win_fan = air.state(1,1).mdot * -FAN.w(1) ;
Win_pmp = water.state(1,1).mdot * -PMP.w(1) ;

Win_tot = Win + Win_fan + Win_pmp ; % Power to be supplied by solar

%%% WEATHER DATA ***
fload    =  '..\..\..\..\Personal\applications\mitsubishi\rothbury_results.csv' ;
dat      = readmatrix(fload,'Range',[2 2]) ;   % Read load file
Ndat     = length(dat(:,1)) ;                  % Number of load periods
Tamb     = dat(:,3) + degC ;
solar_pow = dat(:,1) * 1e3 ;

Tthresh   = 10 + degC ;

Pin_th       = zeros(Ndat,1) ;

for i = 1 : Ndat
   if Tamb(i) < Tthresh
       % Have to use heating. 
       Pin_th(i) = Qtot * (Ti - Tamb(i)) ;
       %Pin_th(i) = 150 * (Ti - Tamb(i)) ;
   end
end

Ein_th = sum(Pin_th) * 3600 ; % Thermal energy input required, W
Ein_e  = Ein_th / COP_para ; % Electrical energy input required, W
%Ein_e  = Ein_th / 3.5 ; % Electrical energy input required, W
Ein_e_solar = sum(solar_pow) * 3600 ; % Electrical energy available from solar, W


toc

