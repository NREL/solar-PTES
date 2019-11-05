% This is a nasty wrapper than can alternately call the steam turbine model
% developed by Kevin Ellingwood and the PTES heat pump and heat exchanger
% models developed by Pau Farres-Antunez. This wrapper calls functions and
% scripts which is admittedly quite horrible and prone to probems arising
% down the line.
%
% This wrapper models a heat pump that charges hot and cold storage. The
% hot storage may be used as the heat source to an existing steam Rankine
% cycle. The cold storage may be used to reduce the heat rejection
% temperature in the steam turbine.
%
% Josh McTigue, JoshuaDominic.McTigue@nrel.gov
% 9 September 2019

clear;
dbstop if error

% First step is to calculate the steam turbine behavior at a particular
% pressure and temperature, with no cooling applied to condenser
cd '.\PTES_KE'

% Determine Operating System
c = computer();

% Addpaths and load CoolProp
switch computer
    case 'GLNXA64' %Linux
        addpath('./Inputs/','./Classes/','./Generic/');
    case 'PCWIN64' %Windows
        addpath('.\Inputs\','.\Classes\','.\Generic\');   
end
load_coolprop
set_graphics

% Input data
test_input

T0design = 25.0 + 273.15 ; % Design ambient temperature
T0actual = 30.0 + 273.15 ; % Actual ambient temperature
P_cond = 0.17 ; %0.2051 is KE's value. 0.08 based on SEGS VI from Patnode for water cooled. For air-cooled get optimal Pcond of 0.17 bar (see p201)
dT1 = 10.; % Set a dT between cooling air and condensing steam

% Call steam turbine function
[~,~,~]	= Steam_fxn(P_cond) ;
test_data % Sort out data constructs

% Now model heat exchanger between steam generation and molten salt
% Divide into two heat exchangers (1) Preheater, Boiler, Superheater
% (2) Reheater. This seems to be the method of Patnode (and SEGS VI -
% Kearney)

% First, move into a new directory and get classes loaded
cd '../PTES_code' ;

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\')
end
load_coolprop
set_graphics

% Now set up fluid classes
% Water
WF1 = fluid_class('Water','WF','CP','TTSE',1,30); % storage fluid       
WF1.state(1,1).p    = BOILER.Cin.P ;
WF1.state(1,1).T    = BOILER.Cin.T - 5 ;
WF1.state(1,1).mdot = 1.0 * mdot_des;
 
% Molten Salt
SF1 = fluid_class('SolarSalt','SF','TAB',NaN,1,30); % working fluid
SF1.state(1,1).p    = 1.e5;
SF1.state(1,1).T    = SUPH.Cout.T + 25 ; %BOILER.Cout.T + 28.0 ;
SF1.state(1,1).mdot = 6 * mdot_des;

% Update initial states
[SF1] = update(SF1,[1,1],1);
[WF1] = update(WF1,[1,1],1);

% Specify HEX performance
eff   = 0.97;
ploss = 0.01;

% Iterate on solar salt mass flow rate to get correct water outlet temp
% Run HEX routine
mdotMIN = 1 * mdot_des;
mdotMAX = 10 * mdot_des;
ii = 0 ;
NN = 100 ;
er = 0 ;
f3 = @(mdotIN) hex_function (mdotIN, SUPH.Cout.T, 0, SF1,[1,1],WF1,[1,1],eff,0,ploss,'hex',4,0) ;

while er < 0.1 && ii<NN
    mdotMAX      = 10. * (NN-ii) * mdot_des / NN ;
    [mdotIN,~,~] = golden_search(f3,mdotMIN,mdotMAX,0.005,'Min') ;
    SF1.state(1,1).mdot = mdotIN ;
    [WF1,SF1,~,~] = hex_TQ_cond(WF1,[1,1],SF1,[1,1],eff,0,ploss,'hex',4,0);
    if ii == 0
        diff0 = (SUPH.Cout.T - WF1.state(1,2).T) ;
    end
    diff1 = (SUPH.Cout.T - WF1.state(1,2).T) ;
    er    = 100. * abs(diff1 - diff0) / diff0 ;
    ii = ii + 1 ;
end

SF1.state(1,1).mdot = mdotIN ;
[WF1,SF1,~,~] = hex_TQ_cond(WF1,[1,1],SF1,[1,1],eff,0,ploss,'hex',4,0);
% Plot T-Q diagram
plot_hex(SF1,[1,1],WF1,[1,1],100,2)

% Now consider reheater. % Water
WF1.state(1,3).p    = REH.Cin.P ;
WF1.state(1,3).T    = REH.Cin.T ;
WF1.state(1,3).mdot = 1.0 * mdot_des ;
 
% Molten Salt
SF1.state(1,3).p    = 1.e5 ;
SF1.state(1,3).T    = REH.Cout.T + 25.0 ;
SF1.state(1,3).mdot = 1.32 * mdot_des;

% Update initial states
[SF1] = update(SF1,[1,3],1);
[WF1] = update(WF1,[1,3],1);

% Iterate on solar salt mass flow rate to get correct water outlet temp
% Iterate on solar salt mass flow rate to get correct water outlet temp
% Run HEX routine

mdotMIN = 1 * mdot_des;
mdotMAX = 2 * mdot_des;
ii = 0 ;
NN = 100 ;
er = 0 ;
f4 = @(mdotIN) hex_function (mdotIN, REH.Cout.T, 0, SF1,[1,3],WF1,[1,3],eff,0,ploss,'hex',4,0) ;

while er < 0.1 && ii<NN
    mdotMAX      = 2. * (NN-ii) * mdot_des / NN ;
    [mdotIN,~,~] = golden_search(f4,mdotMIN,mdotMAX,0.005,'Min') ;
    SF1.state(1,3).mdot = mdotIN ;
    [WF1,SF1,~,~] = hex_TQ_cond(WF1,[1,3],SF1,[1,3],eff,0,ploss,'hex',4,0);
    if ii == 0
        diff0 = (REH.Cout.T - WF1.state(1,4).T) ;
    end
    diff1 = (REH.Cout.T - WF1.state(1,4).T) ;
    er    = 100. * abs(diff1 - diff0) / diff0 ;
    ii = ii + 1 ;
end
SF1.state(1,3).mdot = mdotIN ;
[WF1,SF1,~,~] = hex_TQ_cond(WF1,[1,3],SF1,[1,3],eff,0,ploss,'hex',4,0);
% Plot T-Q diagram
plot_hex(SF1,[1,3],WF1,[1,3],100,3)


% Now mix the outlet solar salt
SF1.state(1,5).p = min(SF1.state(1,2).p,SF1.state(1,4).p) ;
SF1.state(1,5).h = (SF1.state(1,2).mdot * SF1.state(1,2).h + ...
                   SF1.state(1,4).mdot * SF1.state(1,4).h) / ...
                  (SF1.state(1,2).mdot + SF1.state(1,4).mdot) ;
SF1.state(1,5).mdot = SF1.state(1,2).mdot + SF1.state(1,4).mdot ;
[SF1] = update(SF1,[1,5],2); % Use mode 2 as enthalpy and pressure are known


% Now the molten salt operating temperatures are known. Run the PTES code
% to ensure it runs between these limits
% Call golden search routine to get correct compressor inlet temperature
cmpTINmin = 500 + 273.15 ;
cmpTINmax = 700 + 273.15 ;
f1 = @(cmpTIN) heat_pump_function(SF1.state(1,1).T,SF1.state(1,5).T,cmpTIN,T0actual,0) ;  
[cmpTIN,~,~] = golden_search(f1,cmpTINmin,cmpTINmax,0.005,'Min') ;
% Rerun to get pretty pictures
[~, fldH, fldC, gas] = heat_pump_function(SF1.state(1,1).T,SF1.state(1,5).T,cmpTIN,T0actual,1) ; 

% Now need to balance heat input/output of hot store.
% Adjust mass flow rate of steam appropriately
% The molten salt is operating between the correct temps, and the heat pump
% mass flow is fixed by the input file. For this gas mdot, get a certain
% molten salt mdot. Discharge at the same mdot, which sets the steam mdot
mdot_des = mdot_des * fldH.state(1,1).mdot / SF1.state(1,5).mdot ;
mdot_tot = mdot_des ;

% Rerun steam function
cd '..\PTES_KE'
[~,~,~]	= Steam_fxn(P_cond) ; % Run again at design point
test_data
Wout_0 = CYC.WnetD ;
Qin_0 = CYC.Qin ;
eta_0 = CYC.eta ;
cd '..\PTES_code'

% Set up air and condensing steam. Find design air mass flow rate
% Inlet steam
WF1.state(1,5).p    = DEXP3.out.P ;
WF1.state(1,5).Q    = DEXP3.out.q ;
WF1.state(1,5).mdot = 1.0 * DEXP3.out.mdot ;
[WF1] = update(WF1,[1,5],3); 

% Outlet steam is fully condensed
WF1.state(1,6).p    = DEXP3.out.P ;
WF1.state(1,6).Q    = 0.0 ; 
WF1.state(1,6).mdot = 1.0 * DEXP3.out.mdot ;
[WF1] = update(WF1,[1,6],3); 

% Inlet air
AIR = fluid_class('Nitrogen','SF','CP','TTSE',1,30); % Not sure why air won't work
AIR.state(1,1).p    = 1.013e5 ;
AIR.state(1,1).T    = T0design ;
[AIR] = update(AIR,[1,1],1);

% Outlet air
AIR.state(1,2).p    = 1.013e5 ;
AIR.state(1,2).T    = WF1.state(1,5).T - dT1 ;
[AIR] = update(AIR,[1,2],1);

% Find air mass flow rate to achieve this
AIR.state(1,1).mdot = WF1.state(1,5).mdot * ...
                     (WF1.state(1,5).h - WF1.state(1,6).h) / ...
                     (AIR.state(1,2).h - AIR.state(1,1).h) ;

AIR.state(1,2).mdot = AIR.state(1,1).mdot ;

air_mdot0 = AIR.state(1,2).mdot ;

% Calculate the actual condenser pressure, which is at Tsat = T0 + change
% in temp of air flow at design point (assume it remains constant)
WF1.state(1,5).T    = T0actual + dT1 + AIR.state(1,2).T - AIR.state(1,1).T ;
WF1.state(1,5).Q    = 1.0 ;
[WF1] = update(WF1,[1,5],4);

P_cond = WF1.state(1,5).p / 1e5 ;

% Now recalculate steam performance for this condenser pressure
cd '..\PTES_KE'
[~,~,~]	= Steam_fxn(P_cond) ; % Run again at design point
test_data
Wout_nocooling = CYC.WnetD ;
Qin_nocooling = CYC.Qin ;
eta_nocooling = CYC.eta ;
cd '..\PTES_code'

% Recaculate air properties
% Set up air and condensing steam. Find design air mass flow rate
% Inlet steam
WF1.state(1,5).p    = DEXP3.out.P ;
WF1.state(1,5).Q    = DEXP3.out.q ;
WF1.state(1,5).mdot = 1.0 * DEXP3.out.mdot ;
[WF1] = update(WF1,[1,5],3); 

% Outlet steam is fully condensed
WF1.state(1,6).p    = DEXP3.out.P ;
WF1.state(1,6).Q    = 0.0 ; 
WF1.state(1,6).mdot = 1.0 * DEXP3.out.mdot ;
[WF1] = update(WF1,[1,6],3); 

% Inlet air
AIR = fluid_class('Nitrogen','SF','CP','TTSE',1,30); % Not sure why air won't work
AIR.state(1,1).p    = 1.013e5 ;
AIR.state(1,1).T    = T0actual ;
[AIR] = update(AIR,[1,1],1);

% Outlet air
AIR.state(1,2).p    = 1.013e5 ;
AIR.state(1,2).T    = WF1.state(1,5).T - dT1 ;
[AIR] = update(AIR,[1,2],1);

% Find air mass flow rate to achieve this
AIR.state(1,1).mdot = WF1.state(1,5).mdot * ...
                     (WF1.state(1,5).h - WF1.state(1,6).h) / ...
                     (AIR.state(1,2).h - AIR.state(1,1).h) ;

AIR.state(1,2).mdot = AIR.state(1,1).mdot ;

% Call golden search routine to get correct compressor inlet temperature
Pmin = P_cond*0.5 ;
Pmax = P_cond*0.99 ;

f2 = @(Pin) cooling_routine (AIR,WF1,fldC,Pin,eff,ploss,dT1) ;
[Pin,~,iter] = golden_search(f2,Pmin,Pmax,0.005,'Min') ;
[diff, fldCmdot, actmdot, AIR, WF1, SF2] = cooling_routine (AIR,WF1,fldC,Pin,eff,ploss,dT1) ;
Wout_cooling = CYC.WnetD ;
Qin_cooling = CYC.Qin ;
eta_cooling = CYC.eta ;




% ENDING
unloadlibrary coolprop % Not sure why I have to do this, but program crashes if I don't :(
cd '..'



% Function used to match cooling flows between storage and condenser
function [diff, fldCmdot, actmdot, AIR, WF1, SF2] = cooling_routine (AIR,WF1,fldC,P_cond,eff,ploss,dT1)

global degC BAR KIL MEG %unit conversion
degC = 273.15; BAR = 1.e5;	KIL = 1.e3;	MEG = 1.e6 ;
% Rerun steam function
cd '..\PTES_KE'
[~,~,~]	= Steam_fxn(P_cond) ;
test_data
cd '..\PTES_code'

% Reduce the condenser pressure - new steam inlet
WF1.state(1,7).p    = DEXP3.out.P;
WF1.state(1,7).Q    = DEXP3.out.q ; 
WF1.state(1,7).mdot = 1.0 * DEXP3.out.mdot ;
[WF1] = update(WF1,[1,7],3); 

% Air has same inlet properties and mass flow
AIR.state(1,3) = AIR.state(1,1) ;

% Air is still dT colder than exit
AIR.state(1,4).p    = 1.013e5 ;
AIR.state(1,4).T    = WF1.state(1,7).T - dT1 ;
[AIR] = update(AIR,[1,4],1);

% Now find steam outlet properties
WF1.state(1,8).p    = WF1.state(1,7).p ;
WF1.state(1,8).h    = WF1.state(1,7).h - AIR.state(1,3).mdot * ...
                      (AIR.state(1,4).h - AIR.state(1,3).h) / ...
                      WF1.state(1,7).mdot ;
WF1.state(1,8).mdot = 1.0 * DEXP3.out.mdot ;
[WF1] = update(WF1,[1,8],2); 

% The steam is now condensed using the cold storage fluid
% Can estimate cold fluid mass flow since we know the desired steam outlet
WF1.state(1,9).p    = WF1.state(1,8).p;
WF1.state(1,9).Q    = 0.0 ; 
WF1.state(1,9).mdot = WF1.state(1,8).mdot ;
[WF1] = update(WF1,[1,9],3); 

% Find mass flow rate of fluid that can condense the steam 
SF2 = fluid_class('INCOMP::MEG2[0.56]','SF','TAB',NaN,1,30); % working fluid
SF2.state(1,1).p    = 1.e5;
SF2.state(1,1).T    = fldC.state(1,2).T ;
SF2.state(1,2).p    = 1.e5;
SF2.state(1,2).T    = WF1.state(1,9).T - dT1 ;

% Update initial states
[SF2] = update(SF2,[1,1],1);
[SF2] = update(SF2,[1,2],1);

% Estimate cooling fluid mdot
SF2.state(1,1).mdot = WF1.state(1,9).mdot * ...
                     (WF1.state(1,8).h - WF1.state(1,9).h) / ...
                     (SF2.state(1,2).h - SF2.state(1,1).h) ;
SF2.state(1,2).mdot = SF2.state(1,1).mdot ;

% Run HEX routine
%[WF1,SF2,~,~] = hex_TQ_cond(WF1,[1,8],SF2,[1,1],eff,0,ploss,'hex',4,0);
%plot_hex(SF2,[1,1],WF1,[1,8],100,4)


fldCmdot = fldC.state(1,2).mdot ;
actmdot = SF2.state(1,2).mdot ;

diff = abs(fldCmdot - actmdot) ;


end



function [diff,fluidH,fluidC] = hex_function (mdot, Ttarg, mode, fluidH, indH, fluidC, indC, eff, Crat, ploss, stage_type, cond, var)


if mode == 0
    fluidH.state(indH(1),indH(2)).mdot = mdot ;
else
   error('Not implemented')
end

[fluidH, fluidC,~,~] = hex_TQ_cond(fluidH, indH, fluidC, indC, eff, Crat, ploss, stage_type, cond, var);

if mode == 0
    diff = abs(fluidC.state(indC(1),indC(2)+1).T - Ttarg) ;
else
   error('Not implemented') 
end


end