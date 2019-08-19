% Declare the state points of the gas for the PTES cycle
% 'WF' is working fluid
% 'CP' or 'TAB' indicate CoolProp or Tabular reading modes
% Num indicates number of preallocated streams
gas1  = fluid_class(gasName,'WF','CP','HEOS',30); % main cycle
gas2  = fluid_class(gasName,'WF','CP','HEOS',30); % secondary stream 1
gas3  = fluid_class(gasName,'WF','CP','HEOS',30); % secondary stream 2
gasHP = fluid_class(gasName,'WF','CP','HEOS',30); % heat pump cycle

% Declare the streams of storage fluids and the storage tanks
% There is one stream of storage fluid for each compression and expansion
% process during charge. Each stream a 2-by-2 matrix of states (rows for
% charge/discharge, columns for inlet/outlet)
% Double tanks have four states:
% 1=begin charge, 2=end charge, 3=start discharge, 4=end discharge
% Hot fluids and tanks
fluidH(1:Nc_ch) = fluid_class(fluidHname,'SF','TAB',0,2);
HT = double_tank_class(4);  %hot double tank with 4 states
% Cold fluids and tanks
fluidC1(1:1) = fluid_class(fluidC1name,'SF','TAB',0,2);
CT1 = double_tank_class(4); %cold double tank with 4 states
fluidC2(1:1) = fluid_class(fluidC2name,'SF','TAB',0,2);
CT2 = double_tank_class(4); %cold double tank with 4 states
% Air sources and sinks
LA(1:4) = tank_state_class;  %liquid air tank
AA(1:4) = tank_state_class;  %ambient air source/sink

% Declare heat rejection streams
environ = environment_class(T0,p0,10);

% Set number of points to plot each stage
num = 10;