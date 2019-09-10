%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steam Rankine test script
%
% This file is a test script which calls on a steam cycle function. The
% function inputs the condenser pressure of the steam cycle and outputs the
% parameters needed to simulate charging of the Joule heat pump. Other
% potential inputs are defined in the text input file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine Operating System
c = computer();

% Addpaths and load CoolProp
switch computer
    case 'GLNXA64' %Linux
        addpath('./_inputs_/','./_classes_/','./_generic_/');
    case 'PCWIN64' %Windows
        addpath('.\_inputs_\','.\_classes_\','.\_generic_\');   
end
load_coolprop
set_graphics

% Input file reads data from input text file (test_steam.txt) and defines
% global variables used in steam cycle function
test_input

% The function outputs the flowrates for each turbine, the inlet
% temperature of the steam reheater, and the inlet temperature of the
% boiler feedwater preheater. 
P_cond = 0.2051;
[mdot,reh_Tin,pre_Tin]	= Steam_fxn(P_cond);

% Test_data file allocates cycle points in a single array DIS. Cycle
% metrics are allocated in CYC using the Cycle class
test_data

fprintf(1,'\nConditions at Reheater inlet\n');
fprintf(1,'T = %8.2f , p = %8.2f bar',REH.Cin.T - degC,REH.Cin.P/BAR);
fprintf(1,'\nConditions at Reheater outlet\n');
fprintf(1,'T = %8.2f C, p = %8.2f bar',REH.Cout.T - degC,REH.Cout.P/BAR);
