%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steam Rankine w/ Joule PTES cycle
%
% Josh McTigue & Kevin Ellingwood, 7 Aug, 2019
% JoshuaDominic.McTigue@nrel.gov
% joshmctigue@hotmail.co.uk
%
% National Renewable Energy Laboratory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file acts as in the input file for the Rankine-JoulePTES cycle. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% Global Variables
% unit conversions
global degC BAR KIL MEG %unit conversion
degC = 273.15; BAR = 1.e5;	KIL = 1.e3;	MEG = 1.e6 ;
% Rankine cycle 
global data FLD mdot_des mdot_tot DEN DCN HP_sat LP_sat HP_Tin IP_Tin  ...
    DE_eta DC_eta eff_mode TTD xf_des dP dT %#ok<NUSED>
%% Set working directory
%workspace   = 'C:\Users\jmctigue\Documents\solar-PTES\simple_concepts' ;
dir         = '';                   % Directory to input files
fname       = 'test_steam.txt' ; 	% Input file name

% Read in File Data and determine the setup of saturation profile
fIN         = fopen(strcat(dir,fname),'r');
readdum(fIN,3);                                         % This function just skips over lines that aren't of interest

% Fluid data
data.name   = char(readtxt(fIN,':',1,1)) ;                % Fluid name
data.Tref   = str2double(readtxt(fIN,':',1,1)) ;          % Ambient temperature
FLD         = FLUID(data);                              % Call class constructor - for fluid class 
readdum(fIN,3); 
%% Cycle data 
mdot_des    = str2double(readtxt(fIN,':',1,1)) ;        % Design mass flow rate of cycle (kg/s)
mdot_tot    = str2double(readtxt(fIN,':',1,1)) ;        % Actual mass flow rate of cycle (kg/s)
DEN         = str2double(readtxt(fIN,':',1,1)) ;        % Number of expansion stages
DCN         = DEN;                                      % Equal number of compression/pumping stages                            
HP_sat      = str2double(readtxt(fIN,':',1,1)) ;        % high saturation pressure of cycle
LP_sat      = str2double(readtxt(fIN,':',1,1)) ;        % low saturation pressure of cycle
HP_Tin      = str2double(readtxt(fIN,':',1,1)) + degC ; % superheat temperature, K
IP_Tin      = str2double(readtxt(fIN,':',1,1)) + degC ; % reheat temperature, K
DE_eta      = str2double(readtxt(fIN,':',1,1)) ;        % reference efficiency of expander
DC_eta      = str2double(readtxt(fIN,':',1,1)) ;        % reference efficiency of pumps
eff_mode    = str2double(readtxt(fIN,':',1,1)) ;        % method for determining off-design efficinecy for expanders
TTD     	= str2double(readtxt(fIN,':',1,1)) ;        % terminal T difference, K (pump inlets)
xf_des      = str2double(readtxt(fIN,',',1,2)) ;        % Design extraction fractions for FWH's
readdum(fIN,3);
% HX  pressure losses 
dP          = str2double(readtxt(fIN,':',1,1)) ;        % Fractional pressure loss in HX's
fclose(fIN) ;

%% Support Fxns
%Short function that reads a few dummy lines from a file and dumps them
function readdum(fid,nlines)
    for i = 1:nlines
       fgetl(fid); 
    end
end
function txt = readtxt(fid,delim,col1,col2)
    txt = strtrim(strsplit(fgetl(fid),delim));
    txt = txt(col1:col2) ;
end