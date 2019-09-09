%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steam Rankine w/ Joule PTES cycle
%
% Josh McTigue & Kevin Ellingwood, 26 Jun, 2019
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
global rank_data FLD_rank m_dot_des m_dot_tot DEN DCN HP_sat LP_sat HP_Tin IP_Tin  ...
    DE_eta DC_eta eff_ex_mode TTD xf_des dP_rank dT_rank c_data COOL P_cond_cs dT_air ...
    DEXP1 SOL_reh DEXP2 DEXP3 DEXP4 dREJ DCMP1 DCMP2 DCMP3 DCMP4 FWH1 FWH2 ...
    SOL_pre SOL_gen SOL_sup FAN AIR_CS CYC_Rank P_cond_UB P_cond_LB %#ok<*NUSED>
% Joule cycle 
global  joule_data FLD_joule CCMP1 cHS cRCP CEXP1 CEXP2 CEXP3 cCS1 cCS2 cCS3 dP_joule dT_joule ... 
        T_HS_hi T_HS_lo T_amb CYC_Joule CCN CEN m_joule

%% Rankine cycle
% Set working directory
addpath     ('.\CoolProps\') ;
workspace   = 'C:\Users\jmctigue\Documents\solar-PTES\simple_concepts' ;
dir         = '.\inputs\Rank_PTES\' ;              	% Directory to input files
fname       = 'retrofit_3.txt' ; 	% Input file name

% Read in File Data and determine the setup of saturation profile
fIN         = fopen(strcat(dir,fname),'r');
readdum(fIN,3);                                         % This function just skips over lines that aren't of interest

% Fluid data
rank_data.name = char(readtxt(fIN,':',1,1)) ;               % Fluid name
rank_data.Tref = str2double(readtxt(fIN,':',1,1)) ;         % Ambient temperature
FLD_rank       = FLUID(rank_data);                       	% Call class constructor - for fluid class 
% Cycle data 
readdum(fIN,3); 
m_dot_des   = str2double(readtxt(fIN,':',1,1)) ;        % Design mass flow rate of cycle (kg/s)
m_dot_tot   = str2double(readtxt(fIN,':',1,1)) ;        % Actual mass flow rate of cycle (kg/s)
DEN         = str2double(readtxt(fIN,':',1,1)) ;        % Number of expansion stages
DCN         = DEN;                                      % Equal number of compression/pumping stages                            
HP_sat      = str2double(readtxt(fIN,':',1,1)) ;        % hi saturation pressure of cycle
LP_sat      = str2double(readtxt(fIN,':',1,1)) ;        % lo saturation pressure of cycle
HP_Tin      = str2double(readtxt(fIN,':',1,1)) + degC ; % superheat temperature, K
IP_Tin      = str2double(readtxt(fIN,':',1,1)) + degC ; % reheat temperature, K
DE_eta      = str2double(readtxt(fIN,':',1,1)) ;        % reference efficiency of expander
DC_eta      = str2double(readtxt(fIN,':',1,1)) ;        % reference efficiency of pumps
eff_ex_mode = str2double(readtxt(fIN,':',1,1)) ;        % method for determining off-design efficinecy for expanders
TTD     	= str2double(readtxt(fIN,':',1,1)) ;        % terminal T difference, K or C (pump inlets)
xf_des      = str2double(readtxt(fIN,',',1,2)) ;        % Design extraction fractions for FWH's
% HX losses
readdum(fIN,3); 
dP_rank  	= str2double(readtxt(fIN,':',1,1)) ;        % Fractional pressure loss in HX's
dT_rank 	= str2double(readtxt(fIN,':',1,1)) ;        % temperature difference between streams in HX's (°C)
% Variable Condensor information
readdum(fIN,3);
c_data.name = char(readtxt(fIN,':',1,1)) ;              % Cooling medium name
c_data.Tref = str2double(readtxt(fIN,':',1,1)) ;        % Ambient cooling temperature
COOL        = FLUID(c_data) ;
P_cond_cs   = str2double(readtxt(fIN,':',1,1)) ;        % Condenser operating pressure (cold side)
P_cond_UB   = str2double(readtxt(fIN,':',1,1)) ;        % Condenser UB pressure (steam side)
P_cond_LB   = str2double(readtxt(fIN,':',1,1)) ;        % Condenser LB pressure (steam side)
dT_air      = str2double(readtxt(fIN,':',1,1)) ;        % change in air temperature through condenser (°C)

%% Joule-Brayton cycle
readdum(fIN,3);                           % This function just skips over lines that aren't of interest

% Read in fluid data
joule_data.name = char(readtxt(fIN,':',1,1)) ;                      % Fluid name
joule_data.Tref = str2double(readtxt(fIN,':',1,1)) + degC ;         % Ambient temperature
FLD_joule = FLUID(joule_data) ;      % Call class constructor - for fluid class 

% Read in charging data
readdum(fIN,3); 
CCN             = str2double(readtxt(fIN,':',1,1)) ;        % Number of charging compressors
% Pre-allocate array size
CC              = struct('Pin',cell(1,CCN),'eta',cell(1,CCN),'Tin',cell(1,CCN),'beta',cell(1,CCN),'type',cell(1,CCN));
CC(1).Pin       = str2double(readtxt(fIN,':',1,1)) * BAR ;	% Inlet pressure to first compressor
Cpr             = str2double(readtxt(fIN,':',1,1)) ;        % Overall pressure ratio
LCpr            = strcmp(readtxt(fIN,':',1,1),'y') ;        % Logical - Set each compressor stage pressure ratio to be equal
if ~LCpr; Cpr = 1.0; end
for i = 1:CCN
   dat = str2double(readtxt(fIN,',',1,3));      % Reads various bits
   CC(i).eta = dat(1) ;                         % Isentropic efficiency
   CC(i).Tin = dat(2) + degC ;                  % Inlet temperature to compressor
   if LCpr
       CC(i).beta = Cpr ^ (1. / CCN) ;          % Pressure ratio if each pressure ratio is equa;
   else
       CC(i).beta = dat(3) ;                    % Pressure ratio if read in from file
       Cpr = CC(i).beta * Cpr ;
   end
   CC(i).type = "CMP" ;                         % Machine type is compressor
   
   % Class constructors for the charging compressors
   switch i
       case 1
           CCMP1 = MACHINE(CC(i)) ;
       case 2 
           CCMP2 = MACHINE(CC(i)) ;
       case 3
           CCMP3 = MACHINE(CC(i)) ;
       case 4
           CCMP4 = MACHINE(CC(i)) ;
       case 5
           error('Code not written for case where there are > 5 charging compressors');
   end
      
end
CEN = str2double(readtxt(fIN,':',1,1)) ;        % Number of charging expanders
% Pre-allocate array size
CE = struct('Pin',cell(1,CEN),'eta',cell(1,CEN),'Tin',cell(1,CEN),'beta',cell(1,CEN),'type',cell(1,CEN));
LCpr = strcmp(readtxt(fIN,':',1,1),'y') ;       % Logical - Set each expander stage pressure ratio to be equal
if ~LCpr ; Cpr = 1.0 ; end
for i = 1:CEN
   dat = str2double(readtxt(fIN,',',1,3));       % Reads various bits
   CE(i).eta  = dat(1) ;                         % Isentropic efficiency
   CE(i).Tin  = dat(2) + degC ;                  % Inlet temperature to expander
   CE(i).beta = dat(3) ;                         % Pressure ratio
   
   if LCpr
       CE(i).beta = Cpr ^ (1. / CEN) ;          % Pressure ratio if each pressure ratio is equal
   else
       CE(i).beta = dat(3) ;                    % Pressure ratio if read in from file
       Cpr = CE(i).beta * Cpr ;
   end
   CE(i).type = "EXP" ;                         % Machine type is compressor
   % Class constructors for the charging compressors
   switch i
       case 1
           CEXP1 = MACHINE(CE(i)) ;
       case 2 
           CEXP2 = MACHINE(CE(i)) ;
       case 3
           CEXP3 = MACHINE(CE(i)) ;
       case 4
           CEXP4 = MACHINE(CE(i)) ;
       case 5
           error('Code not written for case where there are > 4 charging expanders');
   end
end
% Read in losses in the heat exchangers
readdum(fIN,3); 
dP_joule    = str2double(readtxt(fIN,':',1,1)) ; % Fractional pressure loss
dT_joule    = str2double(readtxt(fIN,':',1,1)) ; % Temperature difference between streams in heat exchanger
% Storage and Condesner info
readdum(fIN,3);
T_HS_hi     = str2double(readtxt(fIN,':',1,1)) + degC; %HS high temperature
T_HS_lo     = str2double(readtxt(fIN,':',1,1)) + degC; %HS low  temperature 
T_amb       = str2double(readtxt(fIN,':',1,1)) + degC; %baseline ambient air temperature
% Close input file
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