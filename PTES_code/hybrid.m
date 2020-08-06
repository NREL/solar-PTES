%%% To run this script, place it inside the PTES_code folder, alongside the
%%% PTES.m file

%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clc
clear;

% Enter debugging mode if an error occurs
dbstop if error
%dbclear all

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/','./Other/','./LIB/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\','.\Other\','.\LIB\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% INPUTS %%%

% Set indices
iL = 1; i1 = 1; i2 = 1; ii=1;

%Temperature range in packed bed
T_gas=581;
T_bed =445;
mdot=78;

% Declare fluids and specify initial conditions for packed bed regenerator
        % Nitrogen (hot, high pressure)
        F1 = fluidhybrid_class('Nitrogen','WF','CP','TTSE',1,5);
        F1.statehybrid(iL,i1).p = 50e5;
        F1.statehybrid(iL,i1).T = T_gas; %in Kelvin
        F1.statehybrid(iL,i1).mdot = mdot;
        
        % Nitrogen (hotstream, outlet)
        Fout = fluidhybrid_class('Nitrogen','WF','CP','TTSE',1,5);
        Fout.statehybrid(iL,i2).p = 50e5;
        Fout.statehybrid(iL,i2).T = T_bed;
        Fout.statehybrid(iL,i2).mdot = mdot;
        
        % Nitrogen (cold, low pressure)
        F2 = fluidhybrid_class('Nitrogen','WF','CP','TTSE',1,5);
        F2.statehybrid(iL,i2).p = 10e5;
        F2.statehybrid(iL,i2).T = T_bed;
        F2.statehybrid(iL,i2).mdot = mdot;

% Update fluid states
[F1] = update_hybrid(F1,[iL,i1],1);
[F2] = update_hybrid(F2,[iL,i2],1);
[Fout] = update_hybrid(Fout,[iL,i2],1);

state_hot = F1.statehybrid(i1,i2);
state_cold = F2.statehybrid(i1,i2);
state_out = Fout.statehybrid(i1,i2);

%Packed bed Design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1=state_hot.h;
h2=state_out.h;
fprintf('                         Hybrid system\n\n');

%Gas Properties hot stram
rho_g_hot=state_hot.rho;
mu_g_hot = state_hot.mu; %Dynamic viscosity of hot nitrogen 50 bar
k_gas_hot = state_hot.k_therm; %thermal conductivity of nitrogen
Cp_hot = state_hot.Cp; %Hot Nitrogen specific heat

%Gas Properties cold stream
rho_g_cold = state_cold.rho;
mu_g_cold = state_cold.mu;
Cp_cold=state_cold.Cp;
k_gas_cold =state_cold.k_therm;


%Bed Properties
Cp_bed=840; %Specific heat capacity of average rock
k_particle=2.15; %Thermal conductivity of marble; k of rock varies from 0.4 to 7
dp = 0.05; %Diameter of the particle in m
rho_particle = 2700; %Density of the particle
phi = 0.4; %Porosity of bed assumed
D = 1.64; %Diameter of the bed in m


%Packed bed Regenerator assumptions
Ph=35; % Assumed Hot cycle duration (in seconds)
Pc=35; % Assumed Cold cycle duration
cycle_time=Ph+Pc;  %Total time period 

%Heat transfer rate in the recuperator section
%h1=0.607*10^6; %h value of gas stream at the inlet of the recuperator from exit of hot storage charge cycle
%h2=0.316*10^6; %h value of gas stream at the outlet of the recuperator 
Qrate=mdot*(h1-h2);

%Heat energy that needs to be stored in the regenerator

Qstore=Qrate*Ph; 
Qstore_kwh=Qstore/3600000;

%Estimating mass of packed bed required
T_bed_initial=445;
T_bed_final = 573;

mbed=1.5*Qstore/(Cp_bed*(T_bed_final-T_bed_initial));

%Calculating L from the mass of the bed
Vol_s=mbed/rho_particle;
V_bed = Vol_s/(1-phi);
L = V_bed/(pi*(D/2)^2);

%Estimating equivalent matrix capacity ratio 

C_max=mdot*Cp_hot;
C_min = mdot*Cp_cold;
CR=C_min / C_max;

C_m=mbed*Cp_bed/(cycle_time*C_min);
C_me = 2*CR*C_m/(1+CR);

%Nusselet number packed bed hot stream
%Reynold's number calculated assuming the pipe diameter is same as packed bed diameter

Area=pi*(D/2)^2;
Re=mdot*dp/(Area*mu_g_hot);
Pr=Cp_hot*mu_g_hot/k_gas_hot;
Nu_PB=2.19*Pr^(1/3)*Re^(1/3)+0.78*Pr^(1/3)*Re^0.62;

%Nusselet number packed bed cold stream
Re2=mdot*dp/(Area*mu_g_cold);
Pr2=Cp_cold*mu_g_cold/k_gas_cold;

Nu_PB2=2.19*Pr2^(1/3)*Re2^(1/3)+0.78*Pr2^(1/3)*Re2^0.62;

%Heat transfer coefficients
h_hot= Nu_PB*k_gas_hot/dp;
h_cold = Nu_PB2*k_gas_cold/dp;
 

%Calculating heat transfer area
Nr =Vol_s/(4/3*pi*(dp/2)^3);
A = Nr*4*pi*(dp/2)^2; %Surface area of sphere

%Estimating NTU
UA = ((1/(h_hot*A))+(1/(h_cold*A)))^-1;
NTU = UA/C_min;
NTUe=2*CR*NTU/(1+CR); %equivalent NTU

%After lookup to get priliminary effectives, net effectiveness of the
%regenerator can be estimated

epsilon = 0.7; %From lookup table
chi_re=(1-CR^2)*(epsilon/(1-epsilon))/(2*CR);
eff_re = (1 - exp(-chi_re))/(1-CR*exp(-chi_re));

%Pressure drop estimation hot stream
fP=(150+3.89*(Re/(1-phi))^0.87)*((1-phi)^2/(phi^3*Re)); %Jones and Krier correlation
delP_hot=fP*mdot^2*L/(rho_g_hot*dp*Area^2);

%Pressure drop estimation cold stream
fP_c=(150+3.89*(Re2/(1-phi))^0.87)*((1-phi)^2/(phi^3*Re2));
delP_cold=fP_c*mdot^2*L/(rho_g_cold*dp*Area^2);

%Packed bed cost correlation 
maxV  = 50e3 ;
Vol_PB = Area * L;
if Vol_PB > maxV
   nA    = floor(Vol_PB / maxV); % Number of maxV tanks                       
   vA    = mod(Vol_PB / maxV,1)  ; % Remaining volume                      
   PBcost = nA * 13673 * maxV^0.557 ;                      
                         
   PBcost = PBcost + 13673 * vA^0.557 ;
else
   PBcost = 13673 * Vol_PB^0.557 ;
end                     
                                  
% Increase cost if pressurized
       
p = max(delP_hot, delP_cold) / 1e5 ;
                     
if p > 7
   Cfact = 0.922 + 0.0335*p - 0.0003*p^2 +1e-6*p^3 ;
   PBcost = PBcost * Cfact ;
end
CEind_curr=619.2;
CEind_2015=556.8;
PBcost = PBcost * CEind_curr / CEind_2015; 

%Thickness of packed bed shell calculation

Pressure = 725; %Allowable pressure in psi (Correlation Taken from engineering tool box)
Dia=D*39.3701; %Outsie pipe diameter in inches
Y=0.4; %Wall thickness coefficient
S=16000; %Allowable tensile strenght psi
E=0.8; %Quality factor
th = (Pressure*Dia/(2*(Pressure*Y+S*E)))/39.3701; %thickness in m

c_limestone_m3=36.96; %Limestone cost per m3
cost_limestone =c_limestone_m3*Vol_s;

TC = PBcost+cost_limestone; %Total material cost

PB=PB_class(L, D, dp, Ph, Pc, Nu_PB, Nu_PB2, C_me, NTUe, eff_re, delP_hot,delP_cold, TC);
PB.Length = L;
PB.Diameter = D;
PB.dp=dp;
PB.hot_time=Ph;
PB.cold_time=Pc;
PB.Nu_hotstream=Nu_PB;
PB.Nu_coldstream=Nu_PB2;
PB.CM=C_me;
PB.NTU_hotstream=NTUe;
PB.eff=eff_re;
PB.pressure_drop_hot=delP_hot;
PB.pressure_drop_cold=delP_cold;
PB.total_cost=TC;

 fprintf('Packed bed summary: T_in = %4.2f to T_intermediate = %4.2f \n', T_gas, T_bed)
 fprintf('Length  Diameter  d_particle  hot_time  cold_time  Nu_hotstream  Nu_coldstream  CM  NTU_hotstream  eff  PB.pressure_drop_hot  PB.pressure_drop_cold  PB.total_cost\n')
 fprintf('%6.2f %7.2f  %10.4f %9.2f %11.2f %11.3f %12.3f %9.2f %10.3f %8.2f %16.4f %18.4f %19.1e \n',...
     PB.Length,PB.Diameter, PB.dp,PB.hot_time,PB.cold_time,PB.Nu_hotstream,PB.Nu_coldstream,PB.CM,PB.NTU_hotstream,PB.eff, PB.pressure_drop_hot/1e5, PB.pressure_drop_cold/1e5, PB.total_cost);

%Heat exchanger design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        % Set hex_mode
        hex_mode = 0;
        par = 0;
        
        F1 = fluid_class('Nitrogen','WF','CP','TTSE',1,5);
        F1.state(iL,i1).p = 50e5;
        F1.state(iL,i1).T = T_bed; %in Kelvin
        F1.state(iL,i1).mdot = mdot;
        
        % Helium (cold, low pressure)
        F2 = fluid_class('Nitrogen','WF','CP','TTSE',1,5);
        F2.state(iL,i2).p = 10e5;
        F2.state(iL,i2).T = 315;
        F2.state(iL,i2).mdot = mdot;
                
% Update fluid states
[F1] = update(F1,[iL,i1],1);
[F2] = update(F2,[iL,i2],1);

state_hot = F1.state(i1,i2);
state_cold = F2.state(i1,i2);

% Specify HX settings
NX   = 100; % Number of sections (grid)
name = 'hot';
stage_type = 'hex';
model = 'geom';
CEind = create_CEindex() ;

costmode = 23;
switch model
    case 'eff'
        eff   = 0.97;
        ploss = 0.01;
        par1  = eff;
        par2  = ploss;
        par3  = [];
        par4  = [];
        
    case 'UA'
        UA    = 1e6;
        ploss = 0.01;
        par1  = UA;
        par2  = ploss;
        par3  = [];
        par4  = [];
                
    case 'geom'
        % Specify HEX geometry based on performance objectives
                % For comparisson with analytical results
                eff   = 0.97;
                ploss = 0.01;
                D1    = 2e-2; %Hydraulic diameter
                shape = 'circular'; 
        
        par1  = eff;
        par2  = ploss;
        par3  = D1;
        par4  = shape;
end

%%% CONSTRUCT HEAT EXCHANGER %%%
HX = hx_class(name, stage_type, costmode, NX, 1, 1, model, par1, par2, par3, par4);

% Specify geometry manually

%%% DESIGN PERFORMANCE %%%

% Run heat exchanger model under design conditions
[HX,~,~,~,~] = hex_func(HX,iL,F1,i1,F2,i2,hex_mode,par);

cap_cost = 0 ;
cap_sens = 0;

 HX = HX_cost(HX, CEind) ;
    cap_cost = cap_cost + HX.hx_cost.COST ;
    cap_sens = cap_sens + cost_sens(HX.hx_cost, 1); 
%% Print
fprintf('\n\n');
if strcmp(HX.model,'geom')    
    fprintf('Heat Exchanger Summary: T_intermediate = %4.2f to T_out = %4.2f\n',F1.state(iL,i1).T,F2.state(iL,i1).T);  
    print_hexs(HX,1);
end
 
function CEindex = create_CEindex()

    CEindex = zeros(2019,1) ;
    CEindex(1947) = 	64.8  ;
    CEindex(1948) = 	70.2  ;
    CEindex(1949) = 	71.4  ;
    CEindex(1950) = 	73.9  ;
    CEindex(1951) = 	80.4  ;
    CEindex(1952) = 	81.3  ;
    CEindex(1953) = 	84.7  ;
    CEindex(1954) = 	86.1  ;
    CEindex(1955) = 	88.3  ;
    CEindex(1956) = 	93.9  ;
    CEindex(1957) = 	98.5  ;
    CEindex(1958) = 	99.7  ;
    CEindex(1959) = 	101.8 ;
    CEindex(1960) = 	102	;
    CEindex(1961) = 	101.5 ;
    CEindex(1962) = 	102 ;
    CEindex(1963) = 	102.4 ;
    CEindex(1964) = 	103.3 ;
    CEindex(1965) = 	104.2 ;
    CEindex(1966) = 	107.2 ;
    CEindex(1967) = 	109.7 ;
    CEindex(1968) = 	113.7 ;
    CEindex(1969) = 	119 ;
    CEindex(1970) = 	125.7 ;
    CEindex(1971) = 	132.2 ;
    CEindex(1972) = 	137.2 ;
    CEindex(1973) = 	144.1 ;
    CEindex(1974) = 	165.4 ;
    CEindex(1975) = 	182.3 ;
    CEindex(1976) = 	192 ;
    CEindex(1977) = 	204.1 ;
    CEindex(1978) = 	218.8 ;
    CEindex(1979) = 	238.7 ;
    CEindex(1980) = 	261.1 ;
    CEindex(1981) = 	297 ;
    CEindex(1982) = 	314 ;
    CEindex(1983) = 	317 ;
    CEindex(1984) = 	322.6 ;
    CEindex(1985) = 	325.3 ;
    CEindex(1986) = 	318.3 ;
    CEindex(1987) = 	323.7 ;
    CEindex(1988) = 	342.4 ;
    CEindex(1989) = 	355.5 ;
    CEindex(1990) = 	357.6 ;
    CEindex(1991) = 	361.3 ;
    CEindex(1992) = 	358.2 ;
    CEindex(1993) = 	359.2 ;
    CEindex(1994) = 	368.1 ;
    CEindex(1995) = 	381.1 ;
    CEindex(1996) = 	381.7 ;
    CEindex(1997) = 	386.5 ;
    CEindex(1998) = 	389.5 ;
    CEindex(1999) = 	390.6 ;
    CEindex(2000) = 	394.1 ;
    CEindex(2001) = 	394.3 ;
    CEindex(2002) = 	395.6 ;
    CEindex(2003) = 	402   ;
    CEindex(2004) = 	444.2 ;
    CEindex(2005) = 	468.2 ;
    CEindex(2006) = 	499.6 ;
    CEindex(2007) = 	525.4 ;
    CEindex(2008) = 	575.4 ;
    CEindex(2009) = 	521.9 ;
    CEindex(2010) = 	550.8 ;
    CEindex(2011) = 	585.7 ;
    CEindex(2012) = 	584.6 ;
    CEindex(2013) = 	567.3 ;
    CEindex(2014) = 	576.1 ;
    CEindex(2015) = 	556.8 ;
    CEindex(2016) = 	541.7 ;
    CEindex(2017) = 	567.5 ;
    CEindex(2018) = 	603.1 ;
    CEindex(2019) = 	619.2 ;

end
