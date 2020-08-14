clc
clear all

%dbstop if error

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

%Printresutls
%1 - Three cases for single material
%2 - Compare different materials for hybrid storage

Printresults=1;
  
% Set indices
iL = 1; i1 = 1; i2 = 1; ii=1;
m_hot_PB = zeros(3,1);
Vol_hot_PB = zeros(3,1);
m_hot_MS = zeros(3,1);
Vol_hot_MS = zeros(3,1);
V_hottank =zeros(3,1);
Tankcost = zeros(3,1);
mu_fluid = zeros(3,1);

%Bed = cell.empty( 100, 0 );

%Temperature assumptions of cold store
T_cold_initial =-50+273;
T_cold_final =10+273;

mdot=100;

%Cold Store properties for liquid storage and hybrid
F1 = fluidhybrid_class('Methanol','WF','CP','TTSE',1,5);
F1.statehybrid(iL,i1).p = 1e5;
F1.statehybrid(iL,i1).T = T_cold_initial; %in Kelvin
F1.statehybrid(iL,i1).mdot = mdot;

[F1] = update_hybrid(F1,[iL,i1],1);
state_cryogen = F1.statehybrid(i1,i2);

rho_cryogen=state_cryogen.rho;
Cp_cryogen = state_cryogen.Cp; % specific heat
mu_fluid(2)=3.26*0.001;
mu_fluid(3)=3.26*0.001;

%Cold Store properties for Packed bed only case
F2 = fluidhybrid_class('Argon','WF','CP','TTSE',1,5);
F2.statehybrid(iL,i1).p = 7e5;
F2.statehybrid(iL,i1).T = T_cold_initial; %in Kelvin
F2.statehybrid(iL,i1).mdot = mdot;

[F2] = update_hybrid(F2,[iL,i1],1);
state_cryogen_PB = F2.statehybrid(i1,i2);

rho_cryogen_PB=state_cryogen_PB.rho;
Cp_cryogen_PB = state_cryogen_PB.Cp; % specific heat
mu_fluid(1) = state_cryogen_PB.mu;

%Bed Properties
 AR=1; %Aspect Ratio
 phi = 0.4; %Porosity of bed assumed
%D = 1.64; %Diameter of the bed in m
dp = 0.05; %Diameter of the particle in m

%Limestone, Al2O3, Fe2O3, Fe3O4, SiO2, TiO2
Bed(1)="Fe3O4";
M=1; %index for the loop

if Printresults ==2
    Bed=["Limestone", "Al2O3", "Fe2O3", "Fe3O4", "SiO2", "TiO2"];
    M=6; %index for the loop
end

for nn=1:M

switch Bed(nn)
    case 'Limestone'
    Cp_bed=840; %Specific heat capacity of average rock
    k_particle=2.15; %Thermal conductivity of marble; k of rock varies from 0.4 to 7
    rho_particle = 2700; %Density of the particle
    
    case 'Al2O3'
    Cp_bed=540; %Specific heat capacity of average rock
    k_particle=56.8; %Thermal conductivity of marble; k of rock varies from 0.4 to 7
    rho_particle = 4040; %Density of the particle
    
    case 'Fe2O3'
    Cp_bed=500; %Specific heat capacity of average rock
    k_particle=6.95; %Thermal conductivity of marble; k of rock varies from 0.4 to 7
    rho_particle = 5242; %Density of the particle
    
    case 'Fe3O4'
    Cp_bed=520; %Specific heat capacity of average rock
    k_particle=3.94; %Thermal conductivity of marble; k of rock varies from 0.4 to 7
    rho_particle = 5175; %Density of the particle
    
    case 'SiO2'
    Cp_bed=570; %Specific heat capacity of average rock
    k_particle=22.5; %Thermal conductivity of marble; k of rock varies from 0.4 to 7
    rho_particle = 2650; %Density of the particle
    
    case 'TiO2'
    Cp_bed=550; %Specific heat capacity of average rock
    k_particle=11; %Thermal conductivity of marble; k of rock varies from 0.4 to 7
    rho_particle = 4230; %Density of the particle
    
end

%Heat requried to be stored in the TES system
Qcoldstore = 2894*10^6*3600; %Heat energy transfered to hot tank (from PTES code) in Joules
Qcoldstore = 863*10^6*3600;


%Volume of hot tank calculations
Qcoldstore_delT=Qcoldstore/(T_cold_final-T_cold_initial); %Qcoldstore / delT

%Volume of tank for a completely packed bed based TES system
V_coldtank(1)=Qcoldstore_delT/(rho_particle*(1-phi)*Cp_bed);

%Volume of tank for a completely cryogen salt based TES system
V_coldtank(2)=Qcoldstore_delT/(rho_cryogen*Cp_cryogen);

%Volume of tank for a hybrid TES system
V_coldtank(3)=Qcoldstore_delT/(rho_particle*(1-phi)*Cp_bed+rho_cryogen*phi*Cp_cryogen);

%Length and Diameter of Bed
Rad_bed=(V_coldtank./(AR.*2.*pi)).^(1./3);
D=2.*Rad_bed;
L=D.*AR;

%Packed bed cost correlation 

maxV  = 50e3 ;
Vol_PB = V_coldtank;
V_hybrid(nn) = V_coldtank(3);

for ii=1:3
    
if ii==1 %% Packed Bed only case
    Vol_cold_PB(1)=V_coldtank(1); %Volume of solid particles
    Vol_cold_MS(1)=0; %Volume of molten salt
    m_cold_PB(1) = rho_particle.*Vol_cold_PB(1);
    m_cold_MS(1) = rho_cryogen.*Vol_cold_MS(1);
    
        %Reynolds Number
     Area=pi.*(D(ii)./2).^2;
     Re=mdot.*dp/(Area.*mu_fluid(1));
    
    %Pressure drop Correlation
    fP(ii)=(150+3.89.*(Re./(1-phi)).^0.87).*((1-phi).^2./(phi.^3.*Re)); %Jones and Krier correlation
    delP_cold(ii)=fP(ii).*mdot.^2.*L(ii)./(rho_cryogen_PB.*dp.*Area.^2);
end

if ii==2 %Cryogen only case
     
     Vol_cold_PB(2)=0;
     Vol_cold_MS(2)=V_coldtank(2);
     m_cold_PB(2) = rho_particle.*Vol_cold_PB(2);
     m_cold_MS(2) = rho_cryogen.*Vol_cold_MS(2);
     
     %Reynolds Number
     Area=pi.*(D(ii)./2).^2;
     Re=mdot.*dp/(Area.*mu_fluid(1));
     
     fP(ii) = 0;
     delP_cold(ii) = 0;      
end
    
if ii==3 %cryogen and packed bed
Vol_cold_PB(3) = (1-phi)*V_coldtank(3);
Vol_cold_MS(3) = phi*V_coldtank(3); %Volume of cryogen in the tank
m_cold_PB(3) = rho_particle*Vol_cold_PB(3);
m_cold_MS(3) = rho_cryogen*Vol_cold_MS(3);

%Reynolds Number
     Area=pi.*(D(ii)./2).^2;
     Re=mdot.*dp/(Area.*mu_fluid(1));
    
    %Pressure drop Correlation
    fP(ii)=(150+3.89.*(Re./(1-phi)).^0.87).*((1-phi).^2./(phi.^3.*Re)); %Jones and Krier correlation
    delP_cold(ii)=fP(ii).*mdot.^2.*L(ii)./(rho_cryogen.*dp.*Area.^2);
end

if Vol_PB(ii) > maxV
   nA    = floor(Vol_PB(ii) / maxV); % Number of maxV tanks                       
   vA    = mod(Vol_PB(ii) / maxV,1)  ; % Remaining volume                      
   Tankcost(ii) = nA .* 13673 .* maxV.^0.557 ;                      
                         
   Tankcost(ii) = Tankcost(ii) + 13673 .* vA.^0.557 ;
else
   Tankcost(ii) = 13673 .* Vol_PB(ii).^0.557 ;
end                     
                                  
% Increase cost if pressurized
if ii==1 
    p=F2.statehybrid(iL,i1).p ./ 1e5 ;
end
if ii~=1
p = F1.statehybrid(iL,i1).p ./ 1e5 ;
end
                     
if p > 7
   Cfact = 0.922 + 0.0335.*p - 0.0003.*p.^2 +1e-6.*p.^3 ;
   Tankcost(ii) = Tankcost(ii) .* Cfact ;
end
CEind_curr=619.2;
CEind_2015=556.8;
Tankcost(ii) = Tankcost(ii) .* CEind_curr ./ CEind_2015; 

c_limestone_m3=36.96; %Limestone cost per m3
cost_limestone(ii) =c_limestone_m3*Vol_cold_PB(ii);

c_cryogen_kg = 0.8;
cost_cryogen(ii) = m_cold_MS(ii)*c_cryogen_kg;

ED(ii)=Qcoldstore/V_coldtank(ii);

end

ED(2)=ED(2)/2; %Taking 2 tanks into consideration
Tankcost(2)=2*Tankcost(2); %Taking 2 tanks into consideration
ED_hybrid(nn) = ED(3);
Tankcost_hybrid(nn) = Tankcost(3); %Total cost hybrid

P_drop_PB(nn) = delP_cold(1); %Storing Pressure drop values of Packed bed only case for different material(nn=1:6)
P_drop_MS(nn) = delP_cold(2); %Storing Pressure drop values (Zero)of Molten salt only case for different material(nn=1:6)
P_drop_hyb(nn) = delP_cold(3); %Storing Pressure drop values (Zero)of Molten salt only case for different material(nn=1:6)

end

if Printresults==1
    
fprintf('Storage based on Packed Bed\n')

fprintf('Porosity  HeatEnergy (MWh)  mass_limestone(kg)  mass_cryogen(kg)  Vol_Tank    Tank Cost   Limestone cost  cryogen salt cost  Totalcost  Energy Density\n');
fprintf('%6.2f %14.2f %18.2e %18.2e%16.2e %12.2e %13.2e %16.2e %14.2e %12.2e\n\n',phi,Qcoldstore/(3600*10^6)...
    ,m_cold_PB(1),m_cold_MS(1), V_coldtank(1) ,Tankcost(1),cost_limestone(1),cost_cryogen(1),Tankcost(1)+cost_limestone(1)+cost_cryogen(1),ED(1));

fprintf('Storage based on tanks completely filled with cryogen\n')

fprintf('Porosity  HeatEnergy (MWh)  mass_limestone(kg)  mass_cryogen(kg)  Vol_Tank    Tank Cost   Limestone cost  cryogen salt cost  Totalcost  Energy Density\n');
fprintf('%6.2f %14.2f %18.2e %18.2e%16.2e %12.2e %13.2e %16.2e %14.2e %12.2e\n\n',1,Qcoldstore/(3600*10^6)...
    ,m_cold_PB(2),m_cold_MS(2), V_coldtank(2) ,Tankcost(2),cost_limestone(2),cost_cryogen(2),Tankcost(2)+cost_limestone(2)+cost_cryogen(2),ED(2));

fprintf('Storage based on hybrid system \n')

fprintf('Porosity  HeatEnergy (MWh)  mass_limestone(kg)  mass_cryogen(kg)  Vol_Tank    Tank Cost   Limestone cost  cryogen salt cost  Totalcost  Energy Density\n');
fprintf('%6.2f %14.2f %18.2e %18.2e%16.2e %12.2e %13.2e %16.2e %14.2e %12.2e\n\n',phi,Qcoldstore/(3600*10^6)...
    ,m_cold_PB(3),m_cold_MS(3), V_coldtank(3) ,Tankcost(3),cost_limestone(3),cost_cryogen(3),Tankcost(3)+cost_limestone(3)+cost_cryogen(3),ED(3));

fprintf('Mode                  Aspect Ratio  Diameter  Length  Pressure Drop\n')
fprintf('Packed Bed only %14.2f %12.2f %8.2f %12.2e\n',AR, D(1), L(1), P_drop_PB(1));
fprintf('Hybrid %23.2f %12.2f %8.2f %12.2e\n',AR,D(1), L(1), P_drop_hyb(1));

end

if Printresults==2
    fprintf('Hybrid system - Cold Storage (Solid Material + Methanol)\n\n');
    fprintf(' Materials   Hybrid Energy Density\n')  
    fprintf(' Limestone   %4.2e\n Al2O3       %4.2e\n Fe2O3       %4.2e\n Fe3O4       %4.2e\n SiO2        %4.2e\n TiO2        %4.2e\n',...
        ED_hybrid(1),ED_hybrid(2),ED_hybrid(3),ED_hybrid(4),ED_hybrid(5),ED_hybrid(6))
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