clc
clear all

mdot=78; %Mass flow rate in kg/s

%Gas Properties hot stram
mu_g_hot = 29.331*10^-6; %Dynamic viscosity of hot nitrogen
k_gas_hot = 0.045; %thermal conductivity of nitrogen
rho_g_hot = 27.48; %Density of nitrogen at a given P and T
Cp_hot = 1070; %Hot Nitrogen specific heat

%Gas Properties cold stream
mu_g_cold = 19*10^-6;
Cp_cold = 1040; %Cold Nitrogen specific heat
k_gas_cold =0.028;
rho_g_cold = 11.23;


%Bed Properties
Cp_bed=840; %Specific heat capacity of average rock
k_particle=2.15; %Thermal conductivity of marble; k of rock varies from 0.4 to 7
dp = 0.003; %Diameter of the particle in m
rho_particle = 2700; %Density of the particle
phi = 0.39; %Porosity of bed assumed
D = 2; %Diameter of the bed in m


%Packed bed Regenerator assumptions
Ph=80; % Assumed Hot cycle duration (in seconds)
Pc=80; % Assumed Cold cycle duration
cycle_time=Ph+Pc;  %Total time period 


%Heat transfer rate in the recuperator section
h1=0.687*10^6; %h value of gas stream at the inlet of the recuperator from exit of hot storage charge cycle
h2=0.318*10^6; %h value of gas stream at the outlet of the recuperator 
Qrate=mdot*(h1-h2)/2;


fprintf('Mass flow rate = %4.2f kg/s\n',mdot);
fprintf('Enthaly change of the hot stream in the recuperator obtained from PTES charge cycle= %4.2f J/kg \n',h1-h2);

%Heat energy that needs to be stored in the regenerator

Qstore=Qrate*Ph; 
Qstore_kwh=Qstore/3600000;

fprintf('Total heat energy transfered into the bed from the hot stream for %4.2f seconds (T_gas_in = 581, T_bed_initial = 315 K) = %4.2f kWh\n',Ph, Qstore_kwh);

%Estimating mass of packed bed required

T_gas=581;
T_bed =315;
mbed=1.5*Qstore/(Cp_bed*(T_gas-T_bed));

fprintf('Mass of Limestone particles of size %4.2f mm required (multiplied by 1.5) = %4.2f kg \n\n',dp*1000, mbed);

%Calculating Geometry of packed bed
%AR=3;
%V_bed = Vol_s/(1-phi);
%Radius=nthroot((V_bed.*2./(AR.*pi)),3);
%D=2*Radius;
%L=AR.*D
%D=2.52;
Vol_s=mbed/rho_particle;

V_bed = Vol_s/(1-phi);
L = V_bed/(pi*(D/2)^2);
%Estimating equivalent matrix capacity ratio 

C_max=mdot*Cp_hot;
C_min = mdot*Cp_cold;
CR=C_min / C_max;

C_m=mbed*Cp_bed/(cycle_time*C_min);
C_me = 2*CR*C_m/(1+CR);
fprintf('Equivalent matrix capacity ratio = %4.2f \n', C_me);

%Nusselet number packed bed hot stream
%Reynold's number calculated assuming the pipe diameter is same as packed bed diameter

Area=pi*(D/2)^2;
Re=mdot*dp/(Area*mu_g_hot);
Pr=Cp_hot*mu_g_hot/k_gas_hot;

%Nu_l=0.664*Pr^(1/3)*(Re/phi)^0.5; %Gnielinksi correlation
%Nu_t= 0.037.*(Re/phi).^0.8.*Pr / (1+2.433.*(Re/phi).^(-0.1).*(Pr.^(2/3)-1));
%Nu_sp=2 + (Nu_t+Nu_l)^0.5;
%fe = 1+1.5*(1-phi);
%Nu_PB=fe*Nu_sp
Nu_PB=2.19*Pr^(1/3)*Re^(1/3)+0.78*Pr^(1/3)*Re^0.62;

%Nusselet number packed bed cold stream
Re2=mdot*dp/(Area*mu_g_cold);
Pr2=Cp_cold*mu_g_cold/k_gas_cold;

%Nu_l2=0.664*Pr2^(1/3)*(Re2/phi)^0.5;
%Nu_t2= 0.037.*(Re2/phi).^0.8.*Pr2 / (1+2.433.*(Re2/phi).^(-0.1).*(Pr2.^(2/3)-1));
%Nu_sp2=2 + (Nu_t2+Nu_l2)^0.5;
%fe2 = 1+1.5*(1-phi);
%Nu_PB2=fe2*Nu_sp2;

Nu_PB2=2.19*Pr2^(1/3)*Re2^(1/3)+0.78*Pr2^(1/3)*Re2^0.62;


%Heat transfer coefficients
h_hot= Nu_PB*k_gas_hot/dp;
h_cold = Nu_PB2*k_gas_cold/dp;

fprintf('Estimated Nusselt number of hot stream= %4.2f \n', Nu_PB);
fprintf('Estimated Nusselt number of cold stream= %4.2f\n', Nu_PB2);


%Calculating heat transfer area
Nr =Vol_s/(4/3*pi*(dp/2)^3);
A = Nr*pi*(dp/2)^2;

%Estimating NTU
UA = ((1/(h_hot*A))+(1/(h_cold*A)))^-1;
NTU = UA/C_min;
NTUe=2*CR*NTU/(1+CR); %equivalent NTU
fprintf('Estimated equivalent NTU = %4.2f \n', NTUe);

%After lookup to get priliminary effectives, net effectiveness of the
%regenerator can be estimated

epsilon = 0.79; %From lookup table
chi_re=(1-CR^2)*(epsilon/(1-epsilon))/(2*CR);
eff_re = (1 - exp(-chi_re))/(1-CR*exp(-chi_re));

fprintf('Effectiveness of the regenerator (look up table not yet coded) = %4.2f \n\n', eff_re);

%Pressure drop estimation hot stream
%Re_m=Re/(1-phi);
%f1L=136/(1-phi)^0.38;
%f1T=29/((1-phi)^1.45*phi^2);
%f2=1.87*phi^0.75/(1-phi)^0.26;
%qq=exp(-1*(phi^2*(1-phi)*Re_m)/12.6);
%fP=((qq*f1L/Re_m)+(1-qq)*(f2+(f1T/Re_m)))*(1-phi)/phi^3;
%delP_hot=fP*mdot^2*L/(rho_g_hot*dp*Area^2);
fP=(150+3.89*(Re/(1-phi))^0.87)*((1-phi)^2/(phi^3*Re)); %Jones and Krier correlation
delP_hot=fP*mdot^2*L/(rho_g_hot*dp*Area^2);

%Pressure drop estimation cold stream
%Re_m2=Re2/(1-phi);
%f1L_c=136/(1-phi)^0.38;
%f1T_c=29/((1-phi)^1.45*phi^2);
%f2_c=1.87*phi^0.75/(1-phi)^0.26;
%qq_c=exp(-(phi^2*(1-phi)*Re_m2)/12.6);
%fP_c=((qq_c*f1L_c/Re_m2)+(1-qq_c)*(f2_c+(f1T_c/Re_m2)))*(1-phi)/phi^3;
fP_c=(150+3.89*(Re2/(1-phi))^0.87)*((1-phi)^2/(phi^3*Re2));
delP_cold=fP_c*mdot^2*L/(rho_g_cold*dp*Area^2);

fprintf('Hot stream pressure drop = %4.2f bar\n',delP_hot/10^5);
fprintf('Cold stream pressure drop = %4.2f bar\n',delP_cold/10^5);

% Insulation thickness and volume
%tinsA = AR * (exp(obj.costdat.ins_k / (AR * UA)) - 1) ;
%obj.ins_volA = pi * tinsA * (2 + tinsA) * AL ; % Side wall insulation

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

%Material Cost correlations

Pressure = 725; %Allowable pressure in psi (Correlation Taken from engineering tool box)
Dia=D*39.3701; %Outsie pipe diameter in inches
Y=0.4; %Wall thickness coefficient
S=16000; %Allowable tensile strenght psi
E=0.8; %Quality factor
th = (Pressure*Dia/(2*(Pressure*Y+S*E)))/39.3701; %thickness in m
rho_metal=7850;
weight=(pi*((D/2)+th)^2*L-pi*(D/2)^2*L)*rho_metal*2.205; %Shell weight in pounds

cost_ss=575*weight^0.609; %Stainless steel cost
c_limestone_m3=36.96; %Limestone cost per m3
cost_limestone =c_limestone_m3*Vol_s;

TC = PBcost+cost_limestone; %Total material cost


fprintf('Length of bed calculated for a diameter of %4.2f m is = %4.2f m \n', D,L);
fprintf('Calculated thickness of stainless for an allowable pressure of %4.2f bar = %4.2f m \n',Pressure/14.504,th);
fprintf('Total cost of two packed bed based periodic regenerators(tank cost + Limestone) = %4.2f USD\n\n',2*TC);

%Heat changer based recuperator
T_gas_out=315;
delT1=T_gas-T_gas_out;
delT2=573-300; %State points 7 and 8 from PTES charge cycle
LMTD=(delT1-delT2)/log(delT1/delT2);
UA_hx=Qrate/LMTD;

Area_HX = 1.76*10^4; %Area of heat exchanger (recuperator) obtained from PTES code
Len_HX = 15.5; %Length of heat exchanger (recuperator) obtained from PTES code
Vol_HX = Area_HX*Len_HX;
weight_HX = rho_metal*Vol_HX*2.205;

cost_ss_HX=575*weight_HX^0.609; %Stainless steel cost
recuperator_cost=2.45e7;
fprintf('Recuperator parameters from PTES cycle, Area: %4.2f Length: %4.2f \n',Area_HX,Len_HX);
fprintf('Recuperator cost from PTES cycle = %4.2f USD\n',recuperator_cost);

%Calculation of viscosity

M_W= 14.0067; %Molecular weight of nitrogen
T_R(1) = T_gas*1.8; %Convert from Kelvin to Rankine units
T_R(2) = T_bed*1.8;
rho_g_gcc(1)=rho_g_hot/1000;
rho_g_gcc(2)=rho_g_cold/1000;

for kk=1:2
X(kk)=-3.41179 + (999.9949/T_R(kk))+0.14457226* M_W;
Y(kk)=1.293365+0.11245*X(kk);
K(kk)=(0.0001*(9.206483+3.550036*M_W)*T_R(kk)^1.217156)/(208.9852+18.70586*M_W+T_R(kk));
mu_g(kk)=K(kk)*exp(X(kk)*rho_g_gcc(kk)^Y(kk))*0.001; %Viscosity in SI, multiplied by 0.001
end

mu_g_hot=mu_g(1);
mu_g_cold=mu_g(2);