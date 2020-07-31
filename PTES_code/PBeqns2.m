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
Ph=40:40:400; % Assumed Hot cycle duration (in seconds)
Pc=40:40:400; % Assumed Cold cycle duration
cycle_time=Ph+Pc;  %Total time period 
N=10;
AR=zeros(1,N);


%Heat transfer rate in the recuperator section
h1=0.687*10^6; %h value of gas stream at the inlet of the recuperator from exit of hot storage charge cycle
h2=0.318*10^6; %h value of gas stream at the outlet of the recuperator 
Qrate=mdot*(h1-h2);


%Heat energy that needs to be stored in the regenerator

for j=1:N
Qstore=Qrate.*Ph(j); 
Qstore_kwh=Qstore./3600000;

%Estimating mass of packed bed required

T_gas=581;
T_bed =315;
mbed(j)=1.5.*Qstore./(Cp_bed.*(T_gas-T_bed));
%Calculating Geometry of packed bed

Vol_s=mbed(j)./rho_particle;
V_bed = Vol_s./(1-phi);
L = V_bed./(pi.*(D./2).^2);
AR(j)=L./D;
%Estimating equivalent matrix capacity ratio 

C_max=mdot.*Cp_hot;
C_min = mdot.*Cp_cold;
CR=C_min ./ C_max;

C_m=mbed(j).*Cp_bed./(cycle_time.*C_min);
C_me = 2.*CR.*C_m./(1+CR);


%Nusselet number packed bed hot stream
%Reynold's number calculated assuming the pipe diameter is same as packed bed diameter

Area=pi.*(D/2).^2;
Re=mdot.*dp/(Area.*mu_g_hot);
Pr=Cp_hot.*mu_g_hot/k_gas_hot;

%Nu_l=0.664*Pr^(1/3)*(Re/phi)^0.5; %Gnielinksi correlation
%Nu_t= 0.037.*(Re/phi).^0.8.*Pr / (1+2.433.*(Re/phi).^(-0.1).*(Pr.^(2/3)-1));
%Nu_sp=2 + (Nu_t+Nu_l)^0.5;
%fe = 1+1.5*(1-phi);
%Nu_PB=fe*Nu_sp
Nu_PB=2.19.*Pr^(1/3).*Re.^(1./3)+0.78.*Pr.^(1./3).*Re.^0.62;

%Nusselet number packed bed cold stream
Re2=mdot.*dp./(Area.*mu_g_cold);
Pr2=Cp_cold.*mu_g_cold./k_gas_cold;

%Nu_l2=0.664*Pr2^(1/3)*(Re2/phi)^0.5;
%Nu_t2= 0.037.*(Re2/phi).^0.8.*Pr2 / (1+2.433.*(Re2/phi).^(-0.1).*(Pr2.^(2/3)-1));
%Nu_sp2=2 + (Nu_t2+Nu_l2)^0.5;
%fe2 = 1+1.5*(1-phi);
%Nu_PB2=fe2*Nu_sp2;

Nu_PB2=2.19.*Pr2.^(1./3).*Re2.^(1./3)+0.78.*Pr2.^(1./3).*Re2.^0.62;


%Heat transfer coefficients
h_hot= Nu_PB.*k_gas_hot./dp;
h_cold = Nu_PB2.*k_gas_cold./dp;


%Calculating heat transfer area
Nr =Vol_s./(4./3.*pi.*(dp./2).^3);
A = Nr.*pi.*(dp./2).^2;

%Estimating NTU
UA = ((1./(h_hot.*A))+(1./(h_cold.*A))).^-1;
NTU = UA./C_min;
NTUe(j)=2.*CR.*NTU/(1+CR); %equivalent NTU

%After lookup to get priliminary effectives, net effectiveness of the
%regenerator can be estimated

epsilon = 0.79; %From lookup table
chi_re=(1-CR.^2)*(epsilon./(1-epsilon))/(2.*CR);
eff_re = (1 - exp(-chi_re))./(1-CR.*exp(-chi_re));


%Pressure drop estimation hot stream

fP=(150+3.89.*(Re./(1-phi)).^0.87).*((1-phi).^2./(phi.^3.*Re)); %Jones and Krier correlation
delP_hot(j)=fP.*mdot.^2.*L./(rho_g_hot.*dp.*Area.^2);
%Pressure drop estimation cold stream

fP_c=(150+3.89.*(Re2./(1-phi)).^0.87).*((1-phi).^2./(phi.^3*Re2));
delP_cold(j)=fP_c.*mdot^2.*L./(rho_g_cold.*dp.*Area.^2);
end

figure(1)
plot(Ph,delP_hot/1e5,'o');
grid on;
xlabel('Time of flow (seconds)');
ylabel('Hot stream pressure drop (bar)');

figure(2)
plot(Ph,delP_cold/1e5,'o');
grid on;
xlabel('Time of flow (seconds)');
ylabel('Cold stream pressure drop (bar)');

figure(3)
plot(Ph,AR,'o');
grid on;
xlabel('Time of flow (seconds)');
ylabel('Calculated Aspect Ratio for D = 2 m');

figure(4)
plot(Ph,NTUe,'o');
grid on;
xlabel('Time of flow (seconds)');
ylabel('Equivalent NTU');


