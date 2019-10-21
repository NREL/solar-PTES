function [A] = create_table(substance)
%   CREATE_TABLE Create table A with thermophysical properties as a
%   function of temperature. Tbot and Ttop define the bottom and top
%   temperatures within which data is correct. However, a broader range
%   between Tmin and Tmax (with constant properties) is allowed for
%   screening purposes (i.e. program does not crash when reading outside
%   Tbot and Ttop).

% Set fileName according to the substance's name
if strncmp(substance,'INCOMP::',8) %skip first part of substance name
    fileName  = strcat('./Data/',substance(9:end),'.dat');
else
    fileName  = strcat('./Data/',substance,'.dat');
end

% Check if file already exists. If it does, return control to invoking
% function. Otherwise, create file and proceed.
if isfile(fileName)    
    A = load(fileName);
    return
else
    fileID = fopen(fileName,'w');
end

% Is material a liquid available in the CoolProp database?
cond1 = any(strcmp(substance,{'Methanol','Propane','Isopentane','Ethanol','Hexane','Pentane','Oxygen','Water'}));
% Is material an incompressible liquid available in the CoolProp database?
cond2 = any(strcmp(substance,{'INCOMP::MEG2[0.56]'}));
% Obtain CoolProp handle
if cond1
    ierr = 0; buffer_size = 10; herr= char((1:1:buffer_size)); backend = 'HEOS';
    handle = calllib('coolprop','AbstractState_factory',backend,substance,ierr,herr,buffer_size);
end

% Set temperature limits
Tmin = 50;    % Not correct. Used only for screening purposes
Tmax = 1300;  % Not correct. Used only for screening purposes
if cond1
    Tbot  = CP1(0,0,0,'T_triple',handle) + 0.1;
    Ttop  = CP1('PQ_INPUTS',1e5,0,'T',handle) - 0.1; % Max pressure set to 1 bar
    
elseif cond2
    Tbot = py.CoolProp.CoolProp.PropsSI('Tmin',' ',0,' ',0,substance) + 0.1;
    Ttop = py.CoolProp.CoolProp.PropsSI('Tmax',' ',0,' ',0,substance) + 0.1;
    
elseif strcmp(substance,'SolarSalt') % Solar salt
    Tbot  = 500;
    Ttop  = 880;

elseif strcmp(substance,'MineralOil') % Mineral oil
    Tbot  = 275;
    Ttop  = 600;
    
elseif strcmp(substance,'Generic') % Imaginary (constant prop.) liquid
    Tbot = Tmin;
    Ttop = Tmax;
    
else
    fprintf(1,'substance == %s',substance);
    error('Unknown substance')
end

% Set arrays
T  = ((Tmin):1:(Tmax))'; % temperature array, K
n  = length(T);          % get dimensions of T array
Cp = zeros(n,1);         % specific heat capacity, J/kg/K
h  = zeros(n,1);         % specific enthalpy, J/kg
rho= zeros(n,1);         % density, kg/m3
s  = zeros(n,1);         % specific entropy, J/kg/K
k  = zeros(n,1);         % conductivity, W/(m.K)
mu = zeros(n,1);         % viscosity, Pa.s

% Obtain thermophysical properties
bot = T<Tbot;
top = T>Ttop;
mid = ~bot & ~top;
if cond1 % Extract data from CoolProp
    
    [Cp(mid),rho(mid),k(mid),mu(mid),~]...
        = CP5('PT_INPUTS',1e5*ones(size(T(mid))),T(mid),'C','D','CONDUCTIVITY','VISCOSITY','T',handle);
    
elseif cond2 % Extract data from CoolProp (requires access to  Python interface)
    
    for i=find(mid,1,'first'):find(mid,1,'last')
        Cp(i)  = py.CoolProp.CoolProp.PropsSI('C','T',T(i),'P',1e5,substance);
        rho(i) = py.CoolProp.CoolProp.PropsSI('D','T',T(i),'P',1e5,substance);
        k(i)   = py.CoolProp.CoolProp.PropsSI('CONDUCTIVITY','T',T(i),'P',1e5,substance);
        mu(i)  = py.CoolProp.CoolProp.PropsSI('VISCOSITY','T',T(i),'P',1e5,substance);
    end
    
elseif strcmp(substance,'SolarSalt')    
    % Data from SQM Thermo-solar Salts
    T_C = T - 273.15;
    Cp(mid)  = 1443 + 0.172.*T_C(mid);
    rho(mid) = 2090 - 0.636.*T_C(mid);
    k(mid)   = 0.443 + 1.9e-4.*T_C(mid);
    mu(mid)  = 1e-3*( 22.714 - 0.120.*T_C(mid) + 2.281e-4.*T_C(mid).^2 - 1.474e-7.*T_C(mid).^3);
    
elseif strcmp(substance,'MineralOil')
    % Data from Shell Heat Transfer Oil S2
    T_dat   = [    0,   20,   40,  100,  150,  200,  250,  300,  340] + 273.15;
    Cp_dat  = [ 1809, 1882, 1954, 2173, 2355, 2538, 2720, 2902, 3048];        
    rho_dat = [  876,  863,  850,  811,  778,  746,  713,  681,  655];
    k_dat   = [0.136 0.134 0.133 0.128 0.125 0.121 0.118 0.114 0.111];
    Pr_dat  = [ 3375,  919,  375,   69,   32,   20,   14,   11,    9];
    mu_dat  = Pr_dat.*k_dat./Cp_dat;
    
    Cp(mid)  = interp1(T_dat,Cp_dat ,T(mid));
    rho(mid) = interp1(T_dat,rho_dat,T(mid));
    k(mid)   = interp1(T_dat,k_dat  ,T(mid));
    mu(mid)  = interp1(T_dat,mu_dat ,T(mid));
    
elseif strcmp(substance,'Generic')
    Cp  = 1500.*ones(size(T));
    rho = 2000.*ones(size(T));
    k   =  0.5.*ones(size(T));
    mu  = 2e-3.*ones(size(T));
end

% Set constant properties outside of normal range
Cp(bot)  = Cp( find(mid,1,'first')).*ones(size(T(bot)));
rho(bot) = rho(find(mid,1,'first')).*ones(size(T(bot)));
k(bot)   = k(  find(mid,1,'first')).*ones(size(T(bot)));
mu(bot)  = mu( find(mid,1,'first')).*ones(size(T(bot)));

Cp(top)  = Cp( find(mid,1,'last')).*ones(size(T(top)));
rho(top) = rho(find(mid,1,'last')).*ones(size(T(top)));
k(top)   = k(  find(mid,1,'last')).*ones(size(T(top)));
mu(top)  = mu( find(mid,1,'last')).*ones(size(T(top)));

% Compute derived parameters
v  = 1./rho;    % specific volume, m3/kg
Pr = Cp.*mu./k; % Prandtl number
h(1) = 0;
s(1) = 0;
for i=1:(n-1)
    h(i+1)  = h(i) + 0.5*(Cp(i) + Cp(i+1))*(T(i+1)-T(i));    % enthalpy increase in isobaric process
    s(i+1)  = s(i) + (h(i+1) - h(i))/(0.5*(T(i+1) + T(i)));  % entropy  increase in isobaric process
end

%Save in a single matrix for printing to file
A = [ T, h, v, s, Cp, k, mu, Pr];

%Print to file
fprintf(fileID,'%%%7s %20s %20s %20s %20s %20s %20s %20s\n',...
    'T[K]','h[J/kg]','v[m3/kg]','s[J/(kg.K)]','Cp[J/kg/K]','k[W/m.K]','mu[Pa.s]','Pr[-]');
fprintf(fileID,' %7.2f %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n',A');
fclose(fileID);

end

