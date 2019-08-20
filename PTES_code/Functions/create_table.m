function [A] = create_table(substance)
%   CREATE_TABLE Create tables with thermophysical properties as a function of temperature.
%   Tmin and Tmax define the temperature range of the table.
%   Return the file ID where data has been written

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

% Generate thermophysical properties
if strcmp(substance,'SolarSalt') % Solar salt
    Tmin  = 10;    % Not correct. Used only for "multy-run" screening
    Tmax  = 2000;  % Not correct. Used only for "multy-run" screening

elseif strcmp(substance,'MineralOil') % Mineral oil
    Tmin  = 10;   % Not correct. Used only for "multy-run" screening
    Tbot  = 275;
    Ttop  = 600;
    Tmax  = 2000; % Not correct. Used only for "multy-run" screening
    
elseif strcmp(substance,'Methanol')
    Tmin  = 125; % Not correct. Used only for "multy-run" screening
    Tbot  = py.CoolProp.CoolProp.PropsSI('T_triple',' ',0,' ',0,'Methanol') + 0.1;
    Ttop  = py.CoolProp.CoolProp.PropsSI('T','P',1e5,'Q',0,'Methanol');
    Tmax  = 1300; % Not correct. Used only for screening purposes
    
elseif strcmp(substance,'Propane')
    Tmin  = 50; % Not correct. Used only for screening purposes
    Tbot  = py.CoolProp.CoolProp.PropsSI('T_triple',' ',0,' ',0,'Propane') + 0.1;
    Ttop  = py.CoolProp.CoolProp.PropsSI('T','P',1e5,'Q',0,'Propane'); 
    Tmax  = 1300; % Not correct. Used only for screening purposes
    
elseif strcmp(substance,'Isopentane')
    Tmin  = 50; % Not correct. Used only for screening purposes
    Tbot  = py.CoolProp.CoolProp.PropsSI('T_triple',' ',0,' ',0,'Isopentane') + 0.1;
    Ttop  = py.CoolProp.CoolProp.PropsSI('T','P',1e5,'Q',0,'Isopentane'); 
    Tmax  = 1300; % Not correct. Used only for screening purposes 
    
elseif strcmp(substance,'Ethanol')
    Tmin  = 50; % Not correct. Used only for screening purposes
    Tbot  = py.CoolProp.CoolProp.PropsSI('T_triple',' ',0,' ',0,'Ethanol') + 0.1;
    Ttop  = py.CoolProp.CoolProp.PropsSI('T','P',1e5,'Q',0,'Ethanol'); 
    Tmax  = 1300; % Not correct. Used only for screening purposes
    
elseif strcmp(substance,'Hexane')
    Tmin  = 50; % Not correct. Used only for screening purposes
    Tbot  = py.CoolProp.CoolProp.PropsSI('T_triple',' ',0,' ',0,'Hexane') + 0.1;
    Ttop  = py.CoolProp.CoolProp.PropsSI('T','P',1e5,'Q',0,'Hexane'); 
    Tmax  = 1300; % Not correct. Used only for screening purposes
    
elseif strcmp(substance,'Pentane')
    Tmin  = 50; % Not correct. Used only for screening purposes
    Tbot  = py.CoolProp.CoolProp.PropsSI('T_triple',' ',0,' ',0,'Pentane') + 0.1;
    Ttop  = py.CoolProp.CoolProp.PropsSI('T','P',1e5,'Q',0,'Pentane'); 
    Tmax  = 1300; % Not correct. Used only for screening purposes
    
elseif strcmp(substance,'EthyleneGlycol')
    Tmin  = 50; % Not correct. Used only for screening purposes
    Tbot  = 270;
    Ttop  = 550; 
    Tmax  = 1300; % Not correct. Used only for screening purposes
    
elseif strcmp(substance,'Oxygen')
    Tmin  = 20; % Not correct. Used only for screening purposes
    Tbot  = 55;
    Ttop  = 120; %10 bar 
    Tmax  = 300; % Not correct. Used only for screening purposes
    
elseif strcmp(substance,'Water')
    Tmin  = 50; % Not correct. Used only for screening purposes
    Tbot  = py.CoolProp.CoolProp.PropsSI('T_triple',' ',0,' ',0,'Water') + 0.1;
    Ttop  = py.CoolProp.CoolProp.PropsSI('T','P',1e5,'Q',0,'Water') - 0.5; 
    Tmax  = 1300; % Not correct. Used only for screening purposes  
    
elseif strcmp(substance,'Generic') %imaginary liquid
    Tmin  = 10;
    Tmax  = 1300;
    
elseif strcmp(substance,'INCOMP::MEG2[0.56]') % Ethylene Glycol solution
    Tmin  = 50; % Not correct. Used only for screening purposes
    Tbot  = py.CoolProp.CoolProp.PropsSI('Tmin',' ',0,' ',0,'INCOMP::MEG2[0.56]')+0.1;
    Ttop  = py.CoolProp.CoolProp.PropsSI('Tmax',' ',0,' ',0,'INCOMP::MEG2[0.56]')-0.1;
    Tmax  = 1300; % Not correct. Used only for screening purposes 
    
else
    fprintf(1,'substance == %s',substance);
    error('Unknown substance')
end

T  = ((Tmin):1:(Tmax))'; % temperature array, K
n  = length(T);          % get dimensions of T array
ID = zeros(n,1);         % fluid ID
Cp = zeros(n,1);         % specific heat capacity, J/kg/K
h  = zeros(n,1);         % specific enthalpy, J/kg
v  = zeros(n,1);         % specific volume, m3/kg
s  = zeros(n,1);         % specific entropy, J/kg/K

if any(strcmp(substance,{'Methanol','Propane','Isopentane','Ethanol','Hexane','Pentane','Oxygen','Water','INCOMP::MEG2[0.56]'}))
    for i=1:n
        if T(i) > Tbot
            if T(i) < Ttop
                Cp(i)  = py.CoolProp.CoolProp.PropsSI('C','T',T(i),'P',1e5,substance);
                v(i)   = 1/py.CoolProp.CoolProp.PropsSI('D','T',T(i),'P',1e5,substance);
            else
                Cp(i)  = py.CoolProp.CoolProp.PropsSI('C','T',Ttop-0.5,'P',1e5,substance);
                v(i)   = 1/py.CoolProp.CoolProp.PropsSI('D','T',Ttop-0.5,'P',1e5,substance);
            end
        else
            Cp(i)  = py.CoolProp.CoolProp.PropsSI('C','T',Tbot,'P',1e5,substance);
            v(i)   = 1/py.CoolProp.CoolProp.PropsSI('D','T',Tbot,'P',1e5,substance);
        end
    end
elseif strcmp(substance,'SolarSalt')
    %Data from Bauer et al. 2013
    for i=1:n
        Cp(i)  = 1550;   % It is essentially constant
        v(i)   = 1/1850; % Density presents only small variation from 2000 to 1750 kg/m3
    end
elseif strcmp(substance,'MineralOil')
    % Data from Shell Heat Transfer Oil S2
    T_dat   = ([  0,   20,  40,   100,  150,  200,  250,  300,  340] + 273.15);
    Cp_dat  = [1809, 1882, 1954, 2173, 2355, 2538, 2720, 2902, 3048];        
    rho_dat = [ 876,  863,  850,  811,  778,  746,  713,  681,  655];
    for i=1:n
        if T(i) > Tbot
            if T(i) < Ttop
                Cp(i) = interp1(T_dat,Cp_dat,T(i));
                v(i)  = 1/(interp1(T_dat,rho_dat,T(i)));
            else
                Cp(i) = interp1(T_dat,Cp_dat,Ttop);
                v(i)  = 1/(interp1(T_dat,rho_dat,Ttop));
            end
        else
            Cp(i) = interp1(T_dat,Cp_dat,Tbot);
            v(i)  = 1/(interp1(T_dat,rho_dat,Tbot));
        end
    end
elseif strcmp(substance,'EthyleneGlycol')
    % Data from MEGlobal 2008
    for i=1:n
        % Specific heat capacity
        a1  = 0.54467;
        a2  = 1.1854e-3;
        % Density
        T_dat   = ([50, 100, 150, 200, 250, 300, 350] + 459.67)*5/9;
        rho_dat = [1.12, 1.10, 1.08, 1.06, 1.04, 1.01, 0.99]*1000;
        % rho = b1*T + b2 = (T 1)(b1; b2),  Y   = X*B
        Y = rho_dat';
        X = [T_dat' ones(size(T_dat'))];
        B = X\Y; % B = (X\) * Y        
        if T(i) > Tbot
            if T(i) < Ttop
                Cp(i) = (a1 + a2*(T(i) - 273.15))*4187;
                v(i)  = 1/(B(1)*T(i) + B(2));
            else
                Cp(i) = (a1 + a2*(Ttop - 273.15))*4187;
                v(i)  = 1/(B(1)*Ttop + B(2));
            end
        else
            Cp(i) = (a1 + a2*(Tbot - 273.15))*4187;
            v(i)  = 1/(B(1)*Tbot + B(2));
        end
    end
elseif strcmp(substance,'Generic')
    for i=1:n
        Cp(i)  = 2000;
        v(i)   = 1/1000;
    end   
end

h(1) = 0;
s(1) = 0;
for i=1:(n-1)
    h(i+1)  = h(i) + 0.5*(Cp(i) + Cp(i+1))*(T(i+1)-T(i));    % enthalpy increase in isobaric process
    s(i+1)  = s(i) + (h(i+1) - h(i))/(0.5*(T(i+1) + T(i)));  % entropy  increase in isobaric process
end

%Save in a single matrix for printing to file
A = [ T, h, v, s, Cp];

%Print to file
fprintf(fileID,'%%%7s %20s %20s %20s %20s\n','T[K]','h[J/kg]','v[m3/kg]','s[J/(kg.K)]','Cp[J/kg/K]');
fprintf(fileID,' %7.2f %20.12e %20.12e %20.12e %20.12e\n',A');
fclose(fileID);

end

