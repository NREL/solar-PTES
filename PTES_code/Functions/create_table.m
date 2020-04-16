function [A] = create_table(Mname)
%CREATE_TABLE Create a table of thermophysical properties.
%
%   The thermophysical properties of a material as a function of
%   temperature are stored in the matrix 'A'. 'Mname' is the name of the
%   material. 'Tbot' and 'Ttop' define the bottom and top temperatures
%   within which data is correct. However, a broader range between 'Tmin'
%   and 'Tmax' (with constant properties) is allowed for screening purposes
%   (i.e. program does not crash when reading outside 'Tbot' and 'Ttop').
%   
%   This function can be used either for incompressible fluids or solid
%   materials. The matrix 'A' can be used as a property of a 'fluid_class'
%   or a 'packbed_class' object. Properties can be read manually from 'A'
%   or using the RP1 function.
%
%   For incompressible mixtures, creation of CoolProp handles to use with
%   the low-level interface is not supported. Therefore, properties are
%   obtained via the Python interface. Once the table is created, the
%   Python interface is not used anywhere else on the PTES code.
%
%   USAGE e.g.:
%   TAB     = create_table('Water');
%   obj.TAB = create_table(Fname);
%   pbH(ii).sld   = create_table(pbH(ii).Sname);

% Set global TABS structure (which contains tabular information on various
% fluids)
global TABS

% Set fileName
fileName  = strcat('./Data/',valid_name(Mname),'.dat');

% Check if file already exists. If it does, load the table, save it into
% the matrix A and into the global TABS structure, and return control to
% invoking function. Otherwise, create file and proceed.
if isfile(fileName)    
    A = load(fileName);
    TABS.(valid_name(Mname)) = A;
    return
else
    fileID = fopen(fileName,'w');
end

% Is material a liquid available in the CoolProp database?
cond1 = any(strcmp(Mname,{'Methanol','Propane','Isopentane','Ethanol',...
    'Hexane','Pentane','Oxygen','Water'}));

% Is material an incompressible liquid available in the CoolProp database?
cond2 = any(strcmp(Mname,{'INCOMP::MEG2[0.56]'}));

% Obtain CoolProp handle
if cond1
    ierr = 0; buffer_size = 10; herr= char((1:1:buffer_size)); backend = 'HEOS';
    handle = calllib('coolprop','AbstractState_factory',backend,Mname,ierr,herr,buffer_size);
end

% Set temperature limits
Tmin = 50;   % Artificial limit. Used only for screening purposes
Tmax = 1300; % Artificial limit. Used only for screening purposes
if cond1
    Tbot  = CP1(0,0,0,'T_triple',handle) + 0.1;
    Ttop  = CP1('PQ_INPUTS',1e5,0,'T',handle) - 0.1; % Max pressure set to 1 bar
    
elseif cond2
    Tbot = py.CoolProp.CoolProp.PropsSI('Tmin',' ',0,' ',0,Mname) + 0.1;
    Ttop = py.CoolProp.CoolProp.PropsSI('Tmax',' ',0,' ',0,Mname) + 0.1;
    
else
    switch Mname
        case 'SolarSalt'  % Solar salt
            Tbot  = 500;
            Ttop  = 880;
            
        case 'MineralOil' % Mineral oil
            Tbot  = 275;
            Ttop  = 600;
            
        case 'Magnetite'  % Magnetite
            Tbot  = 100;
            Ttop  = 1200;
            
        case 'SiO2'       % Silicon Oxide
            Tbot  = 100;
            Ttop  = 1200;
            
        case 'Generic'    % Imaginary material with constant properties
            Tbot = Tmin;
            Ttop = Tmax;
            
        otherwise
            fprintf(1,'substance == %s',Mname);
            error('Unknown substance')
    end
end

% Pre-allocate arrays
T  = ((Tmin):1:(Tmax))'; % temperature array, K
n  = length(T);          % get dimensions of T array
Cp = zeros(n,1);         % specific heat capacity, J/kg/K
h  = zeros(n,1);         % specific enthalpy, J/kg
rho= zeros(n,1);         % density, kg/m3
s  = zeros(n,1);         % specific entropy, J/kg/K
k  = zeros(n,1);         % conductivity, W/(m.K)
mu = zeros(n,1);         % viscosity, Pa.s
Q  = -1.0*ones(n,1);     % vapour quality

% Create logical arrays defining the different 'temperature regions'
bot = T<Tbot;
top = T>Ttop;
mid = ~bot & ~top;

% Obtain thermophysical properties
if cond1 % Extract data from CoolProp. Assume fluid is kept at ambient pressure.
    p = 1e5*ones(size(T));
    [Cp(mid)]   = CP1('PT_INPUTS',p(mid),T(mid),'C',handle);
    [rho(mid)]  = CP1('PT_INPUTS',p(mid),T(mid),'D',handle);
    [k(mid)]    = CP1('PT_INPUTS',p(mid),T(mid),'CONDUCTIVITY',handle);
    [mu(mid)]   = CP1('PT_INPUTS',p(mid),T(mid),'VISCOSITY',handle);
    
elseif cond2 % Extract data from CoolProp (requires access to  Python interface)
    
    for i=find(mid,1,'first'):find(mid,1,'last')
        Cp(i)  = py.CoolProp.CoolProp.PropsSI('C','T',T(i),'P',1e5,Mname);
        rho(i) = py.CoolProp.CoolProp.PropsSI('D','T',T(i),'P',1e5,Mname);
        k(i)   = py.CoolProp.CoolProp.PropsSI('CONDUCTIVITY','T',T(i),'P',1e5,Mname);
        mu(i)  = py.CoolProp.CoolProp.PropsSI('VISCOSITY','T',T(i),'P',1e5,Mname);
    end
    
else
    switch Mname
        case 'SolarSalt'
            % Data from SQM Thermo-solar Salts
            T_C = T - 273.15;
            Cp(mid)  = 1443 + 0.172.*T_C(mid);
            rho(mid) = 2090 - 0.636.*T_C(mid);
            k(mid)   = 0.443 + 1.9e-4.*T_C(mid);
            mu(mid)  = 1e-3*( 22.714 - 0.120.*T_C(mid) + 2.281e-4.*T_C(mid).^2 - 1.474e-7.*T_C(mid).^3);
            
        case 'MineralOil'
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
            
        case 'Magnetite'
            % Data from JANEF Tables and McTigue thesis
            T_dat   = [ 10,15,20,25,30,35,40,45,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,350,400,500,600,700,800,900,1000,1100,1200,1300,1400];
            Cp_dat  = [ 0.487877348,1.825022673,4.46317426,9.341947743,16.76852516,26.38151587,37.62076441,50.39592313,64.29139279,94.52171885,127.2998488,161.9391406,198.4757331,236.0422889,280.384919,312.0337033,343.6824876,370.6060894,395.5420428,419.3938242,442.1614338,463.6641762,484.263442,503.5978406,522.209458,539.5562082,556.1801771,572.0813647,587.2597711,601.8960916,615.8096307,629.0003887,641.4683654,653.3942561,706.6993738,739.106,830.904,918.007,1005.109,1092.212,1179.324,1179.324,1179.324,1179.324,1179.324,1179.324];
            rho_dat = 5175 .* ones(size(T_dat)) ;
            k_dat   = 5.0 .* ones(size(T_dat)) ;
            
            Cp(mid)  = interp1(T_dat,Cp_dat ,T(mid));
            rho(mid) = interp1(T_dat,rho_dat,T(mid));
            k(mid)   = interp1(T_dat,k_dat  ,T(mid));
            mu(mid)  = 0.0 ;
            
        case 'SiO2'
            % Data from JANEF Tables and McTigue thesis
            T_dat   = [0, 100, 200, 298.15, 300, 400, 500, 600, 700, 800, 847];
            Cp_dat  = [0, 261.072, 543.232, 742.123, 745.119, 889.270, 992.677, 1072.134, 1144.55, 1226.653, 1273.388];
            rho_dat = 2660 .* ones(size(T_dat)) ;
            k_dat   = 5.0 .* ones(size(T_dat)) ;
            
            Cp(mid)  = interp1(T_dat,Cp_dat ,T(mid));
            rho(mid) = interp1(T_dat,rho_dat,T(mid));
            k(mid)   = interp1(T_dat,k_dat  ,T(mid));
            mu(mid)  = 0.0 ;
            
        case 'Generic'
            Cp  = 1500.*ones(size(T));
            rho = 2000.*ones(size(T));
            k   =  0.5.*ones(size(T));
            mu  = 2e-3.*ones(size(T));
        otherwise
            error('not implemented')
    end
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
Pr = Cp.*mu./k; % Prandtl number
h(1) = 0;
s(1) = 0;
for i=1:(n-1)
    h(i+1)  = h(i) + 0.5*(Cp(i) + Cp(i+1))*(T(i+1)-T(i));    % enthalpy increase in isobaric process
    s(i+1)  = s(i) + (h(i+1) - h(i))/(0.5*(T(i+1) + T(i)));  % entropy  increase in isobaric process
end

% Save in a single matrix for printing to file
A = [ T, h, rho, s, Cp, k, mu, Pr, Q];

% Print to file
fprintf(fileID,'%%%7s %20s %20s %20s %20s %20s %20s %20s %11s\n',...
    'T[K]','h[J/kg]','rho[m3/kg]','s[J/(kg.K)]','Cp[J/kg/K]','k[W/m.K]','mu[Pa.s]','Pr[-]','Q[-]');
fprintf(fileID,' %7.2f %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %11.6f\n',A');
fclose(fileID);

% Save into global TABS structure
TABS.(valid_name(Mname)) = A;

end

