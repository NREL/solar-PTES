% Plot out characteristic curves of compressors and turbines
%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clear;

% Enter debugging mode if an error occurs
%dbstop if error
%dbclear all

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop

% Set atmospheric conditions and cycle parameters
T0      = 30 + 273.15;  % ambient temp, K
p0      = 1e5;          % ambient pressure, Pa
% Set component parameters
eta0   = 0.90;  % polytropic efficiency

Load.time = [10; 10].*3600;               % time spent in each load period, s
Load.type = ["chg"; "chg"];    % type of load period
Load.mdot = [1; 1];                       % working fluid mass flow rate, kg/s
Load.num  = numel(Load.time);
Load.ind  = 1:Load.num;

gas = fluid_class('Nitrogen','WF','CP','BICUBIC&HEOS',Load.num,30);

CMP = compexp_class('comp', 'isen', 2, eta0, Load.num) ; % Compressors
EXP = compexp_class('exp', 'isen', 2, eta0, Load.num) ; % Expanders

gas.state(1,1).p    = p0; 
gas.state(1,1).T    = T0;
gas.state(1,1).mdot = Load.mdot(1);
[gas]               = update(gas,[1,1],1);

% *** COMPRESSOR *** %
% Design point
[CMP,~,~] = compexp_offdesign (CMP , gas.state(1,1), 1 , 1) ;
Nlo = 0.6 * 2 * pi * 3600;
Nhi = 1.* 2 * pi * 3600;
nj  = 5 ;

% The max and min reduced mass flows depend on the reduced speed
mlo(1) = 0.4;
mhi(1) = 0.6;

mlo(2) = 0.5;
mhi(2) = 0.7;

mlo(3) = 0.6;
mhi(3) = 0.8;

mlo(4) = 0.7;
mhi(4) = 0.9;

mlo(5) = 0.8;
mhi(5) = 1.1;

ni  = 50 ;

Nr  = zeros(ni,nj) ;
mr  = zeros(ni,nj) ;
eta = zeros(ni,nj) ;
pr  = zeros(ni,nj) ;

for j = 1:nj

    n = Nlo + (j-1.)*(Nhi-Nlo)/(nj-1.) ;
    CMP.N(2) = n;
        
    for i = 1:ni
        m = mlo(j) + (i-1.)*(mhi(j)-mlo(j))/(ni-1.) ;
        gas.state(1,1).mdot = m;
        [CMP,Nr(i,j),mr(i,j)] = compexp_offdesign (CMP , gas.state(1,1), 2 , 1) ;
        
        eta(i,j) = CMP.eta(2) ;
        pr(i,j)  = CMP.pr(2) ;
        
    end

end
    
% Plot compressor results
figure(1)
plot(mr,eta/eta0)

xlabel('Reduced mass flow')
ylabel('Reduced isentropic efficiency $$\eta / \eta_0$$')
xlim([0.4 1.1])
ylim([0.8 1.0])
legend({'$$\dot{n}_c = $$ 0.6','$$\dot{n}_c = $$ 0.7','$$\dot{n}_c = $$ 0.8','$$\dot{n}_c = $$ 0.9','$$\dot{n}_c = $$ 1.0'},'Location','Best')


figure(2)
plot(mr,pr)

xlabel('Reduced mass flow')
ylabel('Reduced pressure ratio $$\beta / \beta_0$$')
xlim([0.4 1.1])
ylim([0.0 1.1])
legend({'$$\dot{n}_c = $$ 0.6','$$\dot{n}_c = $$ 0.7','$$\dot{n}_c = $$ 0.8','$$\dot{n}_c = $$ 0.9','$$\dot{n}_c = $$ 1.0'},'Location','Best')
    



% *** EXPANDER *** %
% Design point
[EXP,~,~] = compexp_offdesign (EXP , gas.state(1,1), 1 , 1) ;
EXP.pr0 = 10. ;
Nlo = 0.6 * 2 * pi * 3600;
Nhi = 1.* 2 * pi * 3600;
nj  = 5 ;

% The max and min reduced mass flows depend on the reduced speed
mlo(1) = 0.2;
mhi(1) = 1.2;

mlo(2) = 0.2;
mhi(2) = 1.2;

mlo(3) = 0.2;
mhi(3) = 1.2;

mlo(4) = 0.2;
mhi(4) = 1.2;

mlo(5) = 0.2;
mhi(5) = 1.2;

ni  = 200 ;

for j = 1:nj

    n = Nlo + (j-1.)*(Nhi-Nlo)/(nj-1.) ;
    EXP.N(2) = n;
        
    for i = 1:ni
        m = mlo(j) + (i-1.)*(mhi(j)-mlo(j))/(ni-1.) ;
        gas.state(1,1).mdot = m;
        [EXP,Nr(i,j),mr(i,j)] = compexp_offdesign (EXP , gas.state(1,1), 2 , 1) ;
        
        eta(i,j) = EXP.eta(2) ;
        pr(i,j)  = EXP.pr(2) ;
        
    end

end
    
% Plot compressor results
figure(3)
plot(mr,eta/eta0)

xlabel('Reduced mass flow')
ylabel('Reduced isentropic efficiency $$\eta / \eta_0$$')
xlim([0.4 1.1])
ylim([0.6 1.0])
legend({'$$\dot{n}_t = $$ 0.6','$$\dot{n}_t = $$ 0.7','$$\dot{n}_t = $$ 0.8','$$\dot{n}_t = $$ 0.9','$$\dot{n}_t = $$ 1.0'},'Location','Best')


figure(4)
pr(pr<=0) = nan;
plot(mr,pr/EXP.pr0)

xlabel('Reduced mass flow')
ylabel('Reduced pressure ratio $$\beta / \beta_0$$')
xlim([0.4 1.2])
ylim([0. 1.2])
legend({'$$\dot{n}_t = $$ 0.6','$$\dot{n}_t = $$ 0.7','$$\dot{n}_t = $$ 0.8','$$\dot{n}_t = $$ 0.9','$$\dot{n}_t = $$ 1.0'},'Location','Best')
    
    