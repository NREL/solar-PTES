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

%gas = fluid_class('CarbonDioxide','WF','CP','TTSE',Load.num,30);
gas = fluid_class('Nitrogen','WF','CP','BICUBIC&HEOS',Load.num,30);

CMP = compexp_class('comp', 'isen', 2, eta0, Load.num) ; % Compressors
EXP = compexp_class('exp', 'isen', 2, eta0, Load.num) ; % Expanders

gas.state(1,1).p    = p0; 
gas.state(1,1).T    = T0;
gas.state(1,1).mdot = Load.mdot(1);
[gas]               = update(gas,[1,1],1);

% *** COMPRESSOR *** %
% Design point
[CMP,~,~] = compexp_offdesign (CMP , gas.state(1,1), 1, 1 , 1) ;
Nlo = 0.6 * 2 * pi * 3600;
Nhi = 1.* 2 * pi * 3600;
nj  = 5 ;
CMP.mdot0 = Load.mdot(1);
CMP.Tin = T0;
CMP.Pin = p0;

% The max and min reduced mass flows depend on the reduced speed
mlo(1) = 0.4;
mhi(1) = 0.6;

mlo(2) = 0.4;
mhi(2) = 0.7;

mlo(3) = 0.5;
mhi(3) = 0.8;

mlo(4) = 0.6;
mhi(4) = 0.9;

mlo(5) = 0.7;
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
        [CMP,Nr(i,j),mr(i,j)] = compexp_offdesign (CMP , gas.state(1,1), 2 , 1, 1) ;
        
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
legend boxoff

figure(2)
plot(mr,pr)

xlabel('Reduced mass flow')
ylabel('Reduced pressure ratio $$\beta / \beta_0$$')
xlim([0.4 1.1])
ylim([0.0 1.1])
legend({'$$\dot{n}_c = $$ 0.6','$$\dot{n}_c = $$ 0.7','$$\dot{n}_c = $$ 0.8','$$\dot{n}_c = $$ 0.9','$$\dot{n}_c = $$ 1.0'},'Location','Best')
legend boxoff

figure(3)
N = 5;
for i = 2:N
    p(i)=plot(mr(:,i),pr(:,i));hold on;
    p(i).LineWidth = 2.5;
    %p(i).Color = (0.8 * [1 1 1] * ((N-i) / N)) ;
end
p(2).Color = 0.2 * [1 1 1];
p(3).Color = 0.4 * [1 1 1];
p(4).Color = 0.6 * [1 1 1];
p(5).Color = 0.8 * [1 1 1];
xlim([0.4 1.1])
ylim([0.5 1.2])
m1 = plot(0.99829,1.002,'o');
m1.MarkerEdgeColor = 0.2 * [1 1 1];
m1.MarkerSize = 8;

m3 = plot(0.79592,0.81959,'s');
m3.MarkerEdgeColor = 0.4 * [1 1 1];
m3.MarkerSize = 8;

a2 = annotation('arrow',[0.72 0.563],[0.69 0.51]);
a2.Color = 0.2 * [1 1 1];
a2.LineWidth = 1;
a2.HeadStyle = "vback3";
a2.HeadLength = 8 ;
pbaspect([1 1 1])
hold off;

xlabel('Reduced mass flow $$m_r$$')
%xlabel('Reduced mass flow $$\dot{m} / \dot{m}_o$$')
ylabel('Relative pressure ratio $$\beta / \beta_0$$')
legend({'$$n_r = $$ 0.7','$$n_r = $$ 0.8','$$n_r = $$ 0.9','$$n_r = $$ 1.0'},'Location','northwest')
legend boxoff



figure(4)
N = 5;
for i = 2:N
    p(i)=plot(mr(:,i),eta(:,i)/eta0);hold on;
    p(i).LineWidth = 2.5;
end
p(2).Color = 0.2 * [1 1 1];
p(3).Color = 0.4 * [1 1 1];
p(4).Color = 0.6 * [1 1 1];
p(5).Color = 0.8 * [1 1 1];
xlim([0.4 1.1])
ylim([0.8 1.02])
m1 = plot(1.002,1.,'o');
m1.MarkerEdgeColor = 0.2 * [1 1 1];
m1.MarkerSize = 8;

m3 = plot(0.80204,0.98213,'s');
m3.MarkerEdgeColor = 0.4 * [1 1 1];
m3.MarkerSize = 8;

%a1 = annotation('arrow',([0.99829 0.79796]-0.4)./0.7,([1.002 1.1117]-0.5)./0.7);
a2 = annotation('arrow',[0.715 0.57],[0.85 0.79]);
a2.Color = 0.2 * [1 1 1];
a2.LineWidth = 1;
a2.HeadStyle = "vback3";
a2.HeadLength = 8 ;
pbaspect([1 1 1])
hold off;

xlabel('Reduced mass flow $$m_r$$')
ylabel('Relative efficiency $$\eta / \eta_0$$')
legend({'$$n_r = $$ 0.7','$$n_r = $$ 0.8','$$n_r = $$ 0.9','$$n_r = $$ 1.0'},'Location','southeast')
legend boxoff



% *** EXPANDER *** %
% Design point
gas.state(1,1).p    = 10.*p0; 
gas.state(1,1).T    = 500;
[EXP,~,~] = compexp_offdesign (EXP , gas.state(1,1), 1 , 1, 1) ;
EXP.pr0 = 1.5 ;
Nlo = 0.6 * 2 * pi * 3600;
Nhi = 1.* 2 * pi * 3600;
nj  = 5 ;

EXP.mdot0 = Load.mdot(1);
EXP.Tin = gas.state(1,1).T;
EXP.Pin = gas.state(1,1).p;


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
        [EXP,Nr(i,j),mr(i,j)] = compexp_offdesign (EXP , gas.state(1,1), 2 , 1, 1) ;
        
        eta(i,j) = EXP.eta(2) ;
        pr(i,j)  = EXP.pr(2) ;
        
    end

end
    
% Plot expander results
figure(5)
plot(mr,eta/eta0)

xlabel('Reduced mass flow')
ylabel('Reduced isentropic efficiency $$\eta / \eta_0$$')
xlim([0.4 1.1])
ylim([0.6 1.0])
legend({'$$\dot{n}_t = $$ 0.6','$$\dot{n}_t = $$ 0.7','$$\dot{n}_t = $$ 0.8','$$\dot{n}_t = $$ 0.9','$$\dot{n}_t = $$ 1.0'},'Location','Best')
legend boxoff

figure(6)
pr(pr<=0) = nan;
plot(mr,pr)

xlabel('Reduced mass flow')
ylabel('Reduced pressure ratio $$\beta / \beta_0$$')
xlim([0.4 1.2])
ylim([0. 1.2])
legend({'$$\dot{n}_t = $$ 0.6','$$\dot{n}_t = $$ 0.7','$$\dot{n}_t = $$ 0.8','$$\dot{n}_t = $$ 0.9','$$\dot{n}_t = $$ 1.0'},'Location','Best')
legend boxoff
    
figure(7)
N = 5;
for i = 2:N
    p(i)=plot(mr(:,i),pr(:,i));hold on;
    p(i).LineWidth = 2.5;
end
p(2).Color = 0.2 * [1 1 1];
p(3).Color = 0.4 * [1 1 1];
p(4).Color = 0.6 * [1 1 1];
p(5).Color = 0.8 * [1 1 1];
xlim([0.4 1.1])
ylim([0.5 1.2])
m1 = plot(0.99829,0.999,'o');
m1.MarkerEdgeColor = 0.1 * [1 1 1];
m1.MarkerSize = 8;


m3 = plot(0.79592,0.81959,'s');
m3.MarkerEdgeColor = 0.1 * [1 1 1];
m3.MarkerSize = 8;

%a1 = annotation('arrow',([0.99829 0.79796]-0.4)./0.7,([1.002 1.1117]-0.5)./0.7);
a2 = annotation('arrow',[0.72 0.563],[0.69 0.51]);
a2.Color = 0.2 * [1 1 1];
a2.LineWidth = 1;
a2.HeadStyle = "vback3";
a2.HeadLength = 8 ;
pbaspect([1 1 1])
hold off;

xlabel('Reduced mass flow $$m_r$$')
ylabel('Relative pressure ratio $$\beta / \beta_0$$')
legend({'$$n_r = $$ 0.7','$$n_r = $$ 0.8','$$n_r = $$ 0.9','$$n_r = $$ 1.0'},'Location','northwest')
legend boxoff



figure(8)
N = 5;

for i = 2:N
    p(i)=plot(mr(:,i),eta(:,i)/eta0);hold on;
    p(i).LineWidth = 2.5;
end
p(2).Color = 0.2 * [1 1 1];
p(3).Color = 0.4 * [1 1 1];
p(4).Color = 0.6 * [1 1 1];
p(5).Color = 0.8 * [1 1 1];
xlim([0.4 1.1])
ylim([0.8 1.02])
m1 = plot(1.0,1.0,'o');
m1.MarkerEdgeColor = 0.1 * [1 1 1];
m1.MarkerSize = 8;

m3 = plot(0.79799,0.98071,'s');
m3.MarkerEdgeColor = 0.1 * [1 1 1];
m3.MarkerSize = 8;

a2 = annotation('arrow',[0.715 0.57],[0.85 0.79]);
a2.Color = 0.2 * [1 1 1];
a2.LineWidth = 1;
a2.HeadStyle = "vback3";
a2.HeadLength = 8 ;
pbaspect([1 1 1])
hold off;

xlabel('Reduced mass flow $$m_r$$')
ylabel('Relative efficiency $$\eta / \eta_0$$')
legend({'$$n_r = $$ 0.7','$$n_r = $$ 0.8','$$n_r = $$ 0.9','$$n_r = $$ 1.0'},'Location',[0.5, 0.15, 0.25, 0.25])
legend boxoff

