function [] = plot_hex_TQA(HX,C_or_K)
% PLOT_HEX_TQA Make T-Q, T-A and p-A diagrams of a heat exchanger.
% Use data stored in the HX structure.

% Extract inputs
C  = HX.C;
H  = HX.H;
AS = HX.AS;
QS = HX.QS;

% Set plots in either Celcius or Kelvin
if strcmp(C_or_K,'C')
    K = 273.15;
    ytext = 'Temperature [$^\circ $C]';
elseif strcmp(C_or_K,'K')
    K = 0;
    ytext = 'Temperature [K]';
else
    error('C_or_K must be set to either "C" or "K"');
end

% Plot TQ diagram
figure(20)
plot(QS./QS(end),H.T-K,'r'); hold on;
plot(QS./QS(end),C.T-K,'b'); hold off;
xlabel('Cumulative heat transfer')
ylabel(ytext)
legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],[C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')

% Plot TA diagram
figure(21)
plot(AS./AS(end),H.T-K,'r'); hold on;
plot(AS./AS(end),C.T-K,'b'); hold off;
xlabel('Cumulative heat transfer area')
ylabel(ytext)
legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],[C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')

% Plot pA diagram
figure(22)
plot(AS./AS(end),H.p/H.pin*100,'r'); hold on;
plot(AS./AS(end),C.p/C.pin*100,'b'); hold off;
xlabel('Cumulative heat transfer area')
ylabel('Normalised pressure [$\%$]')
legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],[C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')

% Plot pA diagram
figure(23)
plot(AS./AS(end),H.Re,'r'); hold on;
plot(AS./AS(end),C.Re,'b'); hold off;
xlabel('Cumulative heat transfer area')
ylabel('Reynolds number')
legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],[C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')

% Plot pA diagram
figure(24)
plot(AS./AS(end),H.G,'r'); hold on;
plot(AS./AS(end),C.G,'b'); hold off;
xlabel('Cumulative heat transfer area')
ylabel('Mass flux, kg/s/m2')
legend([H.name,', ',sprintf('%.1f',H.pin/1e5),' bar'],[C.name,', ',sprintf('%.1f',C.pin/1e5),' bar'],'Location','Best')
end