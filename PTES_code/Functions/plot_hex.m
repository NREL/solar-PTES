function [] = plot_hex(HX,iL,fignum,C_or_K,Lsimple)
% PLOT_HEX_TQA Make T-Q, T-A and p-A diagrams of a heat exchanger.
% Use data stored in the HX structure.
% If the AS array (cummulative area) does not exist (i.e. the HX structure
% was created by the hex_TQ function, rather than the hex_TQA function),
% only the T-Q diagram is made.
% Use the Lsimple logical to force a single plot output.

% Extract inputs
C  = HX.C(iL);
H  = HX.H(iL);
QS = HX.QS(iL,:);
AS = HX.AS;

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

% Set legend
L = {[valid_name(H.name,2),', ',sprintf('%.1f',H.pin/1e5),' bar'],[valid_name(C.name,2),', ',sprintf('%.1f',C.pin/1e5),' bar']};

% Plot TQ diagram
figure(fignum)
clf(fignum)
plot(QS./QS(end),H.T-K,'r'); hold on;
plot(QS./QS(end),C.T-K,'b'); hold off;
xlabel('Normalised cumulative heat transfer')
ylabel(ytext)
legend(L,'Location','SouthEast')
dim = [.2 .6 .3 .3];
str = {sprintf('DppH = %.3f',HX.DppH(iL)),sprintf('DppC = %.3f',HX.DppC(iL))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

if ~isempty(AS) && ~Lsimple
    % Plot T-A diagram
    figure(fignum+1)
    plot(AS./AS(end),H.T-K,'r'); hold on;
    plot(AS./AS(end),C.T-K,'b'); hold off;
    xlabel('Normalised cumulative heat transfer area')
    ylabel(ytext)
    legend(L,'Location','Best')
    
    % Plot p-Q diagram
    figure(fignum+2)
    plot(QS./QS(end),H.p/H.pin*100,'r'); hold on;
    plot(QS./QS(end),C.p/C.pin*100,'b'); hold off;
    xlabel('Normalised cumulative heat transfer')
    ylabel('Normalised pressure [$\%$]')
    legend(L,'Location','Best')
    
    % Plot Re-Q diagram
    figure(fignum+3)
    semilogy(QS./QS(end),H.Re,'r'); hold on;
    semilogy(QS./QS(end),C.Re,'b'); hold off;
    xlabel('Normalised cumulative heat transfer')
    ylabel('Reynolds number')
    legend(L,'Location','Best')

    % Plot ht-Q diagram
    figure(fignum+4)
    semilogy(QS./QS(end),H.ht, 'r'); hold on;
    semilogy(QS./QS(end),C.ht, 'b'); hold on;
    %semilogy(QS./QS(end),H.ht/2, 'r--'); hold on;
    semilogy(QS./QS(end),HX.Ul,'k-.'); hold off;
    xlabel('Normalised cumulative heat transfer')
    ylabel('Heat transfer coefficient [$\mathrm{W/m^2/K}$]')
    legend([L(:)',{'Overall'}],'Location','Best')
    %legend([L(:)',{'Minimum','Overall'}],'Location','Best')
    %ylim([40 3000])
    
    % Plot dpdL-Q diagram
    figure(fignum+5)
    semilogy(QS./QS(end),-H.dpdL, 'r'); hold on;
    semilogy(QS./QS(end),-C.dpdL, 'b'); hold off;
    xlabel('Normalised cumulative heat transfer')
    ylabel('Pressure gradient [Pa/m]')
    legend(L,'Location','Best')
end

end