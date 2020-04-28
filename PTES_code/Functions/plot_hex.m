function [] = plot_hex(HX,iL,fignum,C_or_K)
% PLOT_HEX_TQA Make T-Q, T-A and p-A diagrams of a heat exchanger.
% Use data stored in the HX structure.
% If the AS array (cummulative area) does not exist (i.e. the HX structure
% was created by the hex_TQ function, rather than the hex_TQA function),
% only the T-Q diagram is made. 

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
plot(QS./QS(end),H.T-K,'r'); hold on;
plot(QS./QS(end),C.T-K,'b'); hold off;
xlabel('Normalised cumulative heat transfer')
ylabel(ytext)
legend(L,'Location','Best')

if ~isempty(AS)
    % Plot T-A diagram
    figure(fignum+1)
    plot(AS./AS(end),H.T-K,'r'); hold on;
    plot(AS./AS(end),C.T-K,'b'); hold off;
    xlabel('Normalised cumulative heat transfer area')
    ylabel(ytext)
    legend(L,'Location','Best')
    
    % Plot p-A diagram
    figure(fignum+2)
    plot(AS./AS(end),H.p/H.pin*100,'r'); hold on;
    plot(AS./AS(end),C.p/C.pin*100,'b'); hold off;
    xlabel('Normalised cumulative heat transfer area')
    ylabel('Normalised pressure [$\%$]')
    legend(L,'Location','Best')
    
    % Plot Re-A diagram
    figure(fignum+3)
    semilogy(AS./AS(end),H.Re,'r'); hold on;
    semilogy(AS./AS(end),C.Re,'b'); hold off;
    xlabel('Normalised cumulative heat transfer area')
    ylabel('Reynolds number')
    legend(L,'Location','Best')

    % Plot ht-Q diagram
    figure(fignum+4)
    semilogy(QS./QS(end),H.ht, 'r'); hold on;
    semilogy(QS./QS(end),C.ht, 'b'); hold on;
    semilogy(QS./QS(end),H.ht/2, 'r--'); hold on;
    semilogy(QS./QS(end),HX.Ul,'k-.'); hold off;
    xlabel('Normalised cumulative heat transfer')
    ylabel('Heat transfer coefficient [$\mathrm{W/m^2/K}$]')
    legend([L(:)',{'Overall'}],'Location','Best')
    legend([L(:)',{'Minimum','Overall'}],'Location','Best')
    ylim([40 3000])
end

end