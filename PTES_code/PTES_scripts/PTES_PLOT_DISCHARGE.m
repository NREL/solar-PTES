% Plot discharge cycle (PTES)

plot_file = load ('./Outputs/Plot_file.txt');
n1 = stages_ch*num;
n2 = (stages_ch + stages_dis)*num;
figure(1);

% Plot curves
plot(plot_file(1:n1,4)/1000,plot_file(1:n1,1),'k-','LineWidth',2); hold on;
plot(plot_file(n1+1:n2,4)/1000,plot_file(n1+1:n2,1),'k:','LineWidth',2); hold on;

% Plot storage media temperatures
Smin = min(plot_file(1:n2,4))/1000;
Smax = max(plot_file(1:n2,4))/1000;
Spoints = linspace(Smin-0.00*(Smax-Smin),Smax+0.00*(Smax-Smin),2);
Tlevels(1,1:2) = T0;
Tlevels(2,1:2) = HT.B(2).T;
Tlevels(3,1:2) = HT.A(1).T;
Tlevels(6,1:2) = CT.B(2).T;
Tlevels(7,1:2) = CT.A(1).T;
plot(Spoints,Tlevels(1,:),'k--','LineWidth',1.0);
plot(Spoints,Tlevels(2,:),'r-','LineWidth',1.0);
plot(Spoints,Tlevels(3,:),'r-','LineWidth',1.0);
plot(Spoints,Tlevels(6,:),'b-','LineWidth',1.0);
plot(Spoints,Tlevels(7,:),'b-','LineWidth',1.0);

% Plot points
for int = 1:stages_ch
    pl1 = plot([gas.state(1,int).s]/1000,[gas.state(1,int).T],'k-o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5); hold on;
end
for int = 1:stages_dis
    pl2 = plot([gas.state(2,int).s]/1000,[gas.state(2,int).T],'k:s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5); hold on;
end
xlabel('Specific Entropy [kJ/kg.K]');
ylabel('Temperature [K]');
%axis('auto');
ylim([0 1000])
yticks(0:100:1000)
chi = 0;%RTeff(gas,t_ch,t_dis);
title(strcat('$$\eta$$',sprintf('=%.2f',eta),'  $$\epsilon$$',sprintf('=%.3f',eff),'  $$f_p$$',sprintf('=%.3f',ploss), '  $$\chi$$',sprintf('=%.2f',chi)));
%title(strcat('$$\eta$$',sprintf('=%.2f',eta),'  $$\epsilon$$',sprintf('=%.3f',eff),'  $$f_p$$',sprintf('=%.3f',ploss),...
%    '  $$\mathrm{WR}_{dis} $$',sprintf('=%.1f',WR_dis), '  $$\chi$$',sprintf('=%.2f',chi)));
grid off;
lgn = legend([pl1 pl2],'charge','discharge','Location','best'); %lgn.FontSize = 12;
hold off;