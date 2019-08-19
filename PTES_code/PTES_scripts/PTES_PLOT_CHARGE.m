% Plot PTES charge cycle
figure(1);
plot_file = load ('./Outputs/Plot_file.txt');
n1 = stages_ch*num;
plot(plot_file(1:n1,4)/1000,plot_file(1:n1,1),'k-','LineWidth',2); hold on;

% Plot storage media temperatures
Spoints = linspace(min(plot_file(1:n1,4))/1000,max(plot_file(1:n1,4))/1000,2);
Tlevels(1,1:2) = T0;
Tlevels(2,1:2) = HT.B(2).T;
Tlevels(3,1:2) = HT.A(1).T;
Tlevels(6,1:2) = CT.B(2).T;
Tlevels(7,1:2) = CT.A(1).T;
plot(Spoints,Tlevels(1,:),'k--','LineWidth',1.0);
plot(Spoints,Tlevels(2,:),'r-','LineWidth',2.0);
plot(Spoints,Tlevels(3,:),'r-','LineWidth',2.0);
plot(Spoints,Tlevels(6,:),'b-','LineWidth',2.0);
plot(Spoints,Tlevels(7,:),'b-','LineWidth',2.0);
% Plot points
for int = 1:stages_ch
    pl1 = plot([gas.state(1,int).s]/1000,[gas.state(1,int).T],'k-o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7); hold on;
end
xlabel('Specific Entropy [kJ/kg.K]');
ylabel('Temperature [K]');
axis('auto');
switch mode
    case 0 % PTES
        title(strcat('$$\eta$$',sprintf('=%.2f',eta),'  $$\epsilon$$',sprintf('=%.3f',eff),...
            '  $$f_p$$',sprintf('=%.3f',ploss)));
    case 1 % Heat pump only
        title(strcat('$$\eta$$',sprintf('=%.2f',eta),'  $$\epsilon$$',sprintf('=%.3f',eff),...
            '  $$f_p$$',sprintf('=%.3f',ploss),sprintf('  COP=%.2f',COP)));
    case 2 % Heat engine only
        error('not implemented')
end

grid off;
hold off;