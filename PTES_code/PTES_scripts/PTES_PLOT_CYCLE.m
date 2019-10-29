% Plot T-s diagram
plot_file = load ('./Outputs/Plot_file.txt');
figure(1);

Celcius = 1;
switch Celcius
    case 0
        K_C = 0;
    case 1
        K_C = -273.15;
end

switch Load.mode
    case {0,4} % PTES
        n1 = stages_ch*num;
        n2 = (stages_ch + stages_dis)*num;
        
        % Plot curves
        plot(plot_file(1:n1,4)/1000,plot_file(1:n1,1)+K_C,'k-','LineWidth',2); hold on;
        plot(plot_file(n1+1:n2,4)/1000,plot_file(n1+1:n2,1)+K_C,'k:','LineWidth',2); hold on;
        
        % Set levels of storage media temperatures
        Smin = min(plot_file(1:n2,4))/1000;
        Smax = max(plot_file(1:n2,4))/1000;
        Spoints = linspace(Smin-0.00*(Smax-Smin),Smax+0.00*(Smax-Smin),2);
        Tlevels(1,1:2) = T0;
        Tlevels(2,1:2) = HT.B(2).T;
        Tlevels(3,1:2) = HT.A(1).T;
        Tlevels(6,1:2) = CT.B(2).T;
        Tlevels(7,1:2) = CT.A(1).T;
        
        % Plot levels of storage media temperatures
        plot(Spoints,Tlevels(1,:)+K_C,'k--','LineWidth',1.0);
        plot(Spoints,Tlevels(2,:)+K_C,'r-','LineWidth',1.0);
        plot(Spoints,Tlevels(3,:)+K_C,'r-','LineWidth',1.0);
        plot(Spoints,Tlevels(6,:)+K_C,'b-','LineWidth',1.0);
        plot(Spoints,Tlevels(7,:)+K_C,'b-','LineWidth',1.0);
        
    case 1 % Heat pump only
        n1 = stages_ch*num;
        
        % Plot curves
        plot(plot_file(1:n1,4)/1000,plot_file(1:n1,1)+K_C,'k-','LineWidth',2); hold on;
        
        % Set levels of storage media temperatures
        Smin = min(plot_file(1:n1,4))/1000;
        Smax = max(plot_file(1:n1,4))/1000;
        Spoints = linspace(Smin-0.00*(Smax-Smin),Smax+0.00*(Smax-Smin),2);
        Tlevels(1,1:2) = T0;
        Tlevels(2,1:2) = HT.B(2).T;
        Tlevels(3,1:2) = HT.A(1).T;
        Tlevels(6,1:2) = CT.B(2).T;
        Tlevels(7,1:2) = CT.A(1).T;
        
        % Plot levels of storage media temperatures
        plot(Spoints,Tlevels(1,:)+K_C,'k--','LineWidth',1.0);
        plot(Spoints,Tlevels(2,:)+K_C,'r-','LineWidth',1.0);
        plot(Spoints,Tlevels(3,:)+K_C,'r-','LineWidth',1.0);
        plot(Spoints,Tlevels(6,:)+K_C,'b-','LineWidth',1.0);
        plot(Spoints,Tlevels(7,:)+K_C,'b-','LineWidth',1.0);
        
    case 2
        n2 = stages_dis*num;
        
        % Plot curves
        plot(plot_file(1:n2,4)/1000,plot_file(1:n2,1)+K_C,'k-','LineWidth',2); hold on;
        
        % Set levels of storage media temperatures
        Smin = min(plot_file(1:n2,4))/1000;
        Smax = max(plot_file(1:n2,4))/1000;
        Spoints = linspace(Smin-0.00*(Smax-Smin),Smax+0.00*(Smax-Smin),2);
        Tlevels(1,1:2) = T0;
        Tlevels(2,1:2) = HT.B(2).T;
        Tlevels(3,1:2) = HT.A(2).T;
        
        % Plot levels of storage media temperatures
        plot(Spoints,Tlevels(1,:)+K_C,'k--','LineWidth',1.0);
        plot(Spoints,Tlevels(2,:)+K_C,'r-','LineWidth',1.0);
        plot(Spoints,Tlevels(3,:)+K_C,'r-','LineWidth',1.0);
end

switch Load.mode
    case {0,4} % PTES
        % Plot points
        for iL=1:Load.num
            if any(strcmp(Load.type(iL),{'chg','chgCO2'}))
                for int = 1:stages_ch
                    pl1 = plot([gas.state(iL,int).s]/1000,[gas.state(iL,int).T]+K_C,'k-o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5); hold on;
                end
                break
            end
        end
        for iL=1:Load.num
            if any(strcmp(Load.type(iL),{'dis','disCO2'}))
                for int = 1:stages_dis
                    pl2 = plot([gas.state(iL,int).s]/1000,[gas.state(iL,int).T]+K_C,'k:s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5); hold on;
                end
                break
            end
        end
        lgn = legend([pl1 pl2],'charge','discharge','Location','best');
        title(strcat('$$\eta$$',sprintf('=%.2f',eta),'  $$\epsilon$$',sprintf('=%.3f',eff),'  $$f_p$$',sprintf('=%.3f',ploss)));
%         title(strcat('$$\eta$$',sprintf('=%.2f',eta),'  $$\epsilon$$',sprintf('=%.3f',eff),'  $$f_p$$',sprintf('=%.3f',ploss),'  $$\chi$$',sprintf('=%.2f',chi_PTES)));
%         
%     case 1 % Heat pump only
%         % Plot points
%         for int = 1:stages_ch
%             pl1 = plot([gas.state(1,int).s]/1000,[gas.state(1,int).T]+K_C,'k-o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5); hold on;
%         end
%         lgn = legend(pl1,'charge','Location','best');
%         title(strcat('$$\eta$$',sprintf('=%.2f',eta),'  $$\epsilon$$',sprintf('=%.3f',eff),'  $$f_p$$',sprintf('=%.3f',ploss),'  $$\chi_{\mathrm{hot}}$$',sprintf('=%.2f',chi_hot),sprintf('  COP=%.2f',COP)));
%         
%     case 2 % Heat engine only
%         % Plot points
%         for int = 1:stages_dis
%             pl2 = plot([gas.state(2,int).s]/1000,[gas.state(2,int).T]+K_C,'k-s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5); hold on;
%         end
%         lgn = legend(pl2,'discharge','Location','best');
%         title(strcat('$$\eta$$',sprintf('=%.2f',eta),'  $$\epsilon$$',sprintf('=%.3f',eff),'  $$f_p$$',sprintf('=%.3f',ploss),'  $$\chi$$',sprintf('=%.2f',chi_tot),'  EFF',sprintf('=%.2f',EFF)));
end

xlabel('Specific Entropy [kJ/kg.K]');
switch Celcius
    case 0
        ylabel('Temperature [K]');
        ylim([0 1000])
        yticks(0:100:1000)
    case 1
        ylabel('Temperature [$$^{\circ}$$C]');
        ylim([-150 700])
        %yticks(0:100:1000)
end
set(gcf,'color','w')
grid off; hold off;