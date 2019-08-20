% Plot T-s diagram
plot_file1 = load ('./Outputs/Plot_file1.txt');
figure(1);

n1 = stages1_ch*num;
n2 = (stages1_ch + stages1_dis)*num;

% Plot curves
plot(plot_file1(1:n1,4),plot_file1(1:n1,1),'k-','LineWidth',2); hold on;
plot(plot_file1(n1+1:n2,4),plot_file1(n1+1:n2,1),'k:','LineWidth',2);

% Set levels of storage media temperatures
Smin = min(plot_file1(1:n2,4));
Smax = max(plot_file1(1:n2,4));
Spoints = linspace(Smin-0.00*(Smax-Smin),Smax+0.00*(Smax-Smin),2);
Tlevels(1,1:2) = T0;
Tlevels(2,1:2) = HT.B(2).T;
Tlevels(3,1:2) = HT.A(1).T;
Tlevels(6,1:2) = CT1.B(2).T;
Tlevels(7,1:2) = CT1.A(1).T;
Tlevels(8,1:2) = CT2.B(2).T;
Tlevels(9,1:2) = CT2.A(1).T;

% Plot levels of storage media temperatures
plot(Spoints,Tlevels(1,:),'k--','LineWidth',1.0);
plot(Spoints,Tlevels(2,:),'r-','LineWidth',1.0);
plot(Spoints,Tlevels(3,:),'r-','LineWidth',1.0);
plot(Spoints,Tlevels(6,:),'b-','LineWidth',1.0);
plot(Spoints,Tlevels(7,:),'b-','LineWidth',1.0);
plot(Spoints,Tlevels(8,:),'b-','LineWidth',1.0);
plot(Spoints,Tlevels(9,:),'b-','LineWidth',1.0);

% Plot points
for int = 1:stages1_ch+1
    pl1 = plot([gas1.state(1,int).s],[gas1.state(1,int).T],'k-o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
end
for int = 1:stages1_dis+1
    pl2 = plot([gas1.state(2,int).s],[gas1.state(2,int).T],'k:s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
end

xlabel('Specific Entropy [kJ/kg.K]');
ylabel('Temperature [K]');
ylim([0 600])
%yticks(0:100:1000)



% Plot saturation curve
ns1 = 100; %number of plotting points for saturation curve
ns2 = 10;  %number of curves with different qualities
T_vect_c = linspace(62, 132.5, ns1); %from above freezing up to critical temperature
Q0_vect  = linspace(0.0, 0.0, ns1);
Q1_vect  = linspace(1.0, 1.0, ns1);
Qv       = linspace(0.0,1.0,ns2);

sL = CP1('QT_INPUTS',Q0_vect,T_vect_c,'S',gas1.handle);
sG = CP1('QT_INPUTS',Q1_vect,T_vect_c,'S',gas1.handle);
for i01=1:ns1
    for i02=1:ns2
        sT(i02,i01) = CP1('QT_INPUTS',Qv(i02),T_vect_c(i01),'S',gas1.handle);
    end
end
plot(sL,T_vect_c,'-k');
plot(sG,T_vect_c,'-k');
for i0=1:ns2
    plot(sT(i0,:),T_vect_c,'-k','LineWidth',1);
end
hold off;

%title(strcat('$$\eta$$',sprintf('=%.2f',eta),'  $$\epsilon$$',sprintf('=%.3f',eff),'  $$f_p$$',sprintf('=%.3f',ploss),'  $$\chi$$',sprintf('=%.2f',chi_LAES)));
lgn = legend([pl1 pl2],'charge','discharge','Location','best');