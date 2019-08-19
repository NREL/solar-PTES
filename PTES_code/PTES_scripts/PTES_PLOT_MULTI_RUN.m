% fclose all;

% Set paths
addpath('./Classes/')
addpath('./Generic/')
addpath('./Functions/')
addpath('./PTES_scripts/')

% Set properties for plots
set_graphics

fignum = 10;
param  = load('./Outputs/Multi_run.txt');

ID0 = fopen('./Outputs/Multi_run_var.txt','r');
Vpnt = textscan(fgetl(ID0),'%s'); Vpnt = Vpnt{1,1};
Npnt = textscan(fgetl(ID0),'%f'); Npnt = Npnt{1,1};
Apnt = zeros(1,Npnt);
for ipnt=1:Npnt
    Anum = textscan(fgetl(ID0),'%f'); Apnt(ipnt) = Anum{1,1};
end
fgetl(ID0);

Vcrv = textscan(fgetl(ID0),'%s'); Vcrv = Vcrv{1,1};
Ncrv = textscan(fgetl(ID0),'%f'); Ncrv = Ncrv{1,1};
Acrv = zeros(1,Ncrv);
for icrv=1:Ncrv
    Anum = textscan(fgetl(ID0),'%f'); Acrv(icrv) = Anum{1,1};
end

if strcmp(Vpnt,'PRch')
    Lpnt = ' $$ \mathrm{PR_{ch}} $$'; 
    Upnt = ' ';
elseif strcmp(Vpnt,'TC_0')
    Lpnt = '$$ T_3 $$';
    Upnt = ' [K]';
elseif strcmp(Vpnt,'TH_0')
    Lpnt = '$$ T_{\mathrm{bot}} $$';
    Upnt = ' [$$^{\circ}$$C]';
    Apnt = Apnt - 273.15;
elseif strcmp(Vpnt,'eff')
    Lpnt = '$$ \epsilon $$';
    Upnt = ' ';
elseif strcmp(Vpnt,'eta')
    Lpnt = '$$ \eta $$';
    Upnt = ' ';
else
    error('not implemented')
end

if strcmp(Vcrv,'PRch')
    Lcrv = ' $$ \mathrm{PR_{ch}} $$'; 
    Ucrv = ' ';
elseif strcmp(Vcrv,'TC_0')
    Lcrv = '$$ T_3 $$';
    Ucrv = ' [K]';
elseif strcmp(Vcrv,'TH_0')
    Lcrv = '$$ T__{\mathrm{bot}} $$';
    Ucrv = ' [$$^{\circ}$$C]';
    Acrv = Acrv - 273.15;
elseif strcmp(Vcrv,'eff')
    Lcrv = '$$ \epsilon $$';
    Ucrv = ' ';
elseif strcmp(Vcrv,'eta')
    Lcrv = '$$ \eta $$';
    Ucrv = ' ';
elseif strcmp(Vcrv,'Ne_ch')
    Lcrv = '$$ \mathrm{N_{stages}} $$';
    Ucrv = ' ';
else
    error('not implemented')
end

% DO NOT PLOT A SPECIFIC POINT?
%OKpnts = (Apnt<420 | Apnt>460); %remove according to condition
OKpnts = ~isnan(Apnt); %keep all

% EXTRACT DATA INTO ARRAYS
chi_mat   = zeros(Npnt,Ncrv);
COP_mat   = zeros(Npnt,Ncrv);
EFF_mat   = zeros(Npnt,Ncrv);
Tmax_mat  = zeros(Npnt,Ncrv);
Tmin_mat  = zeros(Npnt,Ncrv);
WL_1_mat  = zeros(Npnt,Ncrv);
WL_2_mat  = zeros(Npnt,Ncrv);
WL_3_mat  = zeros(Npnt,Ncrv);
WL_4_mat  = zeros(Npnt,Ncrv);
WL_5_mat  = zeros(Npnt,Ncrv);
WL_6_mat  = zeros(Npnt,Ncrv);
L = cell(1,Ncrv);
for icrv=1:Ncrv
    for ipnt=1:Npnt
        chi_mat(ipnt,icrv)   = param((icrv-1)*Npnt + ipnt,5);
        COP_mat(ipnt,icrv)   = param((icrv-1)*Npnt + ipnt,6);
        EFF_mat(ipnt,icrv)   = param((icrv-1)*Npnt + ipnt,7);
        Tmax_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,8);
        Tmin_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,9);
        WL_1_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,13);
        WL_2_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,14);
        WL_3_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,15);
        WL_4_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,16);
        WL_5_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,17);
        WL_6_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,18);
    end
    L{icrv} = strcat(Lcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv);
end

% Exergetic efficiency
figure(fignum);
set(gcf,'DefaultAxesColorOrder',[0 0 0],...
     'DefaultAxesLineStyleOrder','-|-.|--')
for icrv=1:Ncrv
    plot(Apnt(OKpnts),chi_mat(OKpnts,icrv)*100); hold on;
end
hold off;
xlabel(strcat(Lpnt,Upnt))
ylabel('Exergetic efficiency [$$\%$$]')
%ylim([65 85])
ylim([40 70])
legend(L,'Location','Best')
grid on;

% COP
figure(fignum+1);
set(gcf,'DefaultAxesColorOrder',[0 0 0],...
     'DefaultAxesLineStyleOrder','-|-.|--')
for icrv=1:Ncrv
    plot(Apnt(OKpnts),COP_mat(OKpnts,icrv)); hold on;
end
hold off;
xlabel(strcat(Lpnt,Upnt))
ylabel('COP')
ylim([1 1.5])
legend(L,'Location','Best')
grid on;

% Exergetic efficiency and COP
figure(fignum+2);
L2 = cell(1,2*Ncrv);
yyaxis left
for icrv=1:Ncrv
    plot(Apnt(OKpnts),chi_mat(OKpnts,icrv)*100); hold on;
    L2{icrv} = strcat(Lcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv,' (Tmax)');
end
hold off;
xlabel(strcat(Lpnt,Upnt))
ylabel('Exergetic efficiency  [$$\%$$]')
ylim([65 80])
yyaxis right
for icrv=1:Ncrv
    plot(Apnt(OKpnts),COP_mat(OKpnts,icrv)); hold on;
    L2{icrv+Ncrv} = strcat(Lcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv,' (Tmin)');
end
hold off;
ylabel('COP')
ylim([1.0 1.8])
legend(L,'Location','Best')
grid on;

% Engine efficiency
figure(fignum+3);
set(gcf,'DefaultAxesColorOrder',[0 0 0],...
     'DefaultAxesLineStyleOrder','-|-.|--')
for icrv=1:Ncrv
    plot(Apnt(OKpnts),EFF_mat(OKpnts,icrv)*100); hold on;
end
hold off;
xlabel(strcat(Lpnt,Upnt))
ylabel('Engine efficiency [$$\%$$]')
legend(L,'Location','Best')
grid on;



% figure(fignum+2);
% %set(gcf,'DefaultAxesColorOrder',[0 0 0],...
% %      'DefaultAxesLineStyleOrder','-|--|:')
% %      'DefaultAxesLineStyleOrder','-|--|:|-.')
% L2 = cell(1,2*Ncrv);
% yyaxis left
% for icrv=1:Ncrv
%     plot(Apnt(OKpnts),Tmax_mat(OKpnts,icrv)); hold on;
%     L2{icrv} = strcat(Lcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv,' (Tmax)');
% end
% hold off;
% xlabel(strcat(Lpnt,Upnt))
% ylabel('$$ T_{\mathrm{max}} $$ [K]')
% ylim([650 950])
% yyaxis right
% for icrv=1:Ncrv
%     plot(Apnt(OKpnts),Tmin_mat(OKpnts,icrv)); hold on;
%     L2{icrv+Ncrv} = strcat(Lcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv,' (Tmin)');
% end
% hold off;
% ylabel('$$ T_{\mathrm{min}} $$ [K]')
% ylim([200 300])
% legend(L,'Location','Best')
% grid on;

% figure(fignum+3);
% %set(gcf,'DefaultAxesColorOrder',[0 0 0],...
% %      'DefaultAxesLineStyleOrder','-|--|:')
% %      'DefaultAxesLineStyleOrder','-|--|:|-.')
% for icrv=1:Ncrv
%     plot(Apnt(OKpnts),Tmin_mat(OKpnts,icrv)); hold on;
% end
% hold off;
% xlabel(strcat(Lpnt,Upnt))
% ylabel('$$ T_4 $$ [K]')
% legend(L,'Location','Best')
% grid on;



WL_mat = zeros(Npnt,6,Ncrv);
for icrv=1:Ncrv
    WL_mat(:,:,icrv) = [WL_1_mat(:,icrv),WL_2_mat(:,icrv),WL_3_mat(:,icrv),WL_4_mat(:,icrv),WL_5_mat(:,icrv),WL_6_mat(:,icrv)];
end

for icrv=1:Ncrv
    figure(fignum+4+icrv);
    a = area(Apnt(OKpnts),WL_mat(OKpnts,:,icrv)); hold on;
    a(2).FaceColor = c_pale_blue;
    a(3).FaceColor = c_pale_green;
    a(4).FaceColor = 'yellow';
    a(5).FaceColor = c_dark_orange;
    a(6).FaceColor = c_grey;
    hold off;
    xlabel(strcat(Lpnt,Upnt))
    ylabel('Lost Work [$$\%$$]')
    title(L{icrv})
    ylim([0 40])
    legend({'Compressors','Expanders','Heat exchangers','Heat rejected/absorbed','Mixing (liquid)','Tanks'},'Location','NorthWest')
    grid on;
    
    % Do not show liquid_mixing loss and tank_loss bars if they are not required
    if all(all(WL_5_mat==0))
        delete(a(5))
    end
    if all(all(WL_6_mat==0))
        delete(a(6))
    end
end

save_figs = 0;
if save_figs
    save_fig(fignum,'./Outputs/exergy_eff',0,0,0);
    save_fig(fignum+1,'./Outputs/COP',0,0,0);
    save_fig(fignum+2,'./Outputs/chi_COP',0,0,0);
    save_fig(fignum+3,'./Outputs/EFF',0,0,0);
    for icrv=1:Ncrv
        savename = strcat('./Outputs/WL_mat_',sprintf('%d',icrv));
        save_fig(fignum+4+icrv,savename,0,0,0);
    end
end