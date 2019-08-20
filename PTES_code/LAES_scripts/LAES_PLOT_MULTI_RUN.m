% fclose all;
% clear;
% addpath('/home/pf298/Desktop/Matlab and thermodynamic analysis/useful generic functions/')
% set_graphics
% 
% % Automatically save the figures?
% savefig = 1;

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

if strcmp(Vpnt,'PR')
    Lpnt = ' $$ \mathrm{PR} $$'; 
    Upnt = ' ';
elseif strcmp(Vpnt,'TC_0')
    Lpnt = '$$ T_3 $$';
    Upnt = ' [K]';
elseif strcmp(Vpnt,'TH_0')
    Lpnt = '$$ T_1 $$';
    Upnt = ' [K]';
elseif strcmp(Vpnt,'eff')
    Lpnt = '$$ \epsilon $$';
    Upnt = ' ';
elseif strcmp(Vpnt,'eta')
    Lpnt = '$$ \eta $$';
    Upnt = ' ';
else
    error('not implemented')
end

if strcmp(Vcrv,'PR')
    Lcrv = ' $$ \mathrm{PR} $$'; 
    Ucrv = ' ';
elseif strcmp(Vcrv,'TC_0')
    Lcrv = '$$ T_3 $$';
    Ucrv = ' [K]';
elseif strcmp(Vcrv,'TH_0')
    Lcrv = '$$ T_1 $$';
    Ucrv = ' [K]';
elseif strcmp(Vcrv,'eff')
    Lcrv = '$$ \epsilon $$';
    Ucrv = ' ';
elseif strcmp(Vcrv,'eta')
    Lcrv = '$$ \eta $$';
    Ucrv = ' ';
else
    error('not implemented')
end

% DO NOT PLOT A SPECIFIC POINT?
%OKpnts = (Apnt<420 | Apnt>460); %remove according to condition
OKpnts = ~isnan(Apnt); %keep all

% EXTRACT DATA INTO ARRAYS
chi_mat   = zeros(Npnt,Ncrv);
Tmax_mat  = zeros(Npnt,Ncrv);
Tmin_mat  = zeros(Npnt,Ncrv);
WL_1_mat  = zeros(Npnt,Ncrv);
WL_2_mat  = zeros(Npnt,Ncrv);
WL_3_mat  = zeros(Npnt,Ncrv);
WL_4_mat  = zeros(Npnt,Ncrv);
WL_5_mat  = zeros(Npnt,Ncrv);
L = cell(1,Ncrv);
for icrv=1:Ncrv
    for ipnt=1:Npnt
        chi_mat(ipnt,icrv)   = param((icrv-1)*Npnt + ipnt,5);
        Tmax_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,6);
        Tmin_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,7);
        WL_1_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,11);
        WL_2_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,12);
        WL_3_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,13);
        WL_4_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,14);
        WL_5_mat(ipnt,icrv)  = param((icrv-1)*Npnt + ipnt,15);
    end
    L{icrv} = strcat(Lcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv);
end

% MAKE PLOTS
figure(fignum);
set(gcf,'DefaultAxesColorOrder',[0 0 0],...
     'DefaultAxesLineStyleOrder','--|-.|-')
% %      'DefaultAxesLineStyleOrder','-|--|-.|:')
for icrv=1:Ncrv
    plot(Apnt(OKpnts),chi_mat(OKpnts,icrv)); hold on;
end
hold off;
xlabel(strcat(Lpnt,Upnt))
ylabel('Round trip efficiency [$$\%$$]')
ylim([40 80])
%yticks(55:1:65)
legend(L,'Location','Best')
grid on;



figure(fignum+2);
%set(gcf,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-|--|:')
%      'DefaultAxesLineStyleOrder','-|--|:|-.')
L2 = cell(1,2*Ncrv);
yyaxis left
for icrv=1:Ncrv
    plot(Apnt(OKpnts),Tmax_mat(OKpnts,icrv)); hold on;
    L2{icrv} = strcat(Lcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv,' (Tmax)');
end
hold off;
xlabel(strcat(Lpnt,Upnt))
ylabel('$$ T_{\mathrm{max}} $$ [K]')
ylim([650 950])
yyaxis right
for icrv=1:Ncrv
    plot(Apnt(OKpnts),Tmin_mat(OKpnts,icrv)); hold on;
    L2{icrv+Ncrv} = strcat(Lcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv,' (Tmin)');
end
hold off;
ylabel('$$ T_{\mathrm{min}} $$ [K]')
ylim([200 300])
legend(L,'Location','Best')
grid on;

figure(fignum+3);
%set(gcf,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-|--|:')
%      'DefaultAxesLineStyleOrder','-|--|:|-.')
for icrv=1:Ncrv
    plot(Apnt(OKpnts),Tmin_mat(OKpnts,icrv)); hold on;
end
hold off;
xlabel(strcat(Lpnt,Upnt))
ylabel('$$ T_4 $$ [K]')
legend(L,'Location','Best')
grid on;



WL_mat = zeros(Npnt,5,Ncrv);
for icrv=1:Ncrv
    WL_mat(:,:,icrv) = [WL_1_mat(:,icrv),WL_2_mat(:,icrv),WL_3_mat(:,icrv),WL_4_mat(:,icrv),WL_5_mat(:,icrv)];
end

for icrv=1:Ncrv
    figure(fignum+3+icrv);
    a = area(Apnt(OKpnts),WL_mat(OKpnts,:,icrv)); hold on;
    a(2).FaceColor = c_pale_green;
    a(3).FaceColor = 'yellow';
    a(4).FaceColor = c_dark_orange;
    a(5).FaceColor = c_grey;
    hold off;
    xlabel(strcat(Lpnt,Upnt))
    ylabel('Lost Work [$$\%$$]')
    title(L{icrv})
    ylim([0 70])
    legend({'Turbomachinery','Heat exchangers','Heat rejection','Mixing (liquid)','Tanks'},'Location','NorthWest')
    grid on;
    
    %Do not show Mixing liquid loss (because is non existent)
    delete(a(4))
end

if savefig
    fig_size_save(fignum,'round_trip',0,0,0);
    fig_size_save(fignum+2,'T_ranges',0,0,0);    
    for icrv=1:Ncrv
        savename = strcat('WL_mat_',sprintf('%d',icrv));
        fig_size_save(fignum+3+icrv,savename,0,0,0);
    end
end