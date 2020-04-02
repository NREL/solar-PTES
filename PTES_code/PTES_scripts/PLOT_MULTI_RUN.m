% Set properties for plots
set_graphics
fignum = 20;
save_multi_figs = 0;

% Load Multi_run_var file. This contains information about the variables
% being changed along the multi-run calls
load('./Outputs/Multi_run/Multi_run_var.mat')

%{
% Load Multi_run_var file. This contains information about the variables
% being changed along the multi-run calls
%ID0 = fopen('./Outputs/Multi_run/Multi_run_var.txt','r');

% Obtain information on Vpnt (points along each curve)
Vpnt = textscan(fgetl(ID0),'%s'); Vpnt = Vpnt{1,1}; % Variable name
Npnt = textscan(fgetl(ID0),'%d'); Npnt = Npnt{1,1}; % Number of points
Apnt = zeros(1,Npnt); % Values
for ipnt=1:Npnt
    Anum = textscan(fgetl(ID0),'%f'); Apnt(ipnt) = Anum{1,1};
end
fgetl(ID0);

% Obtain information on Vcrv (variable changed between curves)
Vcrv = textscan(fgetl(ID0),'%s'); Vcrv = Vcrv{1,1}; % Variable name
Ncrv = textscan(fgetl(ID0),'%f'); Ncrv = Ncrv{1,1}; % Number of curves
Acrv = zeros(1,Ncrv); % Values
for icrv=1:Ncrv
    Anum = textscan(fgetl(ID0),'%f'); Acrv(icrv) = Anum{1,1};
end
%}

% Generate string with Vpnt legend
switch Vpnt
    case 'PRch'
        Tpnt = ' $$ \mathrm{PR_{ch}} $$'; 
        Upnt = ' ';
    case 'TC_dis0'
        Tpnt = '$$ T_3 $$';
        Upnt = ' [K]';
    case 'TH_dis0'
        Tpnt = '$$ T_{\mathrm{bot}} $$';
        Upnt = ' [$$^{\circ}$$C]';
        Apnt = Apnt - 273.15;
    case 'Ran_TbotC'
        Tpnt = '$$ T_{\mathrm{condenser}} $$';
        Upnt = ' [$$^{\circ}$$C]';
        Apnt = Apnt - 273.15;
    case 'eff'
        Tpnt = '$$ \epsilon $$';
        Upnt = ' ';
    case 'eta'
        Tpnt = '$$ \eta $$';
        Upnt = ' ';
    otherwise
        error('not implemented')
end
Lpnt = strcat(Tpnt,Upnt);

% Generate strings with Vcrv legends
switch Vcrv
    case 'PRch'
        Tcrv = ' $$ \mathrm{PR_{ch}} $$'; 
        Ucrv = ' ';
    case 'TC_dis0'
        Tcrv = '$$ T_3 $$';
        Ucrv = ' [K]';
    case 'TH_dis0'
        Tcrv = '$$ T__{\mathrm{bot}} $$';
        Ucrv = ' [$$^{\circ}$$C]';
        Acrv = Acrv - 273.15;
    case 'Ran_TbotC'
        Tcrv = '$$ T__{\mathrm{condenser}} $$';
        Ucrv = ' [$$^{\circ}$$C]';
        Acrv = Acrv - 273.15;
    case 'eff'
        Tcrv = '$$ \epsilon $$';
        Ucrv = ' ';
    case 'eta'
        Tcrv = '$$ \eta $$';
        Ucrv = ' ';
    case 'Ne_ch'
        Tcrv = '$$ \mathrm{N_{stages}} $$';
        Ucrv = ' ';
    otherwise
        error('not implemented')
end
Lcrv = cell(1,Ncrv);
for icrv=1:Ncrv
    Lcrv{icrv} = strcat(Tcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv);
end

% EXTRACT DATA INTO ARRAYS
chi_mat   = zeros(Npnt,Ncrv);
WL_1_mat  = zeros(Npnt,Ncrv);
WL_2_mat  = zeros(Npnt,Ncrv);
WL_3_mat  = zeros(Npnt,Ncrv);
WL_4_mat  = zeros(Npnt,Ncrv);
WL_5_mat  = zeros(Npnt,Ncrv);
WL_6_mat  = zeros(Npnt,Ncrv);
WL_7_mat  = zeros(Npnt,Ncrv);
for icrv=1:Ncrv
    for ipnt=1:Npnt
        filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',icrv,ipnt);
        S = load(filename);
        chi_mat(ipnt,icrv)  = S.chi;
        WL_1_mat(ipnt,icrv) = S.WL_comp;
        WL_2_mat(ipnt,icrv) = S.WL_exp;
        WL_3_mat(ipnt,icrv) = S.WL_hexs;
        WL_4_mat(ipnt,icrv) = S.WL_reject;
        WL_5_mat(ipnt,icrv) = S.WL_mix_liq;
        WL_6_mat(ipnt,icrv) = S.WL_mix_gas;
        WL_7_mat(ipnt,icrv) = S.WL_tanks;
    end
end

% Exergetic efficiency
figure(fignum);
set(gcf,'DefaultAxesColorOrder',[0 0 0],...
     'DefaultAxesLineStyleOrder','-s|-.s|--s|:s')
for icrv=1:Ncrv
    plot(Apnt,chi_mat(:,icrv)*100); hold on;
end
hold off;
xlabel(Lpnt)
ylabel('Exergetic efficiency [$$\%$$]')
ylim([55 60])
legend(Lcrv,'Location','Best')
grid on;


WL_mat = zeros(Npnt,7,Ncrv);
for icrv=1:Ncrv
    WL_mat(:,:,icrv) = [WL_1_mat(:,icrv),WL_2_mat(:,icrv),WL_3_mat(:,icrv),WL_4_mat(:,icrv),WL_5_mat(:,icrv),WL_6_mat(:,icrv),WL_7_mat(:,icrv)];
end

for icrv=1:Ncrv
    figure(fignum+4+icrv);
    a = area(Apnt,WL_mat(:,:,icrv)); hold on;
    a(1).FaceColor = c_dark_blue;
    a(2).FaceColor = c_pale_blue;
    a(3).FaceColor = c_pale_green;
    a(4).FaceColor = c_yellow;
    a(5).FaceColor = c_pale_orange;
    a(6).FaceColor = c_dark_orange;
    a(7).FaceColor = c_grey;
    hold off;
    xlabel(strcat(Tpnt,Upnt))
    ylabel('Lost Work [$$\%$$]')
    title(Lcrv{icrv})
    ylim([0 50])
    legend({'Compressors','Expanders','Heat exchangers','Heat in/out env.','Mixing (liquid)','Mixing (gas)','Tanks'},'Location','Best')
    grid on;
    
    %{
    % Do not show liquid_mixing loss and tank_loss bars if they are not required
    if all(all(WL_5_mat==0))
        delete(a(5))
    end
    if all(all(WL_6_mat==0))
        delete(a(6))
    end
    %}
end

switch save_multi_figs
    case 1
        formats = {'epsc'};
        save_fig(fignum,'./Outputs/exergy_eff',formats);
        for icrv=1:Ncrv
            savename = strcat('./Outputs/WL_mat_',sprintf('%d',icrv));
            save_fig(fignum+10+icrv,savename,formats);
        end
end