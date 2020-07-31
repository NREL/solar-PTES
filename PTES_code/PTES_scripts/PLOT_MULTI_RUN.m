% Set properties for plots
set_graphics
stl = {'-s','-.s','--s',':s',... % line styles
    '-^','-.^','--^',':^'};
fignum = 20;
save_multi_figs = 1;

% Load Multi_run_var file. This contains information about the variables
% being changed along the multi-run calls
load('./Outputs/Multi_run/Multi_run_var.mat')

% Generate string with Vpnt legend
switch Vpnt
    case 'PRch'
        Tpnt = ' $$ \mathrm{PR_{ch}} $$'; 
        Upnt = ' ';
    case 'pmax'
        Tpnt = ' $$ \mathrm{P_{max}} $$'; 
        Upnt = ' [bar]';
    case 'PRr'
        Tpnt = ' $$ \mathrm{PR_{r}} $$'; 
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
    case 'ploss'
        Tpnt = 'Fractional pressure loss';
        Upnt = ' ';
    case 'eta'
        Tpnt = '$$ \eta $$';
        Upnt = ' ';
    case 'mdot_off'
        Tpnt = '$$ \dot{m} / \dot{m}_0 $$';
        Upnt = ' ';
    case 'T0_off'
        Tpnt = 'Change in ambient temperature';
        Upnt = 'K' ;
    otherwise
        error('not implemented')
end
Lpnt = strcat(Tpnt,Upnt);

% Generate strings with Vcrv legends
switch Vcrv
    case 'PRch'
        Tcrv = ' $$ \mathrm{PR_{ch}} $$'; 
        Ucrv = ' ';
    case 'pmax'
        Tcrv = ' $$ \mathrm{P_{max}} $$'; 
        Ucrv = ' [bar]';
    case 'PRr'
        Tcrv = ' $$ \mathrm{PR_{r}} $$'; 
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
    case 'ploss'
        Tcrv = 'Fractional pressure loss';
        Ucrv = ' ';
    case 'eta'
        Tcrv = '$$ \eta $$';
        Ucrv = ' ';
    case 'Ne_ch'
        Tcrv = '$$ \mathrm{N_{stages\;heat\;pump}} $$';
        Ucrv = ' ';
    case 'mdot_off'
        Tcrv = '$$ \dot{m} / \dot{m}_0 $$';
        Ucrv = ' ';
    case 'T0_off'
        Tcrv = 'Change in ambient temperature';
        Ucrv = 'K' ;
    otherwise
        error('not implemented')
end
Lcrv = cell(1,Ncrv);
for icrv=1:Ncrv
    Lcrv{icrv} = strcat(Tcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv);
end

% EXTRACT DATA INTO ARRAYS
chi_mat  = var_extract('chi_PTES_para',Npnt,Ncrv);
lcos_mat = var_extract('lcosM',Npnt,Ncrv);
capcost_mat = var_extract('cap_costM',Npnt,Ncrv);
Wdis_mat = var_extract('W_net_dis',Npnt,Ncrv);
tdis_mat = var_extract('t_dis',Npnt,Ncrv);

Wpow_mat = Wdis_mat ./ tdis_mat ./ 1e6 ;

mdot_mat = var_extract('mdot',Npnt,Ncrv);
Ttop_mat = var_extract('Ttop',Npnt,Ncrv);
% WL_1_mat = var_extract('WL_comp',   Npnt,Ncrv);
% WL_2_mat = var_extract('WL_exp',    Npnt,Ncrv);
% WL_3_mat = var_extract('WL_hexs',   Npnt,Ncrv);
% WL_4_mat = var_extract('WL_reject', Npnt,Ncrv);
% WL_5_mat = var_extract('WL_mix_liq',Npnt,Ncrv);
% WL_6_mat = var_extract('WL_mix_gas',Npnt,Ncrv);
% WL_7_mat = var_extract('WL_tanks',  Npnt,Ncrv);
% HEeff_mat   = var_extract('HEeff',  Npnt,Ncrv);
% HEeffRC_mat = var_extract('HEeffRC',Npnt,Ncrv);
% HEeffNC_mat = var_extract('HEeffNC',Npnt,Ncrv);
% t_dis_mat   = var_extract('t_dis',  Npnt,Ncrv);
% tdRC_mat    = var_extract('t_disRC',Npnt,Ncrv);
% tdNC_mat    = var_extract('t_disNC',Npnt,Ncrv);


% Exergetic efficiency
figure(fignum);
for icrv=1:Ncrv
    plot(Apnt,chi_mat(:,icrv)*100,stl{icrv}); hold on;
end
hold off;
xlabel([Tpnt,Upnt])
ylabel('Roundtrip efficiency [$$\%$$]')
ylim([20 80])
legend(Lcrv,'Location','Best')
grid on;

% LCOS
figure(fignum+1);
for icrv=1:Ncrv
    plot(Apnt,lcos_mat(:,icrv),stl{icrv}); hold on;
end
hold off;
xlabel([Tpnt,Upnt])
ylabel('LCOS [$$\$$$/kWh]')
ylim([0 0.5])
legend(Lcrv,'Location','Best')
grid on;

% Capital cost
figure(fignum+2);
for icrv=1:Ncrv
    plot(Apnt,capcost_mat(:,icrv),stl{icrv}); hold on;
end
hold off;
xlabel([Tpnt,Upnt])
ylabel('Capital cost [$$\$$$]')
ylim([0 5e7])
legend(Lcrv,'Location','Best')
grid on;

% Power output during discharging
figure(fignum+3);
for icrv=1:Ncrv
    plot(Apnt,Wpow_mat(:,icrv),stl{icrv}); hold on;
end
hold off;
xlabel([Tpnt,Upnt])
ylabel('Discharging power output [$$MW_e$$]')
ylim([0 200])
legend(Lcrv,'Location','Best')
grid on;



%{
% Efficiency Rankine cycle
figure(fignum+1);
yyaxis left
plot(Apnt,HEeffRC_mat(:,1)*100,stl{1}); hold on;
plot(Apnt,HEeffNC_mat(:,1)*100,stl{2}); hold on;
hold off;
xlabel([Tpnt,' (when using cold tanks)',Upnt])
ylabel('Heat engine efficiency [$$\%$$]')
L = cell(1,2+Ncrv);
L{1} = 'Using cold tanks';
L{2} = 'No cold tanks';
ylim([41 47])
yyaxis right
for icrv=1:Ncrv
    plot(Apnt,tdRC_mat(:,icrv)/3600,stl{4+icrv}); hold on;
    L{2+icrv} = Lcrv{icrv};
end
hold off;
ylabel('Discharge time of cold tanks [h]')
ylim([0 8])
legend(L,'Location','North')
grid on;
%}

%{
%Exergetic efficiency and COP
figure(fignum+2);
yyaxis left
for icrv=1:Ncrv
    plot(Apnt(OKpnts),chi_mat(OKpnts,icrv)*100); hold on;
end
hold off;
xlabel(strcat(Lpnt,Upnt))
ylabel('Exergetic efficiency  [$$\%$$]')
%ylim([65 80])
yyaxis right
for icrv=1:Ncrv
    plot(Apnt(OKpnts),COP_mat(OKpnts,icrv)); hold on;
end
hold off;
ylabel('COP')
ylim([1.0 1.8])
legend(L,'Location','Best')
grid on;
%}
%{
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
    xlabel(Lpnt)
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
%}
switch save_multi_figs
    case 1
        formats = {'epsc'};
        save_fig(fignum,'./Outputs/exergy_eff',formats);
        save_fig(fignum+1,'./Outputs/HEeff_and_time',formats);
        %{
        for icrv=1:Ncrv
            savename = strcat('./Outputs/WL_mat_',sprintf('%d',icrv));
            save_fig(fignum+10+icrv,savename,formats);
        end
        %}
end


%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_mat = var_extract(var_name,Npnt,Ncrv)

% Create a matrix by extracting and arranging the values of the variable
% 'var_name' from the different files inside the Multi_run folder.
var_mat   = zeros(Npnt,Ncrv);
for icrv=1:Ncrv
    for ipnt=1:Npnt
        filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',icrv,ipnt);
        S = load(filename,var_name);
        var_mat(ipnt,icrv)  = S.(var_name);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%