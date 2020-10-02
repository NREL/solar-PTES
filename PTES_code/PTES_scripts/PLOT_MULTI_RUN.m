% Set properties for plots
set_graphics
stl = {'-s','-.s','--s',':s',... % line styles
    '-^','-.^','--^',':^'};
fignum = 20;
save_multi_figs = 0;

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
    case 'Ran_Tbot0'
        Tpnt = '$$ T_{\mathrm{condenser}} $$';
        Upnt = ' [$$^{\circ}$$C]';
        Apnt = Apnt - 273.15;
    case 'Ran_ptop'
        Tpnt = '$$ p_{\mathrm{top,\; Rankine}} $$';
        Upnt = ' [bar]';
        Apnt = Apnt/1e5;
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
    case 'Wdis_req'
        Tpnt = 'Discharging power output';
        Upnt = ' MW$$_e$$';
    case 'stH'
        Tpnt = 'Discharging duration';
        Upnt = ' h';
    case 'unbalanced'
        Tpnt = 'Charge duration / discharge duration' ;
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
    case 'Ran_Tbot0'
        Tcrv = '$$ T__{\mathrm{condenser}} $$';
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
    case 'Wdis_req'
        Tcrv = 'Discharging power output';
        Ucrv = ' MW$$_e$$';
    case 'stH'
        Tcrv = 'Discharging duration';
        Ucrv = ' h';
    otherwise
        error('not implemented')
end
Lcrv = cell(1,Ncrv);
for icrv=1:Ncrv
    Lcrv{icrv} = strcat(Tcrv,sprintf(' = %.3g ',Acrv(icrv)),Ucrv);
end

% EXTRACT DATA INTO ARRAYS
chi_mat      = var_extract('chi_PTES',Npnt,Ncrv);
chi_para_mat = var_extract('chi_PTES_para',Npnt,Ncrv);
lcos_mat     = var_extract('lcosM',Npnt,Ncrv);
capcost_mat  = var_extract('cap_costM',Npnt,Ncrv);
Wdis_mat = var_extract('E_net_dis',Npnt,Ncrv);
tdis_mat = var_extract('t_dis',Npnt,Ncrv);
rhoE_mat = var_extract('rhoE',Npnt,Ncrv);

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
HEeff_mat   = var_extract('HEeff',  Npnt,Ncrv);
HEeffRC_mat = var_extract('HEeffRC',Npnt,Ncrv);
HEeffNC_mat = var_extract('HEeffNC',Npnt,Ncrv);
t_chg_mat   = var_extract('t_chg',  Npnt,Ncrv);
t_dis_mat   = var_extract('t_dis',  Npnt,Ncrv);
tdRC_mat    = var_extract('t_disRC',Npnt,Ncrv);
tdNC_mat    = var_extract('t_disNC',Npnt,Ncrv);
COP_mat     = var_extract('COP', Npnt,Ncrv);
HEexergy_eff_mat = var_extract('HEexergy_eff',Npnt,Ncrv);
HEexergy_effRC_mat = var_extract('HEexergy_effRC',Npnt,Ncrv);
HEexergy_effNC_mat = var_extract('HEexergy_effNC',Npnt,Ncrv);
HPexergy_eff_mat = var_extract('HPexergy_eff',Npnt,Ncrv);


% Exergetic efficiency
figure(fignum);
for icrv=1:Ncrv
    plot(Apnt,chi_mat(:,icrv)*100,stl{icrv}); hold on;
end
for icrv=1:Ncrv
    plot(Apnt,chi_para_mat(:,icrv)*100,stl{icrv}); hold on;
end
hold off;
title('No cold store')
xlabel([Tpnt,Upnt])
ylabel('Roundtrip efficiency [$$\%$$]')
ylim([52 62])
legend([Lcrv(:)',{[Lcrv{:},', + parasitics']}],'Location','Best')
%grid on;

% LCOS
figure(fignum+1);
for icrv=1:Ncrv
    plot(Apnt,lcos_mat(:,icrv),stl{icrv}); hold on;
end
hold off;
xlabel([Tpnt,Upnt])
ylabel('LCOS [$$\$$$/kWh]')
%ylim([0 0.5])
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
%Exergetic efficiency heat pump and COP
figure(fignum+4);
yyaxis left
for icrv=1:Ncrv
    plot(Apnt,HPexergy_eff_mat(:,icrv)*100); hold on;
end
hold off;
xlabel(Lpnt)
ylabel('Exergetic efficiency heat pump  [$$\%$$]')
%ylim([65 80])
yyaxis right
L = cell(1,Ncrv);
for icrv=1:Ncrv
    plot(Apnt,COP_mat(:,icrv)); hold on;
    L{icrv} = Lcrv{icrv};
end
hold off;
ylabel('COP')
ylim([1.0 1.8])
legend(L,'Location','Best')
grid on;
%}

%%{
%Exergetic efficiency heat pump and COP
figure(fignum+4);
subplot(1,2,1);
for icrv=1:Ncrv
    plot(Apnt,HPexergy_eff_mat(:,icrv)*100,stl{icrv}); hold on;
end
hold off;
title('Brayton heat pump')
xlabel(Lpnt)
ylabel('Exergy efficiency [$$\%$$]')
%ylim([65 80])
xlim([-Inf Inf])
legend(Lcrv,'Location','South')
subplot(1,2,2);
for icrv=1:Ncrv
    plot(Apnt,COP_mat(:,icrv),stl{icrv}); hold on;
end
hold off;
title('Brayton heat pump')
xlabel(Lpnt)
ylabel('COP')
xlim([-Inf Inf])
%ylim([1.0 1.8])
legend(Lcrv,'Location','South')
%}

%{
% Exergetic efficiency and heat engine efficiency Rankine cycle
figure(fignum+5);
yyaxis left
plot(Apnt,HEexergy_effRC_mat(:,1)*100,stl{1}); hold on;
plot(Apnt,HEexergy_effNC_mat(:,1)*100,stl{2}); hold on;
hold off;
xlabel([Tpnt,Upnt])
ylabel('Exergetic efficiency heat engine [$$\%$$]')
L = cell(1,2);
L{1} = 'Using cold tanks';
L{2} = 'No cold tanks';
ylim([60 76])
yyaxis right
plot(Apnt,HEeffRC_mat(:,1)*100,stl{1}); hold on;
plot(Apnt,HEeffNC_mat(:,1)*100,stl{2}); hold on;
hold off;
ylabel('Heat engine efficiency [$$\%$$]')
ylim([36 50])
legend(L,'Location','South')
grid on;
%}

%%{
% Exergetic efficiency and heat engine efficiency Rankine cycle
figure(fignum+5);
subplot(1,2,1);
plot(Apnt,HEexergy_eff_mat(:,1)*100,stl{1}); hold on;
hold off;
title('Rankine cycle')
xlabel([Tpnt,Upnt])
ylabel('Exergy efficiency[$$\%$$]')
%ylim([66,76])
%L = cell(1,2);
%L{1} = 'Using cold tanks';
%L{2} = 'No cold tanks';
%legend(L,'Location','South')
subplot(1,2,2);
plot(Apnt,HEeff_mat(:,1)*100,stl{1}); hold on;
hold off;
title('Rankine cycle')
xlabel([Tpnt,Upnt])
ylabel('Efficiency [$$\%$$]')
%ylim([36,44])
%}

% Energy density
figure(fignum+6);
for icrv=1:Ncrv
    plot(Apnt,rhoE_mat(:,icrv),stl{icrv}); hold on;
end
hold off;
xlabel([Tpnt,Upnt])
ylabel('Energy density [$$\mathrm{kWh/m^3}$$]')
%ylim([45 60])
legend(Lcrv,'Location','Best')
%grid on;

% Cooling duration
figure(fignum+7);
for icrv=1:Ncrv
    plot(Apnt,tdRC_mat(:,icrv),stl{icrv}); hold on;
end
hold off;
xlabel([Tpnt,Upnt])
ylabel('Cooling duration [h]')
%ylim([45 60])
legend(Lcrv,'Location','Best')
%grid on;

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
        formats = {'fig','epsc'};
        save_fig(fignum,'./Results/Roundtrip_eff',formats);
        save_fig(fignum+4,'./Results/Exergy_heat_pump',formats,[20 10]);
        save_fig(fignum+5,'./Results/Exergy_heat_engine',formats,[20 10]);
        save_fig(fignum+6,'./Results/Energy_density',formats);
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