% Set properties for plots
set_graphics
stl = {'-s','-.s','--s',':s',... % line styles
    '-^','-.^','--^',':^'};
fignum = 20;
save_multi_figs = 0;

% Load Multi_run_var file. This contains information about the variables
% being changed along the multi-run calls
outdir = './Outputs/Multi_run/';
%outdir = '';
load([outdir,'Multi_run_var.mat'])

% Save copy of file in Outputs folder
copyfile('./PTES_scripts/PLOT_MULTI_RUN.m',[outdir,'PLOT_MULTI_RUN.m']);

% Generate string with Vpnt legend
switch Vpnt
    case 'PRch'
        Tpnt = ' $$ \mathrm{PR_{ch}} $$';
        Upnt = ' ';
    case 'pmax'
        Tpnt = ' $$ \mathrm{P_{max}} $$';
        Upnt = ' [bar]';
        Apnt = Apnt/1e5;
    case 'pmax_LA'
        Tpnt = ' $$ p_{\mathrm{max,LA}} $$';
        Upnt = ' [bar]';
        Apnt = Apnt/1e5;
    case 'PRr'
        Tpnt = ' $$ \mathrm{PR_{r}} $$';
        Upnt = ' ';
    case 'TC_dis0'
        Tpnt = '$$ T_3 $$';
        Upnt = ' [K]';
    case 'TH_dis0'
        Tpnt = '$$ T_1 $$';
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
        Tpnt = 'Top pressure';
        Upnt = ' [bar]';
        Apnt = Apnt/1e5;
    case 'eff'
        Tpnt = 'Heat exchanger effectiveness, $$ \epsilon $$';
        Upnt = ' ';
    case 'ploss'
        Tpnt = 'Fractional pressure loss';
        Upnt = ' ';
    case 'eta'
        Tpnt = '$$ \eta $$';
        Upnt = ' ';
    case 'Ne_ch'
        Tpnt = 'Heat pump expansions';
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
        Acrv = Acrv/1e5;
    case 'pmax_LA'
        Tcrv = ' $$ p_{\mathrm{max,LA}} $$';
        Ucrv = ' [bar]';
        Acrv = Acrv/1e5;
    case 'PRr'
        Tcrv = ' $$ \mathrm{PR_{r}} $$';
        Ucrv = ' ';
    case 'TC_dis0'
        Tcrv = '$$ T_3 $$';
        Ucrv = ' [K]';
    case 'TH_dis0'
        Tcrv = '$$ T_1 $$';
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
    case 'Ran_ptop'
        Tcrv = '$$ p_{\mathrm{top,\; Rankine}} $$';
        Ucrv = ' [bar]';
        Acrv = Acrv/1e5;
    case 'eff'
        Tcrv = '$$ \epsilon $$';
        Ucrv = ' ';
    case 'ploss'
        Tcrv = 'Fractional pressure loss';
        Ucrv = ' ';
    case 'eta'
        Tcrv = '$$ \eta_\mathrm{c/t} $$';
        Ucrv = ' ';
    case 'Ne_ch'
        %Tcrv = 'Heat pump expansions';
        Tcrv = 'N$_\mathrm{exp}$';
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
Load_mode    = var_extract('Load_mode',outdir,Npnt,Ncrv); Load_mode=Load_mode(1);
chi_mat      = var_extract('chi_PTES',outdir,Npnt,Ncrv);
chi_para_mat = var_extract('chi_PTES_para',outdir,Npnt,Ncrv);
lcos_mat     = var_extract('lcosM',outdir,Npnt,Ncrv);
capcost_mat  = var_extract('cap_costM',outdir,Npnt,Ncrv);
E_net_chg_mat = var_extract('E_net_chg',outdir,Npnt,Ncrv);
E_net_dis_mat = var_extract('E_net_dis',outdir,Npnt,Ncrv);
tdis_mat = var_extract('t_dis',outdir,Npnt,Ncrv);
rhoE_mat = var_extract('rhoE',outdir,Npnt,Ncrv);

Wpow_mat = E_net_dis_mat ./ tdis_mat ./ 1e6 ;

mdot_mat = var_extract('mdot',outdir,Npnt,Ncrv);
Ttop_mat = var_extract('Ttop',outdir,Npnt,Ncrv);
TLA_mat  = var_extract('TLA',outdir,Npnt,Ncrv);
Te3_mat  = var_extract('Te3',outdir,Npnt,Ncrv);
THmax_mat  = var_extract('THmax',outdir,Npnt,Ncrv);
THPmin_mat = var_extract('THPmin',outdir,Npnt,Ncrv);
HEeff_mat   = var_extract('HEeff',outdir,  Npnt,Ncrv);
HEeffRC_mat = var_extract('HEeffRC',outdir,Npnt,Ncrv);
HEeffNC_mat = var_extract('HEeffNC',outdir,Npnt,Ncrv);
t_chg_mat   = var_extract('t_chg',outdir,  Npnt,Ncrv);
t_dis_mat   = var_extract('t_dis',outdir,  Npnt,Ncrv);
tdRC_mat    = var_extract('t_disRC',outdir,Npnt,Ncrv);
tdNC_mat    = var_extract('t_disNC',outdir,Npnt,Ncrv);
COP_mat     = var_extract('COP',outdir, Npnt,Ncrv);
Qstr_W_mat  = var_extract('Qstr_W_ratio',outdir, Npnt,Ncrv);
QT_W_mat    = var_extract('QT_W_ratio',outdir, Npnt,Ncrv);
HEexergy_eff_mat = var_extract('HEexergy_eff',outdir,Npnt,Ncrv);
HEexergy_effRC_mat = var_extract('HEexergy_effRC',outdir,Npnt,Ncrv);
HEexergy_effNC_mat = var_extract('HEexergy_effNC',outdir,Npnt,Ncrv);
HPexergy_eff_mat = var_extract('HPexergy_eff',outdir,Npnt,Ncrv);
WL_matrix_mat = mat_extract('WL_matrix',outdir,Npnt,Ncrv);

switch Load_mode
    case 3
        
        % Round-trip efficiency
        figure(fignum);
        %for icrv=1:Ncrv
        %    plot(Apnt,chi_mat(:,icrv)*100,stl{icrv}); hold on;
        %end
        for icrv=1:Ncrv
            plot(Apnt,chi_para_mat(:,icrv)*100,stl{icrv}); hold on;
        end
        hold off;
        %title('No cold store')
        xlabel([Tpnt,Upnt])
        ylabel('Roundtrip efficiency [$$\%$$]')
        xlim([-Inf Inf])
        %ylim([52 62])
        legend(Lcrv(:),'Location','Best')
        %legend([Lcrv(:)',{[Lcrv{:},', + parasitics']}],'Location','Best')
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
        
        %Exergetic efficiency heat pump and COP
        %{
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
        
        % Set positions for subplots
        pos1 = [0.14,0.40,0.8,0.52];
        pos2 = [0.14,0.10,0.8,0.25];
        
        %Exergetic efficiency heat pump and COP
        %%{
        figure(fignum+4);
        subplot('Position',pos1);
        for icrv=1:Ncrv
            plot(Apnt,HPexergy_eff_mat(:,icrv)*100,stl{icrv}); hold on;
        end
        hold off;
        title('Brayton heat pump')
        ylabel('Exergy efficiency [$$\%$$]')
        xlim([-Inf Inf])
        xticklabels('');
        legend(Lcrv,'Location','Best')
        subplot('Position',pos2)
        for icrv=1:Ncrv
            plot(Apnt,COP_mat(:,icrv),stl{icrv}); hold on;
        end
        hold off;
        xlim([-Inf Inf])
        xlabel(Lpnt)
        ylabel('COP')
        %}
        
        % Exergetic efficiency and heat engine efficiency Rankine cycle
        %{
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
        subplot('Position',pos1);
        plot(Apnt,HEexergy_eff_mat(:,1)*100,stl{1}); hold on;
        hold off;
        title('Rankine cycle')
        ylabel('Exergy efficiency [$$\%$$]')
        xlim([-Inf Inf])
        xticklabels('');
        %L = cell(1,2);
        %L{1} = 'Using cold tanks';
        %L{2} = 'No cold tanks';
        %legend(L,'Location','South')
        subplot('Position',pos2);
        plot(Apnt,HEeff_mat(:,1)*100,stl{1}); hold on;
        hold off;
        xlim([-Inf Inf])
        ylim([36 41])
        xlabel([Tpnt,Upnt])
        ylabel('Efficiency [$$\%$$]')
        %}
        
        % Energy density
        figure(fignum+6);
        for icrv=1:Ncrv
            plot(Apnt,rhoE_mat(:,icrv),stl{icrv}); hold on;
        end
        hold off;
        xlabel([Tpnt,Upnt])
        ylabel('Energy density [$$\mathrm{kWh/m^3}$$]')
        xlim([-Inf Inf])
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
        
        
        % Distribution of exergetic losses
        %%{
        WL_HP_mat = zeros(Npnt,9,Ncrv);
        WL_HE_mat = zeros(Npnt,9,Ncrv);
        for icrv=1:Ncrv
            for ipnt=1:Npnt
                WL = WL_matrix_mat{ipnt,icrv}/(-E_net_chg_mat(ipnt,icrv))*100;
                WL_HP_mat(ipnt,:,icrv) = WL(1,:);
                WL_HE_mat(ipnt,:,icrv) = WL(2,:);
            end
        end        
        for icrv=1:Ncrv
            figure(fignum+10+icrv);
            %a = bar(Apnt,WL_HP_mat(:,:,icrv),'stacked'); hold on;
            a = area(Apnt,WL_HP_mat(:,:,icrv)); hold on;
            a(1).FaceColor = c_dark_blue;
            a(2).FaceColor = c_pale_blue;
            a(3).FaceColor = c_pale_orange;
            a(4).FaceColor = c_yellow;
            a(5).FaceColor = c_dark_green;
            a(6).FaceColor = c_pale_green;
            a(7).FaceColor = c_dark_orange;
            a(8).FaceColor = c_grey;
            a(9).FaceColor = c_dark_grey;
            hold off;
            xlabel(Lpnt)
            ylabel('Lost work [$$\%$$]')
            title('Exergy losses in heat pump')
            %title(Lcrv{icrv})
            %ylim([0 50])
            legend({'Compressors','Expanders','Heat exchangers','Heat in/out env.','Mixing (tanks)','Mixing (internal)','Tanks','Parasitics','Motor/Gen'},'Location','EastOutside')
            grid on;
            
            figure(fignum+20+icrv);
            %a = bar(Apnt,WL_HE_mat(:,:,icrv),'stacked'); hold on;
            a = area(Apnt,WL_HE_mat(:,:,icrv)); hold on;
            a(1).FaceColor = c_dark_blue;
            a(2).FaceColor = c_pale_blue;
            a(3).FaceColor = c_pale_orange;
            a(4).FaceColor = c_yellow;
            a(5).FaceColor = c_dark_green;
            a(6).FaceColor = c_pale_green;
            a(7).FaceColor = c_dark_orange;
            a(8).FaceColor = c_grey;
            a(9).FaceColor = c_dark_grey;
            hold off;
            xlabel(Lpnt)
            ylabel('Lost work [$$\%$$]')
            title('Exergy losses in heat engine')
            %title(Lcrv{icrv})
            %ylim([0 50])
            legend({'Compressors','Expanders','Heat exchangers','Heat in/out env.','Mixing (tanks)','Mixing (internal)','Tanks','Parasitics','Motor/Gen'},'Location','EastOutside')
            grid on;
            
            % Do not show liquid_mixing loss and tank_loss bars if they are not required
            %{
            if all(all(WL_5_mat==0))
                delete(a(5))
            end
            if all(all(WL_6_mat==0))
                delete(a(6))
            end
            %}
        end
        %}
        
        % Use to delete a bar from an existing figure (e.g. bar number 3)
        %{
        fig     = gcf;
        axObj   = fig.Children(2);
        dataObj = axObj.Children;
        delete(dataObj(3));
        %}
        
    case {22,24}
        
        % Round-trip efficiency
        %%{
        figure(fignum+2);
        for icrv=Ncrv:-1:1
            plot(Apnt,chi_para_mat(:,icrv)*100,sty{Ncrv+1-icrv}); hold on;
        end
        hold off;
        ylim([50 65])
        %ylim([55 70])
        xlabel([Tpnt,Upnt])
        ylabel('Roundtrip efficiency [$$\%$$]')
        legend(Lcrv(Ncrv:-1:1),'Location','Best')
        grid on;
        %}
        
        % Round-trip efficiency and heat-to-work ratio (in different axes)
        %{
        figure(fignum);
        yyaxis left
        for icrv=1:Ncrv
            plot(Apnt,chi_para_mat(:,icrv)*100,stl{icrv}); hold on;
        end
        hold off;
        xlabel([Tpnt,Upnt])
        ylabel('Roundtrip efficiency [$$\%$$]')
        ylim([50 58])
        yyaxis right
        for icrv=1:Ncrv
            plot(Apnt,QT_W_mat(:,icrv),stl{icrv}); hold on;
        end
        for icrv=1:Ncrv
            plot(Apnt,Qstr_W_mat(:,icrv),stl{icrv+1}); hold on;
        end
        ylabel('$$Q / W_\mathrm{net}$$')
        hold off;
        legend('$\eta_{RT}$','$Q_\mathrm{tot}/W_\mathrm{net}$','$Q_{\mathrm{str}}/W_\mathrm{net}$','Location','Best')
        %}
        
        % Round-trip efficiency and heat-to-work ratio (in different subplots)
        %{
        figure(fignum);
        pos1 = [0.12 0.55 0.77 0.40];
        pos2 = [0.12 0.12 0.77 0.38];
        subplot('Position',pos1);
        plot(Apnt,chi_para_mat*100);
        ylabel('Efficiency [$$\%$$]')
        legend('$$\eta_\mathrm{RT}$$ ','Location','Best')
        ylim([50 58])
        xticklabels({})
        subplot('Position',pos2);
        plot(Apnt,QT_W_mat,'-.'); hold on;
        plot(Apnt,Qstr_W_mat,'--'); hold off;
        xlabel([Tpnt,Upnt])
        ylabel('$$Q / W_\mathrm{net}$$')
        legend('$Q_\mathrm{tot}/W_\mathrm{net}$','$Q_{\mathrm{str}}/W_\mathrm{net}$','Location','Best')
        %}
        
        % Temperatures
        %{
        figure(fignum+1);
        subplot('Position',pos1);
        plot(Apnt,THmax_mat-273.15);
        xticklabels({})
        ylabel('Temperature [$$^\circ$$C]')
        legend('Hot stores','Location','Best')
        subplot('Position',pos2);
        plot(Apnt,TLA_mat-273.15,'-.'); hold on;
        plot(Apnt,Te3_mat-273.15,'--'); hold off;
        xlabel([Tpnt,Upnt])
        ylabel('Temperature [$$^\circ$$C]')
        legend('Liquid air (point b7)','Turbine outlet (point e3)','Location','NorthEast')
        %}

        
        %{
        formats = {'fig','svg'};
        outdir = ['~/Dropbox/Work/Publications/Working publications/',...
            'PTES-LAES encyclopedia chapter/Chapter/Figs/'];
        save_fig(fignum,[outdir,'Roundtrip_eff'],formats);
        save_fig(fignum+1,[outdir,'Temperatures'],formats);
        save_fig(fignum+2,[outdir,'Roundtrip_eff_param'],formats);
        %}
end


switch save_multi_figs
    case 1
        formats = {'fig','epsc'};
        save_fig(fignum,[outdir,'Roundtrip_eff'],formats);
        save_fig(fignum+4,[outdir,'Exergy_heat_pump'],formats);
        save_fig(fignum+5,[outdir,'Exergy_heat_engine'],formats);
        save_fig(fignum+6,[outdir,'Energy_density'],formats);
        %%{
        for icrv=1:Ncrv
            savenameHP = strcat(outdir,'WL_HP_mat_',sprintf('%d',icrv));
            save_fig(fignum+10+icrv,savenameHP,formats);
            savenameHE = strcat(outdir,'WL_HE_mat_',sprintf('%d',icrv));
            save_fig(fignum+20+icrv,savenameHE,formats);
        end
        %}
end


%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_mat = var_extract(var_name,outdir,Npnt,Ncrv)

% Create a matrix by extracting and arranging the values of the variable
% 'var_name' from the different files inside the Multi_run folder.
var_mat   = zeros(Npnt,Ncrv);
for icrv=1:Ncrv
    for ipnt=1:Npnt
        filename = sprintf([outdir,'Crv_%d_Pnt_%d.mat'],icrv,ipnt);
        contents = whos('-file',filename);
        if any(strcmp({contents.name},var_name))
            S = load(filename,var_name);
            var_mat(ipnt,icrv)  = S.(var_name);
        else
        end
    end
end

end

function var_mat = mat_extract(mat_name,outdir,Npnt,Ncrv)

% Create a 'matrix of matrices' by extracting and arranging the values of
% the matrix 'mat_name' from the different files inside the Multi_run
% folder.

var_mat   = cell(Npnt,Ncrv);
for icrv=1:Ncrv
    for ipnt=1:Npnt
        filename = sprintf([outdir,'Crv_%d_Pnt_%d.mat'],icrv,ipnt);
        S = load(filename,mat_name);
        var_mat{ipnt,icrv}  = S.(mat_name);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%