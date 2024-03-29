figure(8)
switch Load.mode
    case {0,3,4}
        names = {'$$\mathrm{PTES_{ch}}$$','$$\mathrm{PTES_{dis}}$$'};
        b = bar(WL_matrix./W_in_chg*100,'stacked');
        set(gca, 'XTick', 1:2, 'XTickLabel', names, 'TickLabelInterpreter', 'latex')
    case 1
        names = {'Heat pump'};
        b = bar(WL_matrix./W_in_chg*100,'stacked');
        set(gca, 'XTick', 1, 'XTickLabel', names, 'TickLabelInterpreter', 'latex')
    case {2,5}
        names = {'Heat engine'};
        b = bar(WL_matrix./Exergy_from_tanks*100,'stacked');
        set(gca, 'XTick', 2, 'XTickLabel', names, 'TickLabelInterpreter', 'latex')
        %xlim([0 3.0])
    case 6
        names = {'$$\mathrm{PTES_{ch}}$$','$$\mathrm{PTES_{dis}}$$'};
        b = bar(WL_matrix./(W_in_chg + EX_sol)*100,'stacked');
        set(gca, 'XTick', 1:2, 'XTickLabel', names, 'TickLabelInterpreter', 'latex')
end
ylabel(strcat('Lost Work [$$ \% $$]'))
b(1).FaceColor = c_dark_blue;
b(2).FaceColor = c_pale_blue;
b(3).FaceColor = c_pale_green;
b(4).FaceColor = c_yellow;
b(5).FaceColor = c_pale_orange;
b(6).FaceColor = c_dark_orange;
b(7).FaceColor = c_grey;
legend({'Compressors','Expanders','Heat exchangers','Heat in/out env.','Mixing (liquid)','Mixing (gas)','Tanks'},'Location','Best')

% Do not show liquid_mixing loss and tank_loss bars if they are not required
% if WL_mix_liq==0
%     delete(b(5))
% end
% if WL_tanks==0
%     delete(b(6))
% end


% % Assuming cold reject in heat pump mode
% if mode == 1
%     figure(9)
% 
%     names = {'Heat pump'};
%     b = bar(WL_matrix./W_in_ch*100,'stacked');
%     set(gca, 'XTick', 1, 'XTickLabel', names, 'TickLabelInterpreter', 'latex')
%     b(2).FaceColor = c_pale_blue;
%     b(3).FaceColor = c_pale_green;
%     b(4).FaceColor = 'yellow';
%     b(5).FaceColor = c_dark_orange;
%     b(6).FaceColor = c_grey;
%     ylabel(strcat('Lost Work [$$ \% $$]'))
%     legend({'Compressors','Expanders','Heat exchangers','Cold rejection','Mixing (liquid)','Tanks'},'Location','Best')
% 
% end