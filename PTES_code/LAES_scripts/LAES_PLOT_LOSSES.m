figure(8)
switch mode
    case 0
        names = {'$$\mathrm{PTES_{ch}}$$','$$\mathrm{PTES_{dis}}$$'};
        b = bar(WL_matrix./W_in_ch*100,'stacked');
        set(gca, 'XTick', 1:2, 'XTickLabel', names, 'TickLabelInterpreter', 'latex')
    case 1
        names = {'Heat pump'};
        b = bar(WL_matrix./W_in_ch*100,'stacked');
        set(gca, 'XTick', 1, 'XTickLabel', names, 'TickLabelInterpreter', 'latex')
    case 2
        names = {'Heat engine'};
        b = bar(WL_matrix./Exergy_from_tanks*100,'stacked');
        set(gca, 'XTick', 2, 'XTickLabel', names, 'TickLabelInterpreter', 'latex')
        %xlim([0 3.0])
end
b(2).FaceColor = c_pale_blue;
b(3).FaceColor = c_pale_green;
b(4).FaceColor = 'yellow';
b(5).FaceColor = c_dark_orange;
b(6).FaceColor = c_grey;
ylabel(strcat('Lost Work [$$ \% $$]'))
legend({'Compressors','Expanders','Heat exchangers','Heat rejection','Mixing (liquid)','Tanks'},'Location','Best')


% Do not show liquid_mixing loss bar if it is not required
if all([Nc_ch,Ne_ch]==1)
    delete(b(5))
end