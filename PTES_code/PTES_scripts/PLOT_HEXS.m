% PLOT HEX T-Q's
%plot_hex(HX,iL,20,'C')

if Load.mode == 3
    %plot_hex(HX_REHEAT,iL,21,'C')
    %plot_hex(HX_BOILER,iL,22,'C')
    %plot_hex(HX_CONDEN,iL,23,'C')
    %plot_hex(HX_ACC,iL,23,'C')
    plot_hex(HX(8),3,30,'C')
    figure(30)
    ylim([0 45])
end