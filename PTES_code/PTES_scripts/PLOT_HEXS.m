% PLOT HEX T-Q's
plot_hex(HX,20,'C')

if Load.mode == 3
    plot_hex(HX_REHEAT,21,'C')
    plot_hex(HX_BOILER,22,'C')
    plot_hex(HX_CONDEN,23,'C')
    plot_hex(HX_ACC,23,'C')
end