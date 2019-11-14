% PLOT HEX T-Q's
plot_hex_TQA(HX1,20,'C')

if Load.mode == 3
    plot_hex_TQA(HX_REHEAT,21,'C')
    plot_hex_TQA(HX_BOILER,22,'C')
end