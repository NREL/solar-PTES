% PLOT HEX T-Q's

if Load.mode == 3
    %%{
    plot_hex(HX(1),1,31,'C',true)
    plot_hex(HX(2),1,32,'C',true)
    plot_hex(HX(3),1,33,'C',true)
    plot_hex(HX(4),1,34,'C',true)
    plot_hex(HX(5),1,35,'C',true)
    plot_hex(HX(6),3,36,'C',true)
    plot_hex(HX(7),3,37,'C',true)
    plot_hex(HX(8),4,38,'C',true)
    plot_hex(HX(9),3,39,'C',true)
    %%}
    %%{
    save_fig(31,'./Results/TQ_JB_hot'  ,'epsc')
    save_fig(32,'./Results/TQ_JB_regen','epsc')
    save_fig(33,'./Results/TQ_JB_rej'  ,'epsc')
    save_fig(34,'./Results/TQ_JB_cold' ,'epsc')
    save_fig(35,'./Results/TQ_JB_htf'  ,'epsc')
    save_fig(36,'./Results/TQ_RK_hot'  ,'epsc')
    save_fig(37,'./Results/TQ_RK_cold' ,'epsc')
    save_fig(38,'./Results/TQ_RK_rej'  ,'epsc')
    save_fig(39,'./Results/TQ_RK_boil' ,'epsc')
    %%}
end