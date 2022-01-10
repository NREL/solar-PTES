% PLOT HEX T-Q's

if Load.mode == 3
    %%{
    plot_hex(HX(ihx_hot(1)),1,31,'C',true) % Molten salt
    plot_hex(HX(ihx_reg(1)),1,32,'C',true) % Regenerator
    if ~any(Load.options.useCold)
        plot_hex(HX(ihx_hin(1)),1,33,'C',true) % Heat intake unit
    else
        plot_hex(HX(ihx_cld(1)),1,34,'C',true) % Cold hex
        plot_hex(HX(ihx_htf(1)),1,35,'C',true) % HTF-water
    end
    if ~Load.options.superRank
        plot_hex(HX(ihx_JB+1),3,36,'C',true) % Reheater
        if Load.options.useCold(3)
            plot_hex(HX(ihx_JB+2),3,37,'C',true) % Condenser
        end
        plot_hex(HX(ihx_JB+3),4,38,'C',true) % Air-cooled condenser
        plot_hex(HX(ihx_JB+4),3,39,'C',true) % Boiler
    else
        plot_hex(HX(ihx_JB+1),3,36,'C',true) % First Reheater
        plot_hex(HX(ihx_JB+2),3,36,'C',true) % Second Reheater
        if Load.options.useCold(3)
            plot_hex(HX(ihx_JB+3),3,37,'C',true) % Condenser
        end
        plot_hex(HX(ihx_JB+4),4,38,'C',true) % Air-cooled condenser
        plot_hex(HX(ihx_JB+5),3,39,'C',true) % Boiler
    end
    %%}
    %{
    formats = {'epsc','fig'};
    save_fig(31,'./Results/TQ_JB_hot'  ,formats)
    save_fig(32,'./Results/TQ_JB_regen',formats)
    if ~any(Load.options.useCold)
        save_fig(33,'./Results/TQ_JB_rej'  ,formats)
    else
        save_fig(34,'./Results/TQ_JB_cold' ,formats)
        save_fig(35,'./Results/TQ_JB_htf'  ,formats)
        save_fig(37,'./Results/TQ_RK_cold' ,formats)
    end
    save_fig(36,'./Results/TQ_RK_hot'  ,formats)
    save_fig(38,'./Results/TQ_RK_rej'  ,formats)
    save_fig(39,'./Results/TQ_RK_boil' ,formats)
    %}
end

if Load.mode == 20
    plot_hex(HX(4),1,41,'K',true) % High Coupler
    plot_hex(HX(5),1,42,'K',true) % Cryo Coupler
    plot_hex(HX(8),1,43,'K',true) % Regenerator
end