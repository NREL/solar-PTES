switch Load.mode
    case {0,3}
        chi = chi_PTES;
        COP=0;
        EFF=0;
    case 1
        chi = chi_hot;
        EFF=0;
    case 2
        chi = chi_tot;
        COP=0;
end

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',icrv,ipnt);
save(filename,'PRch','eta','eff','ploss','chi','COP','EFF','rhoE',...
    'WL_comp','WL_exp','WL_hexs','WL_reject','WL_mix_liq','WL_mix_gas',...
    'WL_tanks'...
    ,'t_chg','t_dis','t_disRC','t_disNC',...
    'HEeff','HEeffRC','HEeffNC');