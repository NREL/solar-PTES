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

fprintf(ID1,'%15.5g', PRch , eta, eff, ploss, chi,COP,EFF,...
    HT.B(2).T, CT.B(2).T, rhoE,...
    WL_comp, WL_exp, WL_hexs, WL_reject, WL_mix_liq, WL_mix_gas, WL_tanks);
fprintf(ID1,'\n');