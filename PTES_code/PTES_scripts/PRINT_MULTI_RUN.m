filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',icrv,ipnt);

switch Load.mode
    case 0
        chi = chi_PTES;
        
        cap_costM    = Cdata.cap_costM ;
        cap_costSD   = Cdata.cap_costSD ;
        cap_cost_pow = Cdata.cap_cost_pow ;
        cap_cost_en  = Cdata.cap_cost_en ;
        lcosM        = Cdata.lcosM ;
        lcosSD       = Cdata.lcosSD ;
        mdot         = Load.mdot(1) ;
        Ttop         = gas.state(1,2).T ;
        Tbot         = gas.state(1,6).T ;
        
        
        save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','rhoE',...
            'W_in_chg','W_net_chg','W_out_dis','W_net_dis',...
            'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected',...
            'WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
            'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',...
            'Ttop','Tbot');
    case 1
        chi = chi_hot;
        EFF=0;
    case 2
        chi = chi_tot;
        COP=0;
    case 3
        chi = chi_PTES;
        COP=0;
        EFF=0;
        
        save(filename,'PRch','eta','eff','ploss','chi','COP','EFF','rhoE',...
            'WL_comp','WL_exp','WL_hexs','WL_reject','WL_mix_liq','WL_mix_gas',...
            'WL_tanks'...
            ,'t_chg','t_dis','t_disRC','t_disNC',...
            'HEeff','HEeffRC','HEeffNC');
end

