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
        
        %W_net_chg    = 0;
        %W_net_dis    = 0;
        
        
        save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','rhoE',...
            'W_in_chg','W_out_dis',...
            'E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
            'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected',...
            'WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
            'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',...
            'Ttop','Tbot','First_law_error','Second_law_error');
    case 1
        chi = chi_hot;
        EFF=0;
    case 2
        chi = chi_tot;
        COP=0;
    case 3
        
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
        
        W_net_chg    = 0;
        W_net_dis    = 0;
        
        
        save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','rhoE',...
            'W_in_chg','W_out_dis',...
            'E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
            'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected',...
            'WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
            'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',...
            'Ttop','Tbot','First_law_error','Second_law_error',...
            'HEeff','HEeffRC','HEeffNC');
end

