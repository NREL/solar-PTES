filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',icrv,ipnt);

mdot         = Load.mdot(1) ;
if exist('Cdata','var')
    cap_costM    = Cdata.cap_costM ;
    cap_costSD   = Cdata.cap_costSD ;
    cap_cost_pow = Cdata.cap_cost_pow ;
    cap_cost_en  = Cdata.cap_cost_en ;
    lcosM        = Cdata.lcosM ;
    lcosSD       = Cdata.lcosSD ;
else
    cap_costM    = 0 ;
    cap_costSD   = 0 ;
    cap_cost_pow = 0 ;
    cap_cost_en  = 0 ;
    lcosM        = 0 ;
    lcosSD       = 0 ;
end
Ttop         = gas.state(1,2).T ;
Tbot         = gas.state(1,6).T ;
Load_mode = Load.mode;
switch Load_mode
    case 0
        chi = chi_PTES;
                
        save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','chi_PTES_true',...
            'rhoE','W_in_chg','W_out_dis',...
            'E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
            'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected',...
            'WL_matrix','t_chg','t_dis','mdot',...
            'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',...
            'THPmin','Ttop','Tbot','First_law_error','Second_law_error','Load_mode');
    case 1
        chi = chi_hot;
        EFF=0;
    case 2
        chi = chi_tot;
        COP=0;
    case 3
        chi = chi_PTES;        
        
        save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','chi_PTES_true',...
            'rhoE','W_in_chg','W_out_dis',...
            'E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
            'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected',...
            'WL_matrix','t_chg','t_dis','t_disRC','t_disNC','mdot',...
            'COP','HEeff','HEeffRC','HEeffNC',...
            'HPexergy_eff','HEexergy_eff','HEexergy_effRC','HEexergy_effNC',...
            'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',...
            'THPmin','Ttop','Tbot','First_law_error','Second_law_error','Load_mode');
        
    case {22,24}
        chi = chi_PTES;
        THmax = HT.B(2).T;
        
        save(filename,'eta','eff','ploss','pmax_LA','chi_PTES','chi_PTES_para','chi_PTES_true',...
            'rhoE','W_in_chg','W_out_dis',...
            'E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
            'QH_chg','QH_dis','Heat_rejected',...
            'WL_matrix','t_chg','t_dis','mdot',...
            'First_law_error','Second_law_error','Load_mode','Qstr_W_ratio','QT_W_ratio','TLA','Te3','THmax');
       
end

