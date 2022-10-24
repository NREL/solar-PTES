filename      = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',icrv,ipnt);
JSONFILE_name = sprintf('./Outputs/Multi_run_json/PTES_output_for_SAM_%s_%g_%s_%g.json',Vcrv,Acrv(icrv),Vpnt,Apnt(ipnt));

switch Load.mode
    case {0,4}

        PRINT_JSON

        if ~Loffdesign
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


            save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','chi_PTES_true',...
                'rhoE','W_in_chg','W_out_dis',...
                'E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
                'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected',...
                'WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
                'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',...
                'Ttop','Tbot','First_law_error','Second_law_error');

        else
            filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',icrv,ipnt);

            CCMPeta = CCMP.eta(i_chg) ;
            CCMPpr  = CCMP.pr(i_chg) ;
            CEXPeta = CEXP.eta(i_chg) ;
            CEXPpr  = CEXP.pr(i_chg) ;

            DCMPeta = DCMP.eta(i_dis) ;
            DCMPpr  = DCMP.pr(i_dis) ;
            DEXPeta = DEXP.eta(i_dis) ;
            DEXPpr  = DEXP.pr(i_dis) ;

            lastci = fluidH.Nstg(i_chg)+1;
            lastdi = fluidH.Nstg(i_dis)+1;

            flHmdotC = fluidH.state(i_chg,1).mdot;
            flHT1C   = fluidH.state(i_chg,1).T;
            flHT2C   = fluidH.state(i_chg,lastci).T;

            flHmdotD = fluidH.state(i_dis,1).mdot;
            flHT1D   = fluidH.state(i_dis,1).T;
            flHT2D   = fluidH.state(i_dis,lastdi).T;

            lastci = fluidC.Nstg(i_chg)+1;
            lastdi = fluidC.Nstg(i_dis)+1;

            flCmdotC = fluidC.state(i_chg,1).mdot;
            flCT1C   = fluidC.state(i_chg,1).T;
            flCT2C   = fluidC.state(i_chg,lastci).T;

            flCmdotD = fluidC.state(i_dis,1).mdot;
            flCT1D   = fluidC.state(i_dis,1).T;
            flCT2D   = fluidC.state(i_dis,lastdi).T;

            hSOC_final = HT.A(i_dis+1).M / HT.A(i_chg).M;
            cSOC_final = CT.A(i_dis+1).M / CT.A(i_chg).M;

            save(filename,'chi_PTES','chi_PTES_para','chi_PTES_true',...
                'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
                'QH_chg','QH_dis',...
                'CCMPeta','CCMPpr','CEXPeta','CEXPpr','DCMPeta','DCMPpr','DEXPeta','DEXPpr',...
                'flHmdotC','flHT1C','flHT2C','flHmdotD','flHT1D','flHT2D',...
                'flCmdotC','flCT1C','flCT2C','flCmdotD','flCT1D','flCT2D',...
                't_chg','t_dis','hSOC_final','cSOC_final',...
                'First_law_error','Second_law_error');
        end
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
        
        
        save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','chi_PTES_true',...
            'rhoE','W_in_chg','W_out_dis',...
            'E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
            'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected',...
            'WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','t_disRC','t_disNC','mdot',...
            'HEeff','HEeffRC','HEeffNC',...
            'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',...
            'Ttop','Tbot','First_law_error','Second_law_error');
       
end

