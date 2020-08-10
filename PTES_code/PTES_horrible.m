
% #1
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [20.*fac;20.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,1);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

%}

% #2
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [25.*fac;25.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,2);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #3
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [30.*fac;30.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,3);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #4
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [35.*fac;35.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,4);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #5
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [40.*fac;40.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,5);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #6
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [45.*fac;45.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,6);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #7
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [50.*fac;50.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,7);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #8
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [55.*fac;55.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,8);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #9
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [60.*fac;60.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,9);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #10
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [65.*fac;65.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,10);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #11
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [70.*fac;70.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,11);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #12
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [75.*fac;75.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,12);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #13
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [80.*fac;80.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,13);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #14
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [85.*fac;85.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,14);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');



% #15
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [90.*fac;90.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,15);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #16
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [95.*fac;95.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,16);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #17
clear
clear global

T0 = 25 + 273.15 ; fac = 10 ;
mdotIN = [100.*fac;100.*fac]; T0IN =[T0-30;T0-30] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,17);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

