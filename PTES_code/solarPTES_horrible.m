
% #1
clear
clear global

T0 = 40 + 273.15 ; fac = 25 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES
if ~exist('./Outputs/Multi_run','dir')
    mkdir('./Outputs/Multi_run')
end
filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,1);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');



% #2
clear
clear global

T0 = 40 + 273.15 ; fac = 30 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,2);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #3
clear
clear global

T0 = 40 + 273.15 ; fac = 35 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,3);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #4
clear
clear global

T0 = 40 + 273.15 ; fac = 40 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,4);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #5
clear
clear global

T0 = 40 + 273.15 ; fac = 45 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,5);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #6
clear
clear global

T0 = 40 + 273.15 ; fac = 50 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,6);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #7
clear
clear global

T0 = 40 + 273.15 ; fac = 55 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,7);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #8
clear
clear global

T0 = 40 + 273.15 ; fac = 60 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,8);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #9
clear
clear global

T0 = 40 + 273.15 ; fac = 63 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,9);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #10
clear
clear global

T0 = 40 + 273.15 ; fac = 70 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,10);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #11
clear
clear global

T0 = 40 + 273.15 ; fac = 75 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,11);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #12
clear
clear global

T0 = 40 + 273.15 ; fac = 80 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,12);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');

% #13
clear
clear global

T0 = 40 + 273.15 ; fac = 85 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,13);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #14
clear
clear global

T0 = 40 + 273.15 ; fac = 90 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,14);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');



% #15
clear
clear global

T0 = 40 + 273.15 ; fac = 95 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,15);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


% #16
clear
clear global

T0 = 40 + 273.15 ; fac = 100 ;
mdotIN = [10*fac;0;1*fac;1*fac]; T0IN =[T0-0;T0-0;T0-0;T0-0] ;
PTES

filename = sprintf('./Outputs/Multi_run/Crv_%d_Pnt_%d.mat',1,16);

cap_costM = Cdata.cap_costM ; cap_costSD = Cdata.cap_costSD ; cap_cost_pow = Cdata.cap_cost_pow ;
cap_cost_en = Cdata.cap_cost_en ; lcosM = Cdata.lcosM ; lcosSD = Cdata.lcosSD ; mdot = Load.mdot(1) ;
Ttop = gas.state(1,2).T ; Tbot = gas.state(1,6).T ;
COP_para  = QH_chg/E_net_chg;
HE_para = E_net_dis/QH_dis ;

save(filename,'PRch','PRr','eta','eff','ploss','pmax','chi_PTES','chi_PTES_para','COP_para','HE_para','rhoE',...
    'HEeff','HEeffRC','HEeffNC',...
    'W_in_chg','W_out_dis','E_in_chg','E_net_chg','E_out_dis','E_net_dis',...
    'QH_chg','QH_dis','QC_chg','QC_dis','Heat_rejected','WL_PTES_chg','WL_PTES_dis','t_chg','t_dis','mdot',...
    'cap_costM','cap_costSD','cap_cost_pow','cap_cost_en','lcosM','lcosSD',    'Ttop','Tbot','First_law_error','Second_law_error');


load gong.mat
sound(y)