fid=fopen(JSONFILE_name,'w');

i_chg = Design_Load.ind(any(Design_Load.type == {'chg'},2));
i_dis = Design_Load.ind(any(Design_Load.type == {'dis'},2));

HEcost = mean(compMAT(:,2)+compMAT(:,4));
HPcost = mean(compMAT(:,1)+compMAT(:,3)+compMAT(:,5)+compMAT(:,6)+compMAT(:,7)+compMAT(:,8)+compMAT(:,9)+compMAT(:,15));

HTEScost = mean(compMAT(:,11)+compMAT(:,13));
CTEScost = mean(compMAT(:,12)+compMAT(:,14));

s.tshours             = t_dis/3600 ;
s.T_HT_cold_htf_des   = fluidH.state(i_chg(1),1).T-273.15;
s.W_dot_pc_thermo_des = W_out_dis/t_dis/1e6;
s.cop_hp_thermo_des   = QH_chg/W_in_chg;
s.eta_pc_thermo_des   = W_out_dis/QH_dis*100;
s.T_HT_hot_htf_des    = fluidH.state(i_chg(1),fluidH(1).Nstg(1)+1).T-273.15;
s.T_CT_hot_htf_des    = fluidC.state(i_chg(1),1).T-273.15;
s.T_CT_cold_htf_des   = fluidC.state(i_chg(1),fluidC(1).Nstg(1)+1).T-273.15;

s.f_hp_parasitic_des         = (E_in_chg-W_in_chg+W_fan_chg)/E_in_chg;
s.heat_pump_HT_HTF_pump_coef = -CPMP(1).W(i_chg(1))/t_chg/1e3/fluidH.state(i_chg(1),1).mdot;
s.heat_pump_CT_HTF_pump_coef = -CPMP(2).W(i_chg(1))/t_chg/1e3/fluidC.state(i_chg(1),1).mdot;

s.f_pc_parasitic_des = (W_out_dis-E_out_dis+W_fan_dis + heater_in)/W_out_dis;
s.pb_pump_coef       = -DPMP(2).W(i_dis(1))/t_chg/1e3/fluidH.state(i_dis(1),1).mdot;
s.CT_pb_pump_coef    = -DPMP(1).W(i_dis(1))/t_chg/1e3/fluidC.state(i_dis(1),1).mdot;

s.heat_pump_spec_cost = HPcost/abs(QH_chg/t_chg/1e3);
s.cycle_spec_cost     = HEcost/abs(W_out_dis/t_dis/1e3);

s.tes_spec_cost    = HTEScost/abs(QH_chg/3600e3);
s.CT_tes_spec_cost = CTEScost/abs(QC_chg/3600e3);

s.heater_multiple = HP_mult;

encodedJSON = jsonencode(s); 
fprintf(fid, encodedJSON);
fclose('all');

