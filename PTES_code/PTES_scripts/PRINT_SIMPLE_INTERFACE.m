fprintf(1,'\n\n');
fprintf(1,'-----------------------------\n');
fprintf(1,'SIMPLE INTERFACE RESULTS\n');
fprintf(1,'Input these numbers into SAM.\n');
fprintf(1,'-----------------------------\n');

fprintf(1,'\nSystem Design Page\n');
fprintf(1,'Heat pump COP:                      %8.2f \n',QH_chg/W_in_chg);
fprintf(1,'Cycle thermodynamic power:          %8.2f MWe\n',W_out_dis/t_dis/1e6);
fprintf(1,'Cycle thermodynamic efficiency:     %8.2f %%\n\n',W_out_dis/QH_dis*100);

fprintf(1,'Hot storage hot temperature:        %8.2f C\n',fluidH.state(1,fluidH(1).Nstg(1)+1).T-273.15);
fprintf(1,'Hot storage cold temperature:       %8.2f C\n',fluidH.state(1,1).T-273.15);
fprintf(1,'Cold storage hot temperature:       %8.2f C\n',fluidC.state(1,1).T-273.15);
fprintf(1,'Cold storage cold temperature:      %8.2f C\n',fluidC.state(1,fluidC(1).Nstg(1)+1).T-273.15);

fprintf(1,'\nHeat Pump Page\n');

fprintf(1,'Parastics (non-pumping):            %8.3f \n',(E_in_chg-W_in_chg+W_fan_chg)/E_in_chg);
fprintf(1,'Pumping power rate through hot HX:  %8.2f kW/kg/s\n',-CPMP(1).W(1)/t_chg/1e3/fluidH.state(1,1).mdot);
fprintf(1,'Pumping power rate through cold HX: %8.2f kW/kg/s\n',-CPMP(2).W(1)/t_chg/1e3/fluidC.state(1,1).mdot);

fprintf(1,'\nPower Cycle Page\n');

fprintf(1,'Parastics (non-pumping):            %8.3f \n',(W_out_dis-E_out_dis+W_fan_dis + heater_in)/W_out_dis);
fprintf(1,'Pumping power rate through hot HX:  %8.2f kW/kg/s\n',-DPMP(2).W(2)/t_chg/1e3/fluidH.state(2,1).mdot);
fprintf(1,'Pumping power rate through cold HX: %8.2f kW/kg/s\n',-DPMP(1).W(2)/t_chg/1e3/fluidC.state(2,1).mdot);

fprintf(1,'\nSystem Costs Page\n');
HPcost = mean(compMAT(:,1)+compMAT(:,3)+compMAT(:,5)+compMAT(:,6)+compMAT(:,7)+compMAT(:,8)+compMAT(:,9)+compMAT(:,15));
fprintf(1,'Heat pump cost:  %8.2f $/kWt\n',HPcost/abs(QH_chg/t_chg/1e3));
fprintf(1,'N.B. Heat pump cost includes all heat exchangers, heat rejection equip., motor/generator, pumps, fans, as well as charging turbomachinery.\n');
HEcost = mean(compMAT(:,2)+compMAT(:,4));
fprintf(1,'Power cycle cost:  %8.2f $/kWe\n',HEcost/abs(W_out_dis/t_dis/1e3));
fprintf(1,'N.B. Power cycle cost just includes the turbomachinery.\n');

HTEScost = mean(compMAT(:,11)+compMAT(:,13));
fprintf(1,'Hot thermal energy storage cost:  %8.2f $/kWht\n',HTEScost/abs(QH_chg/3600e3));

CTEScost = mean(compMAT(:,12)+compMAT(:,14));
fprintf(1,'Cold thermal energy storage cost:  %8.2f $/kWht\n',CTEScost/abs(QC_chg/3600e3));
fprintf(1,'-----------------------------\n\n');




%%%%
% Write the same output into json format
% Names correspond to variables in SAM
JSONFILE_name= 'Outputs/PTES_output_for_SAM.json'; 
fid=fopen(JSONFILE_name,'w');

s.tshours = t_dis/3600 ;
s.T_HT_cold_htf_des       = fluidH.state(1,1).T-273.15;
s.W_dot_pc_thermo_des      = W_out_dis/t_dis/1e6;
s.cop_hp_thermo_des                  = QH_chg/W_in_chg;
s.eta_pc_thermo_des = W_out_dis/QH_dis*100;
s.T_HT_hot_htf_des        = fluidH.state(1,fluidH(1).Nstg(1)+1).T-273.15;
s.T_CT_hot_htf_des       = fluidC.state(1,1).T-273.15;
s.T_CT_cold_htf_des      = fluidC.state(1,fluidC(1).Nstg(1)+1).T-273.15;


s.f_hp_parasitic_des            = (E_in_chg-W_in_chg+W_fan_chg)/E_in_chg;
s.heat_pump_HT_HTF_pump_coef  = -CPMP(1).W(1)/t_chg/1e3/fluidH.state(1,1).mdot;
s.heat_pump_CT_HTF_pump_coef = -CPMP(2).W(1)/t_chg/1e3/fluidC.state(1,1).mdot;


s.f_pc_parasitic_des            = (W_out_dis-E_out_dis+W_fan_dis + heater_in)/W_out_dis;
s.pb_pump_coef  = -DPMP(2).W(2)/t_chg/1e3/fluidH.state(2,1).mdot;
s.CT_pb_pump_coef = -DPMP(1).W(2)/t_chg/1e3/fluidC.state(2,1).mdot;

s.heat_pump_spec_cost  = HPcost/abs(QH_chg/t_chg/1e3);
s.cycle_spec_cost = HEcost/abs(W_out_dis/t_dis/1e3);

s.tes_spec_cost = HTEScost/abs(QH_chg/3600e3);
s.CT_tes_spec_cost = CTEScost/abs(QC_chg/3600e3);

s.heater_multiple = 1;


encodedJSON = jsonencode(s); 
fprintf(fid, encodedJSON);
fclose('all');





