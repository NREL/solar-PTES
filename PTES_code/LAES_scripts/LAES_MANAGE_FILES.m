%Open the output files and print the headers
if multi_run
    ID0 = fopen('./Outputs/Multi_run_var.txt','w');
    
    ID1 = fopen('./Outputs/Multi_run.txt','w');
    fprintf(ID1,'%%');
    fprintf(ID1,'%15s','PR [-]','eta [-]','eff [-]','ploss [-]','X_PTES [%]','COP [-]','EFF [%]','Max T [K]','Min T [K]','WR_dis [-]','rhoE [kWh/m3]','rhoP [MW/(m3/s)]','WL_comp [%]','WL_exp [%]','WL_hexs [%]','WL_reject [%]','WL_mix_liq [%]','WL_tanks [%]');
    fprintf(ID1,'\n');
    
    LAES_SET_MULTI_VAR;
else
    ID0 = NaN;
    ID1 = NaN;
end

plotFile1 = fopen('./Outputs/Plot_file1.txt','w+');
fprintf(plotFile1,'%%Things to plot\n%%');
fprintf(plotFile1,'%12s %12s %12s %12s %12s\n','T(K)','T(°C)','p(bar)','s(J/kg/K)','Stage');

plotFile2 = fopen('./Outputs/Plot_file2.txt','w+');
fprintf(plotFile2,'%%Things to plot\n');
fprintf(plotFile2,'%12s %12s %12s %12s %12s\n','T(K)','T(°C)','p(bar)','s(J/kg/K)','Stage');

plotFile3 = fopen('./Outputs/Plot_file3.txt','w+');
fprintf(plotFile3,'%%Things to plot\n');
fprintf(plotFile3,'%12s %12s %12s %12s %12s\n','T(K)','T(°C)','p(bar)','s(J/kg/K)','Stage');