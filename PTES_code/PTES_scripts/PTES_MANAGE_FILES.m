%Open the output files and print the headers
if multi_run
    ID0 = fopen('./Outputs/Multi_run_var.txt','w');
    
    ID1 = fopen('./Outputs/Multi_run.txt','w');
    fprintf(ID1,'%%');
    fprintf(ID1,'%15s','PRch [-]','eta [-]','eff [-]','ploss [-]','X_PTES [%]','COP [-]','EFF [%]','Max T [K]','Min T [K]','WR_dis [-]','rhoE [kWh/m3]','rhoP [MW/(m3/s)]','WL_comp [%]','WL_exp [%]','WL_hexs [%]','WL_reject [%]','WL_mix_liq [%]','WL_tanks [%]');
    fprintf(ID1,'\n');
    
    PTES_SET_MULTI_VAR;
else
    ID0 = NaN;
    ID1 = NaN;
end

plotFile = fopen('./Outputs/Plot_file.txt','w+');
fprintf(plotFile,'%%Things to plot\n');
fprintf(plotFile,'%%T(K)\t   T(°C)\t p(bar)\t s(kJ/kg/K)\t Stage\n');