function [ ID0, ID1, ID2] = manage_files(multi_run)
%Open the output files and print the headers

if multi_run
    ID0 = fopen('./Outputs/Multi_run_var.txt','w');
    
    ID1 = fopen('./Outputs/Multi_run.txt','w');
    fprintf(ID1,'%%');
    fprintf(ID1,'%15s','PR [-]','eta [-]','eff [-]','ploss [-]','X_PTES [%]','COP [-]','Max T [K]','Min T [K]','WR_dis [-]','rhoE [kWh/m3]','rhoP [MW/(m3/s)]','WL_comp [%]','WL_exp [%]','WL_hexs [%]','WL_reject [%]','WL_mix_liq [%]','WL_tanks [%]');
    fprintf(ID1,'\n');
    
    PTES_SET_MULTI_VAR;
else
    ID0 = NaN;
    ID1 = NaN;
end

ID2 = fopen('./Outputs/Plot_file.txt','w+');
fprintf(ID2,'%%Things to plot\n');
fprintf(ID2,'%%T(K)\t   T(°C)\t p(Pa)\t s(J/K)\n');

end

