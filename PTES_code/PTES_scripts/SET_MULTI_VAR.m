fprintf(ID0,'%s\n',Vpnt);
fprintf(ID0,'%d\n',Npnt);
fprintf(ID0,'%f\n',Apnt);
fprintf(ID0,'\n');

fprintf(ID0,'%s\n',Vcrv);
fprintf(ID0,'%d\n',Ncrv);
fprintf(ID0,'%f\n',Acrv);
fprintf(ID0,'\n');


% Make sure that enough liquid streams are declared
if strcmp(Vpnt,'Ne_ch')
    Ne_ch = max(Apnt);
end
if strcmp(Vcrv,'Ne_ch')
    Ne_ch = max(Acrv);
end