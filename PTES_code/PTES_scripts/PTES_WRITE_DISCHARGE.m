% Write discharge cycles
for iL=1:Load.num
    if strcmp(Load.type(iL),'dis')
        stages_dis = gas.Nstg(iL); 
        for ind = 1:stages_dis
            write_file(gas,[iL,ind],plotFile,num);
        end
        break
    end
    if strcmp(Load.type(iL),'ran')
        stages_dis = steam(1).Nstg(iL);
        for ind = 1:stages_dis
            write_file(steam(1),[iL,ind],plotFile,num);
        end
        break
    end
end