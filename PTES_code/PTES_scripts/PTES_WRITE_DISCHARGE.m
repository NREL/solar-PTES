% Write discharge cycles
for iL=1:Load.num
    if any(strcmp(Load.type(iL),{'dis','disCO2'}))
        stages_dis = gas.Nstg(iL); 
        for ind = 1:stages_dis
            write_file(gas,[iL,ind],plotFile,num);
        end
        break
    end
    if strcmp(Load.type(iL),'ran')
        stages_dis = steam.Nstg(iL);
        for ind = 1:stages_dis
            write_file(steam,[iL,ind],plotFile,num);
        end
        break
    end
end