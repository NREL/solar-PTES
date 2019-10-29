% Write discharge cycles
for iL=1:Load.num
    if any(strcmp(Load.type(iL),{'dis','disCO2'}))
        stages_dis = find(~strcmp({gas.stage(iL,:).type},'0'),1,'last')-1;
        for ind = 1:stages_dis
            write_file(gas,[iL,ind],plotFile,num);
        end
        break
    end
end