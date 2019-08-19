% Write discharge cycles

stages_dis = find(~strcmp({gas.stage(2,:).type},'0'),1,'last')-1;
for ind = 1:stages_dis
    write_file(gas,[2,ind],plotFile,num);
end