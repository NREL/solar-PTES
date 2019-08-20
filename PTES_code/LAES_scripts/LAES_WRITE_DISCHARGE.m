% Write discharge cycles
stages1_dis = find(~strcmp({gas1.stage(2,:).type},'0'),1,'last');
for ind = 1:stages1_dis
    write_file(gas1,[2,ind],plotFile1,num);
end