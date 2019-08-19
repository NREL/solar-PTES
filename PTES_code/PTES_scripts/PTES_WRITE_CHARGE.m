% Write charge cycles

stages_ch = find(~strcmp({gas.stage(1,:).type},'0'),1,'last')-1;
for ind = 1:stages_ch
    write_file(gas,[1,ind],plotFile,num);
end