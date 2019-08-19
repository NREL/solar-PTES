% Write charge cycle
stages1_ch = find(~strcmp({gas1.stage(1,:).type},'0'),1,'last');
for ind = 1:stages1_ch
    write_file(gas1,[1,ind],plotFile1,num);
end
stages2_ch = find(~strcmp({gas2.stage(1,:).type},'0'),1,'last');
for ind = 1:stages2_ch
    write_file(gas2,[1,ind],plotFile2,num);
end
stages3_ch = find(~strcmp({gas3.stage(1,:).type},'0'),1,'last');
for ind = 1:stages3_ch
    write_file(gas3,[1,ind],plotFile3,num);
end