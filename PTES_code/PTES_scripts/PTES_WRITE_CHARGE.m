% Write charge cycles
for iL=1:Load.num
    if strcmp(Load.type(iL),'chg')
        stages_ch = find(~strcmp({gas.stage(iL,:).type},'0'),1,'last')-1;
        for ind = 1:stages_ch
            write_file(gas,[iL,ind],plotFile,num);
        end
        break
    end
end