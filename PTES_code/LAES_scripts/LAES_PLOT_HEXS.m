% PLOT HEX T-Q's (charge)
% Hot HEX
ip = 2;
plot_hex(gas,[1,ip],fluidH(1),[1,1],100,5);
% Cold HEX
for i=1:stages_ch
    if strcmp(gas.stage(1,i).type,'exp')
        ip = i + 1;
        plot_hex(gas,[1,ip],fluidC(1),[1,1],100,7);
        break
    end
end