% PLOT HEX T-Q's (charge)
% Hot HEX
ip = 2;
plot_hex(gas,[1,ip],fluidH(1),[1,1],100,14);
% Cold HEX
for i=1:stages_ch
    if strcmp(gas.stage(1,i).type,'exp')
        ip = i + 1;
        plot_hex(gas,[1,ip],fluidC(1),[1,1],100,16);
        break
    end
end

% PLOT HEX T-Q's (discharge)
% Hot HEX
for i=1:stages_dis
    if strcmp(gas.stage(4,i).type,'exp')
        ip = i - 1;
        plot_hex(gas,[4,ip],fluidH(1),[4,1],100,17);
        break
    end
end