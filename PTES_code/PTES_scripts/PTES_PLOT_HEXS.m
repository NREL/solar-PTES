% PLOT HEX T-Q's (charge)
% Hot HEX
ip = 2;
plot_hex(gas,[1,ip],fluidH(1),[1,1],100,4);
% Cold HEX
for i=1:stages_ch
    if strcmp(gas.stage(1,i).type,'exp')
        ip = i + 1;
        plot_hex(gas,[1,ip],fluidC(1),[1,1],100,6);
        break
    end
end

% PLOT HEX T-Q's (charge)
% % Hot HEX
% ip = 2;
% plot_hex(gas,[2,ip],fluidH(1),[1,1],100,5);
% Hot HEX
for i=1:stages_dis
    if strcmp(gas.stage(2,i).type,'exp')
        ip = i - 1;
        plot_hex(gas,[2,ip],fluidH(1),[2,1],100,7);
        break
    end
end