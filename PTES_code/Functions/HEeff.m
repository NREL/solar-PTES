function [ chi ] = HEeff( gas, BH, t_dis)
%HEeff Heat engine efficiency
%   Adds up the contributions of the several stages to compute the heat
%   engine exergetic efficiency

Work_out_dis = 0;

[~,n1] = size(gas.stage);
for i=1:n1
    Work_out_dis = Work_out_dis + gas.stage(2,i).w*gas.state(2,i).mdot*t_dis;
end

chi = Work_out_dis/BH;

end

