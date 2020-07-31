function [ chi ] = RTeff( gas ,t_ch,t_dis)
%RTeff Round-trip efficiency
%   Adds up the contributions of the several stages to compute the round
%   trip efficiency

Work_in_ch = 0;
Work_out_dis = 0;

[~,n1] = size(gas.stage);
for i=1:n1
    Work_in_ch = Work_in_ch - gas.stage(1,i).w*gas.state(1,i).mdot*t_ch;
    Work_out_dis = Work_out_dis + gas.stage(2,i).w*gas.state(2,i).mdot*t_dis;
end

chi = Work_out_dis/Work_in_ch;

end

