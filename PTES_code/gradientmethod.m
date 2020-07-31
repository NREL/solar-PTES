clc
clear all

x0=[723  1.6  0.97]; % Initial guess of the Parameters [TH_dis0  PRr  eff]
    
xl=[523  0.8  0.8];  % lower bound parameter values
xu=[723  1.6  0.97]; % upper bound parameter values

eq1=[];%No equality or inequality constraints
eq2=[];
eq3=[];
eq4=[];

ff=@PTES_opt %PTES cycle
[Parameters, objective]=fmincon(ff,x0,eq1,eq2,eq3,eq4,xl,xu)

