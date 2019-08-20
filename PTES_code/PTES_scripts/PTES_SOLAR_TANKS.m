% Set charged conditions of fluid tanks (i.e. from unspecified solar or
% electrical source)

% Hot tank
if setTmax
else
    error('***Heat engine mode requires setTmax=1 under current implementation***')
end

HT.B(1).T = Tmax;
HT.B(1).p = p0;
HT.B(1).M = 0;
HT.B(1)   = update_tank(HT.B(1),fluidH(1),T0,1);
HT.B(2)   = HT.B(1);
HT.B(2).M = 0.26e6;
HT.B(2)   = update_tank(HT.B(2),fluidH(1),T0,1);
HT.B(3)   = HT.B(2);

HT.A(1).T = TH_0;
HT.A(1).p = p0;
HT.A(1).M = 0.26e6;
HT.A(1)   = update_tank(HT.A(1),fluidH(1),T0,1);
HT.A(2)   = HT.A(1);
HT.A(2).M = 0;
HT.A(2)   = update_tank(HT.A(2),fluidH(1),T0,1);
HT.A(3)   = HT.A(2);

% Set PRch used for reference to set PR_dis range
PRch = PR_estim;