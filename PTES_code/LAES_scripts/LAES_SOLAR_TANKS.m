% Set charged conditions of fluid tanks (i.e. from unspecified solar or
% electrical source)

% Hot tank
HT.B(3).T = Tmax;
HT.B(3).p = p0;
HT.B(3).M = 0.26e6;
HT.B(3)   = update_tank(HT.B(3),fluidH(1),T0,1);