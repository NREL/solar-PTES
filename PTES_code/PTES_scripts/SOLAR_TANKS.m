% Set charged conditions of fluid tanks (i.e. from unspecified solar or
% electrical source)

if setTmax
else
    error('***Heat engine mode requires setTmax=1 under current implementation***')
end

HT = reset_tanks(HT,273+250,p0,0,Tmax,p0,1e6,T0);
HT.A(2) = HT.A(1);
HT.B(2) = HT.B(1);

% Set PRch used for reference to set PR_dis range
PRch = PR_estim;