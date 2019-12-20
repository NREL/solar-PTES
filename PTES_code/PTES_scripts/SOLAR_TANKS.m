% Set charged conditions of fluid tanks (i.e. from unspecified solar or
% electrical source)

if setTmax
    HT = reset_tanks(HT,273+250,p0,0,Tmax,p0,1e6,T0);
    % Set PRch used for reference to set PR_dis range
    PRch = PR_estim;
else
    HT = reset_tanks(HT,TH_dis0,p0,0,TH_chg0,p0,1e6,T0);
end

HT.A(2) = HT.A(1);
HT.B(2) = HT.B(1);

