% Set charged conditions of fluid tanks (i.e. from unspecified solar or
% electrical source)

if setTmax
    HT = reset_tanks(HT,273+250,p0,0,Tmax,p0,1e6,T0);
    % Set PRch used for reference to set PR_dis range
    %PRch = PR_estim;
    HT.A(2) = HT.A(1);
    HT.B(2) = HT.B(1);
else
    switch Load.mode
        case {5,6}
            HT(1) = reset_tanks(HT(1),TH_dis0(1),p0,0,TH_chg0(1),p0,1e6,T0);
    
            HT(1).A(2) = HT(1).A(1);
            HT(1).B(2) = HT(1).B(1);
            
        case {0,1,2,3,4}
            error('Not implemented')
    end
end



