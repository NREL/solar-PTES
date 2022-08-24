% At this point, assume no storage losses

% Hot tanks
for ii = 1 : Nhot

    [HT(ii), HT(ii).A]  = store_tanks(HT(ii) , HT(ii).A , iL , Load.HT_A(iL) , T0); % hot source

    [HT(ii), HT(ii).B]  = store_tanks(HT(ii) , HT(ii).B , iL , Load.HT_B(iL) , T0); % hot sink
    
end

% Cold tanks
for ii = 1 : Ncld

    [CT(ii), CT(ii).A]  = store_tanks(CT(ii) , CT(ii).A , iL , Load.CT_A(iL) , T0); % cold source

    [CT(ii), CT(ii).B]  = store_tanks(CT(ii) , CT(ii).B , iL , Load.CT_B(iL) , T0); % cold sink

end

% Atmospheric tanks
AT.A(iL+1) = AT.A(iL); % charge source
AT.B(iL+1) = AT.B(iL); % charge sink

