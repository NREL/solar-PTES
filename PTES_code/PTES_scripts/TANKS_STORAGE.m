% At this point, assume no storage losses

% Hot tanks
for ii = 1 : Nhot
    HT(ii).A(iL+1) = HT(ii).A(iL); % charge source

    HT(ii).A(iL+1).T = HT(ii).A(iL).T + Load.HT_A(iL);
    HT(ii).A(iL+1)   = update_tank_state(HT(ii),HT(ii).A(iL+1),T0,1);

    HT(ii).B(iL+1) = HT(ii).B(iL); % charge sink

    HT(ii).B(iL+1).T = HT(ii).B(iL).T + Load.HT_B(iL);
    HT(ii).B(iL+1)   = update_tank_state(HT(ii),HT(ii).B(iL+1),T0,1);
    
end

% Cold tanks
for ii = 1 : Ncld
    CT(ii).A(iL+1) = CT(ii).A(iL); % charge source

    CT(ii).A(iL+1).T = CT(ii).A(iL).T + Load.CT_A(iL);
    CT(ii).A(iL+1)   = update_tank_state(CT(ii),CT(ii).A(iL+1),T0,1);

    CT(ii).B(iL+1) = CT(ii).B(iL); % charge sink

    CT(ii).B(iL+1).T = CT(ii).B(iL).T + Load.CT_B(iL);
    CT(ii).B(iL+1)   = update_tank_state(CT(ii),CT(ii).B(iL+1),T0,1);

end

% Atmospheric tanks
AT.A(iL+1) = AT.A(iL); % charge source
AT.B(iL+1) = AT.B(iL); % charge sink

