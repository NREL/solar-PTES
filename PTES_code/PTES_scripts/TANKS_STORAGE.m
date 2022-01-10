% At this point, assume no storage losses

% Hot tanks
for ii = 1 : Nhot
    HT(ii).A(iL+1) = HT(ii).A(iL); % charge source
    HT(ii).B(iL+1) = HT(ii).B(iL); % charge sink
end

% Middle tanks
for ii = 1 : Nmid
    MT(ii).A(iL+1) = MT(ii).A(iL); % charge source
    MT(ii).B(iL+1) = MT(ii).B(iL); % charge sink
end

% Cold tanks
for ii = 1 : Ncld
    CT(ii).A(iL+1) = CT(ii).A(iL); % charge source
    CT(ii).B(iL+1) = CT(ii).B(iL); % charge sink
end

% Atmospheric tanks
AT.A(iL+1) = AT.A(iL); % charge source
AT.B(iL+1) = AT.B(iL); % charge sink

