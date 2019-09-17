% At this point, assume no storage losses

% Hot tanks
HT.A(iL+1) = HT.A(iL); % charge source
HT.B(iL+1) = HT.B(iL); % charge sink

if Nhot == 2
    HT2.A(iL+1) = HT2.A(iL); % charge source
    HT2.B(iL+1) = HT2.B(iL); % charge sink    
end

% Cold tanks
CT.A(iL+1) = CT.A(iL); % charge source
CT.B(iL+1) = CT.B(iL); % charge sink

if Ncld == 2
    CT2.A(iL+1) = CT2.A(iL); % charge source
    CT2.B(iL+1) = CT2.B(iL); % charge sink    
end