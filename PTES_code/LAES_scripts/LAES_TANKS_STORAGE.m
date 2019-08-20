% At this point, assume no storage losses

% Hot tanks
HT.A(3)  = HT.A(2);  % charge source
HT.B(3)  = HT.B(2);  % charge sink

% Cold tanks
CT1.A(3) = CT1.A(2); % charge source
CT1.B(3) = CT1.B(2); % charge sink
CT2.A(3) = CT2.A(2); % charge source
CT2.B(3) = CT2.B(2); % charge sink

% Air sources/sinks
LA(3)    = LA(2);
AA(3)    = AA(2);