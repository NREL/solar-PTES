% SET_DESIGN
% The point of this script is to extract key information once the design
% point load cycle has been run. This information is used to set the design
% point of the compressors, expanders, and heat exchangers. (This
% information is also used to calculate the capital cost of the system).

% TURN OFF SETTMAX
setTmax = 0 ; % Instead everything is based off the pressure ratio
pmax    = CCMP.pr0 * CCMP.Pin ;

% Worth resetting a few components
% Reset cold tanks
for ir = 1 : Ncld
    CT(ir) = reset_tanks(CT(ir),TC_dis0(ir),10.*p0,MC_dis0(ir),TC_chg0(ir),10.*p0,MC_chg0(ir),T0);
end

% Reset hot tanks
for ir = 1 : Nhot
    HT(ir) = reset_tanks(HT(ir),TH_dis0(ir),10.*p0,MH_dis0(ir),TH_chg0(ir),10.*p0,MH_chg0(ir),T0);
end

% Reset atmospheric tanks
AT = reset_tanks(AT,T0,p0,huge,T0,p0,huge,T0);

% Reset the Load parameter to be the off-design load cycle
Load = Load0 ;

design_mode = 0 ; % Logical - no longer in design mode