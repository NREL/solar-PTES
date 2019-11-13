function [eff,DppH,DppC] = hex_analytic(fluidH, indH, fluidC, indC, HX)
% RESOLVE HEX T-Q DIAGRAM FOR A GIVEN EFFECTIVENESS

% Set number of sections for hex_core algorithm
NX = HX.NX;

% Set inlet temperatures (nomenclature: cold inlet is position 1, hot inlet
% is position 2)
TH2 = fluidH.state(indH(1),indH(2)).T;
TC1 = fluidC.state(indC(1),indC(2)).T;

% Check which one is fluidH and which is fluidC and swap them if necessary
if TC1 > TH2 % swap needed
    fluidH0 = fluidH;
    fluidH  = fluidC;
    fluidC  = fluidH0;
    indH0 = indH;
    indH  = indC;
    indC  = indH0;
end

% Import fluid.state and fluid.stage
stateH = fluidH.state(indH(1),indH(2));
stateC = fluidC.state(indC(1),indC(2));

% Set inlet pressures, enthalpies, entropies and mass flow rates
TH2 = stateH.T;
pH2 = stateH.p;
hH2 = stateH.h;
TC1 = stateC.T;
pC1 = stateC.p;
hC1 = stateC.h;
mH  = stateH.mdot;
mC  = stateC.mdot;

% Declare the two fluid streams
H = stream; H.mdot = mH; H.name = fluidH.name;
C = stream; C.mdot = mC; C.name = fluidC.name;

% Obtain temperature arrays as a function of the enthalpy arrays
[hvH,TvH] = get_h_T(fluidH,TC1-1,TH2+1,pH2,NX);
[hvC,TvC] = get_h_T(fluidC,TC1-1,TH2+1,pC1,NX);

% Obtain minimum and maximum enthalpy outlets (hot outlet cannot be colder
% than cold inlet, and vice-versa) and average specific heat capacity
hH1_min = rtab1(TvH,hvH,TC1,1);
hC2_max = rtab1(TvC,hvC,TH2,1);

% Compute average specific heat capacity estimate
CpHmean = (hH2 - hH1_min)/(TH2-TC1);
CpCmean = (hC2_max - hC1)/(TH2-TC1);

% Set enthalpy to estimate average value, and pressure to initial value
H.h = 0.5*(hH2 + hH1_min);
C.h = 0.5*(hC1 + hC2_max);
H.p = pH2;
C.p = pC1;

% Import values from HX structure into C and H structures
if H.p > C.p
    % Hot fluid flows inside the tubes
    H.D = HX.D1;
    H.G = HX.G1;
    H.A = HX.A1;
    C.D = HX.D2;
    C.G = HX.G2;
    C.A = HX.A2;
else
    % Hot fluid flows inside the shell side
    H.D = HX.D2;
    H.G = HX.G2;
    H.A = HX.A2;
    C.D = HX.D1;
    C.G = HX.G1;
    C.A = HX.A1;
end

% UPDATE PROPERTIES
% Cold stream
C = stream_update(fluidC,C,1);
% Hot stream
H = stream_update(fluidH,H,1);

% COMPUTE HEAT TRANSFER COEFFICIENTS
% Cold stream
C.Re = C.D*C.G./C.mu;
[C.Cf,C.St] = developed_flow(C.Re,C.Pr,HX.shape);
C.ht  = C.G*C.Cp.*C.St;
% Hot stream
H.Re = H.D*H.G./H.mu;
[H.Cf,H.St] = developed_flow(H.Re,H.Pr,HX.shape);
H.ht  = H.G*H.Cp.*H.St;
% Overall heat transfer coefficient (based on cold side heat transfer
% area). Neglects wall thermal resistance and axial conduction.
UlC  = 1./(C.A./(H.A*H.ht) + 1./C.ht);

% Compute the global Number of Transfer Units (NTU = UA/Cmin) and
% effectiveness
Cmin = min([CpCmean*mC,CpHmean*mH]);
Cmax = max([CpCmean*mC,CpHmean*mH]);
Cr   = Cmin/Cmax;
NTU  = UlC*C.A/Cmin;
if Cr < 0.999
    eff = ( 1 - exp(-NTU*(1-Cr)) ) / ( 1 - Cr*exp(-NTU*(1-Cr)) );
else % assume Cr=1
    eff = NTU/(1+NTU);
end

% Compute fractional pressure losses
DppH = 2*H.G^2*H.Cf*H.v*HX.L/(H.D*H.p);
DppC = 2*C.G^2*C.Cf*C.v*HX.L/(C.D*C.p);

% % Update outlet conditions and recompute average temperature and
% % properties!
% warning('update properties')
% keyboard

end



%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ hv, Tv ] = get_h_T( fluid, T1, T2, pressure, n )
% Obtain the hv and Tv arrays of a given fluid for the hex subroutines.
% Data ordered in regular intervals of h. T as a function of h.

if strcmp(fluid.read,'CP') %read from CoolProp
    
    h1  = CP1('PT_INPUTS',pressure,T1,'H',fluid.handle);
    h2  = CP1('PT_INPUTS',pressure,T2,'H',fluid.handle);
    hv  = linspace(h1,h2,n)';       % enthalpy array between TC1 and TH2
    pv  = ones(size(hv)).*pressure; % pressure array
    Tv  = CP1('HmassP_INPUTS',hv,pv,'T',fluid.handle); % temperature array
    
elseif strcmp(fluid.read,'TAB') %read from table
    
    Tx  = fluid.TAB(:,1);
    hy  = fluid.TAB(:,2);
    h1  = rtab1(Tx,hy,T1,0);
    h2  = rtab1(Tx,hy,T2,0);
    hv  = linspace(h1,h2,n)'; % enthalpy array between TC1 and TH2
    Tv  = rtab1(hy,Tx,hv,1);  % temperature array between TC1 and TH2
else
    error('not implemented')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%