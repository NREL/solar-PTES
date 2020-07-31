function [TH1,TC2,pH1,pC2] = hex_analytic(HX, iL, fluidH, iH, fluidC, iC)
% Solve T-Q diagram using analytical solutions

% Set inlet temperatures (nomenclature: cold inlet is position 1, hot inlet
% is position 2)
TH2 = fluidH.state(iL,iH).T;
TC1 = fluidC.state(iL,iC).T;

% Check which one is fluidH and which is fluidC and swap them if necessary
if TC1 > TH2 % swap needed
    fluidH0 = fluidH;
    fluidH  = fluidC;
    fluidC  = fluidH0;
    indH0 = iH;
    iH  = iC;
    iC  = indH0;
end

% Import fluid.state and fluid.stage
stateH = fluidH.state(iL,iH);
stateC = fluidC.state(iL,iC);

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
H.read = fluidH.read; H.handle = fluidH.handle; H.HEOS = fluidH.HEOS;
C.read = fluidC.read; C.handle = fluidC.handle; C.HEOS = fluidC.HEOS;
H.shape = HX.shape;
C.shape = HX.shape;

% Obtain minimum and maximum enthalpy outlets (hot outlet cannot be colder
% than cold inlet, and vice-versa) and average specific heat capacity
hH1_min = RPN('PT_INPUTS',pH2,TC1,'H',fluidH);
hC2_max = RPN('PT_INPUTS',pC1,TH2,'H',fluidC);

% Set enthalpy to estimate average value, and pressure to initial value
H.h = 0.5*(hH2 + hH1_min);
C.h = 0.5*(hC1 + hC2_max);
H.p = pH2;
C.p = pC1;

% Import values from HX structure into C and H structures
if H.p > C.p
    % Hot fluid flows inside the tubes
    H.D  = HX.D1;
    H.Af = HX.Af1;
    H.A  = HX.A1;
    C.D  = HX.D2;
    C.Af = HX.Af2;
    C.A  = HX.A2;
else
    % Hot fluid flows inside the shell side
    H.D  = HX.D2;
    H.Af = HX.Af2;
    H.A  = HX.A2;
    C.D  = HX.D1;
    C.Af = HX.Af1;
    C.A  = HX.A1;
end

% Compute mass fluxes
H.G = mH/H.Af;
C.G = mC/C.Af;

% Initial guess
hH1 = hH1_min;
TH1 = TC1;
hC2 = hC2_max;
TC2 = TH2;

for i=1:3
    % UPDATE PROPERTIES
    % Cold stream
    C = stream_update(C,2);
    % Hot stream
    H = stream_update(H,2);
    
    % COMPUTE HEAT TRANSFER COEFFICIENTS
    % Cold stream
    [ C ] = developed_flow( C, 'heating' );
    % Hot stream
    [ H ] = developed_flow( H, 'cooling' );
    % Overall heat transfer coefficient (based on cold side heat transfer
    % area). Neglects wall thermal resistance and axial conduction.
    UlC  = 1./(C.A./(H.A*H.ht) + 1./C.ht);
    
    % Compute the global Number of Transfer Units (NTU = UA/Cmin). Start by
    % computing the average specific heat capacity and the minimum and maximum
    % Cp*mdot product.
    CpHmean = (hH2 - hH1)/(TH2-TH1);
    CpCmean = (hC2 - hC1)/(TC2-TC1);
    Cmin = min([CpCmean*mC,CpHmean*mH]);
    Cmax = max([CpCmean*mC,CpHmean*mH]);
    Cr   = Cmin/Cmax;
    NTU  = UlC*C.A/Cmin;
    
    % Compute the heat exchanger effectiveness using analytical expressions for
    % a pure counter-flow heat exchanger.
    if Cr < 0.999
        eff = ( 1 - exp(-NTU*(1-Cr)) ) / ( 1 - Cr*exp(-NTU*(1-Cr)) );
    else % assume Cr=1
        eff = NTU/(1+NTU);
    end
    
    % Compute fractional pressure losses
    DppH = 2*H.G^2*H.Cf*H.v*HX.L/(H.D*H.p);
    DppC = 2*C.G^2*C.Cf*C.v*HX.L/(C.D*C.p);
        
    % Compute QMAX and the actual total heat transfer, QT, based on the
    % predicted heat exchanger effectiveness. Use QT to obtain the
    % outlet enthalpies and temperatures
    QMAX = min([mC*(hC2_max - hC1),mH*(hH2 - hH1_min)]);
    QT   = eff*QMAX;
    hC2  = hC1 + QT/mC;
    hH1  = hH2 - QT/mH;
    pC2  = pC1*(1 - DppC);
    pH1  = pH2*(1 - DppH);
    TC2 = RPN('HmassP_INPUTS',hC2,pC2,'T',fluidC);
    TH1 = RPN('HmassP_INPUTS',hH1,pH1,'T',fluidH);
    
    % Update average properties
    H.h = 0.5*(hH2 + hH1);
    C.h = 0.5*(hC1 + hC2);
    H.p = pH2 - 0.5*DppH;
    C.p = pC1 - 0.5*DppC;
end

end