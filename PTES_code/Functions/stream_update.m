function [S] = stream_update(S,mode)
%STREAM_UPDATE Updates the thermophysical properties along to HEX streams.
%
%   E.g. Given two arrays of pressure and enthalpy (mode==2), it obtains
%   the corresponding arrays of T, rho, k, mu and Pr. STREAM_UPDATE calls
%   the function RP1.

% Obtain main thermophysical properties
switch mode
    case 1 %Temperature and pressure are known
        [S.h]   = RPN('PT_INPUTS',S.p,S.T,'H',S);
        [S.rho] = RPN('PT_INPUTS',S.p,S.T,'D',S);
        [S.k]   = RPN('PT_INPUTS',S.p,S.T,'CONDUCTIVITY',S);
        [S.mu]  = RPN('PT_INPUTS',S.p,S.T,'VISCOSITY',S);
        [S.Pr]  = RPN('PT_INPUTS',S.p,S.T,'PRANDTL',S);
        [S.x]   = RPN('PT_INPUTS',S.p,S.T,'Q',S);
        error('not implemented');
        
    case 2 %Enthalpy and pressure are known
        
        % Preallocate arrays
        sz    = size(S.p);
        S.T   = zeros(sz);
        S.rho = zeros(sz);
        S.k   = zeros(sz);
        S.mu  = zeros(sz);
        S.Pr  = zeros(sz);
        
        % Detect points falling inside the two-phase region
        S.x = RPN('HmassP_INPUTS',S.h,S.p,'Q',S);
        tp  = 0<=S.x & S.x<=1; % two-phase index
        sp  = ~tp; % single-phase index
        
        % Obtain properties. Obtain values of T and rho for all points.
        [S.T,S.rho] = RPN('HmassP_INPUTS',S.h,S.p,{'T','D'},S);
        
        % Obtain values of k, mu and Pr only for points in single-phase
        % regions. Set other points to NaN.
        [S.k(sp),S.mu(sp),S.Pr(sp)]  = ...
            RPN('HmassP_INPUTS',S.h(sp),S.p(sp),{'CONDUCTIVITY','VISCOSITY','PRANDTL'},S);
        S.k(tp)  = NaN;
        S.mu(tp) = NaN;
        S.Pr(tp) = NaN;
        
    otherwise
        error('not implemented')
end

% For any values falling inside the two-phase region, extract properties
% corresponding to the saturated liquid and gas conditions at that
% pressure.
if any(tp)
    % Set size of arrays
    sz = size(S.p(tp));
    
    % Saturated liquid
    [S.rhoL,S.kL,S.muL,S.PrL] =...
        RPN('PQ_INPUTS',S.p(tp),0.0*ones(sz),{'D','CONDUCTIVITY','VISCOSITY','PRANDTL'},S);
    
    % Saturated gas
    [S.rhoG,S.kG,S.muG,S.PrG] =...
        RPN('PQ_INPUTS',S.p(tp),1.0*ones(sz),{'D','CONDUCTIVITY','VISCOSITY','PRANDTL'},S);
    
    % Enthalpy of vaporisation
    S.hLG  = RPN('PQ_INPUTS',S.p(tp),1.0*ones(sz),'H',S)...
        -RPN('PQ_INPUTS',S.p(tp),0.0*ones(sz),'H',S);
end

% CORRECT NUL COOLPROP VALUES (this seems to happen close to the
% saturation curve)
nul = S.k == 0;
if any(nul)
    % Determine phase
    Pcrit = RPN(0,0,0,'Pcrit',S);
    Tsat  = RPN('PQ_INPUTS',S.p,0.0*ones(size(S.p)),'T',S);
    scrit = S.p > Pcrit;
    liq   = ~scrit & S.T < Tsat;
    gas   = ~scrit & S.T > Tsat;
    
    % Fill null values with closest saturated value
    if any(nul&liq)
        [S.k(nul&liq),S.mu(nul&liq),S.Pr(nul&liq)] =...
            RPN('PQ_INPUTS',S.p(nul&liq),0.0*ones(size(S.p(nul&liq))),{'CONDUCTIVITY','VISCOSITY','PRANDTL'},S);
    end
    if any(nul&gas)
        [S.k(nul&gas),S.mu(nul&gas),S.Pr(nul&gas)] =...
            RPN('PQ_INPUTS',S.p(nul&gas),1.0*ones(size(S.p(nul&gas))),{'CONDUCTIVITY','VISCOSITY','PRANDTL'},S);
    end
    if any(nul&scrit)
        error('not implemented')
    end
end

% Compute derived properties
S.v  = 1./S.rho;
S.Cp = S.Pr.*S.k./S.mu;
end

