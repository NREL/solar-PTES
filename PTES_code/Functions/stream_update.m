function [S] = stream_update(S,mode)
%STREAM_UPDATE Updates the thermophysical properties along to HEX streams.
%
%   E.g. Given two arrays of pressure and enthalpy (mode==2), it obtains
%   the corresponding arrays of T, rho, k, mu and Pr. STREAM_UPDATE calls
%   the function RP1.

% Obtain main thermophysical properties
switch mode
    case 1 %Temperature and pressure are known
        [S.h]   = RP1('PT_INPUTS',S.p,S.T,'H',S);
        [S.rho] = RP1('PT_INPUTS',S.p,S.T,'D',S);
        [S.k]   = RP1('PT_INPUTS',S.p,S.T,'CONDUCTIVITY',S);
        [S.mu]  = RP1('PT_INPUTS',S.p,S.T,'VISCOSITY',S);
        [S.Pr]  = RP1('PT_INPUTS',S.p,S.T,'PRANDTL',S);
        [S.x]   = RP1('PT_INPUTS',S.p,S.T,'Q',S);
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
        S.x = RP1('HmassP_INPUTS',S.h,S.p,'Q',S);
        tp  = 0<=S.x & S.x<=1; % two-phase index
        sp  = ~tp; % single-phase index
        
        % Obtain properties. Obtain values of T and rho for all points.
        S.T   = RP1('HmassP_INPUTS',S.h,S.p,'T',S);
        S.rho = RP1('HmassP_INPUTS',S.h,S.p,'D',S);
        
        % Obtain values of k, mu and Pr only for points in single-phase
        % regions. Set other points to NaN.
        S.k(sp)  = RP1('HmassP_INPUTS',S.h(sp),S.p(sp),'CONDUCTIVITY',S);
        S.mu(sp) = RP1('HmassP_INPUTS',S.h(sp),S.p(sp),'VISCOSITY',S);
        S.Pr(sp) = RP1('HmassP_INPUTS',S.h(sp),S.p(sp),'PRANDTL',S);
        S.k(tp)  = NaN;
        S.mu(tp) = NaN;
        S.Pr(tp) = NaN;
        
    otherwise
        error('not implemented')
end

% Any any values falling inside the two-phase region, extract properties
% corresponding to the saturated liquid and gas conditions at that
% pressure.
if any(tp)
    % Set size of arrays
    sz = size(S.p(tp));
    
    % Saturated liquid
    S.rhoL = RP1('PQ_INPUTS',S.p(tp),0.0*ones(sz),'D',S);
    S.kL   = RP1('PQ_INPUTS',S.p(tp),0.0*ones(sz),'CONDUCTIVITY',S);
    S.muL  = RP1('PQ_INPUTS',S.p(tp),0.0*ones(sz),'VISCOSITY',S);
    S.PrL  = RP1('PQ_INPUTS',S.p(tp),0.0*ones(sz),'PRANDTL',S);
    
    % Saturated gas
    S.rhoG = RP1('PQ_INPUTS',S.p(tp),1.0*ones(sz),'D',S);
    S.kG   = RP1('PQ_INPUTS',S.p(tp),1.0*ones(sz),'CONDUCTIVITY',S);
    S.muG  = RP1('PQ_INPUTS',S.p(tp),1.0*ones(sz),'VISCOSITY',S);
    S.PrG  = RP1('PQ_INPUTS',S.p(tp),1.0*ones(sz),'PRANDTL',S);
    
    % Enthalpy of vaporisation
    S.hLG  = RP1('PQ_INPUTS',S.p(tp),1.0*ones(sz),'H',S)...
        -RP1('PQ_INPUTS',S.p(tp),0.0*ones(sz),'H',S);
end

% CORRECT NULL COOLPROP VALUES (this seems to happen close to the
% saturation curve)
nul = S.k == 0;
if any(nul)
    % Determine phase
    Pcrit = RP1(0,0,0,'Pcrit',S);
    Tsat  = RP1('PQ_INPUTS',S.p,0.0*ones(size(S.p)),'T',S);
    scrit = S.p > Pcrit;
    liq   = ~scrit & S.T < Tsat;
    gas   = ~scrit & S.T > Tsat;
    
    % Fill null values with closest saturated value
    if any(nul&liq)
        S.k(nul&liq)  = RP1('PQ_INPUTS',S.p(nul&liq),0.0*ones(size(S.p(nul&liq))),'CONDUCTIVITY',S);
        S.mu(nul&liq) = RP1('PQ_INPUTS',S.p(nul&liq),0.0*ones(size(S.p(nul&liq))),'VISCOSITY',S);
        S.Pr(nul&liq) = RP1('PQ_INPUTS',S.p(nul&liq),0.0*ones(size(S.p(nul&liq))),'PRANDTL',S);
    end
    if any(nul&gas)
        S.k(nul&gas)  = RP1('PQ_INPUTS',S.p(nul&gas),1.0*ones(size(S.p(nul&gas))),'CONDUCTIVITY',S);
        S.mu(nul&gas) = RP1('PQ_INPUTS',S.p(nul&gas),1.0*ones(size(S.p(nul&gas))),'VISCOSITY',S);
        S.Pr(nul&gas) = RP1('PQ_INPUTS',S.p(nul&gas),1.0*ones(size(S.p(nul&gas))),'PRANDTL',S);
    end
    if any(nul&scrit)
        error('not implemented')
    end
end

% Compute derived properties
S.v  = 1./S.rho;
S.Cp = S.Pr.*S.k./S.mu;
end

