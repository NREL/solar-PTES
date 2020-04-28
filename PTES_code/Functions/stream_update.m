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
        
    case 2 %Enthalpy and pressure are known
        [S.T]   = RP1('HmassP_INPUTS',S.h,S.p,'T',S);
        [S.rho] = RP1('HmassP_INPUTS',S.h,S.p,'D',S);
        [S.k]   = RP1('HmassP_INPUTS',S.h,S.p,'CONDUCTIVITY',S);
        [S.mu]  = RP1('HmassP_INPUTS',S.h,S.p,'VISCOSITY',S);
        [S.Pr]  = RP1('HmassP_INPUTS',S.h,S.p,'PRANDTL',S);
        [S.x]   = RP1('HmassP_INPUTS',S.h,S.p,'Q',S);
        
    otherwise
        error('not implemented')
end

% Check if any values fall inside the two-phase region. If so,
% extract properties corresponding to the saturated liquid
% condition for that pressure.
itp = 0<=S.x & S.x<=1;
if any(itp)
    % Allocate arrays
    S.hLG  = zeros(size(S.p));
    S.rhoL = zeros(size(S.p));
    S.rhoG = zeros(size(S.p));
    S.kL   = zeros(size(S.p));
    S.muL  = zeros(size(S.p));
    S.PrL  = zeros(size(S.p));
    
    % Enthalpy of vaporisation
    [S.hLG(itp)] = RP1('PQ_INPUTS',S.p(itp),1.0*ones(size(S.p(itp))),'H',S)...
        -RP1('PQ_INPUTS',S.p(itp),0.0*ones(size(S.p(itp))),'H',S);
    
    % Density at saturated conditions
    [S.rhoL(itp)] = RP1('PQ_INPUTS',S.p(itp),0.0*ones(size(S.p(itp))),'D',S);
    [S.rhoG(itp)] = RP1('PQ_INPUTS',S.p(itp),1.0*ones(size(S.p(itp))),'D',S);
    
    % Conductivity, viscosity and Prandtl numbers at saturated
    % liquid conditions
    [S.kL(itp)]   = RP1('PQ_INPUTS',S.p(itp),0.0*ones(size(S.p(itp))),'CONDUCTIVITY',S);
    [S.muL(itp)]  = RP1('PQ_INPUTS',S.p(itp),0.0*ones(size(S.p(itp))),'VISCOSITY',S);
    [S.PrL(itp)]  = RP1('PQ_INPUTS',S.p(itp),0.0*ones(size(S.p(itp))),'PRANDTL',S);
end

% Compute derived properties
S.v  = 1./S.rho;
S.Cp = S.Pr.*S.k./S.mu;

%%% TO BE IMPROVED
%%% ------------->
% CORRECT OUTLIERS FROM mu, k, Cp AND Pr ARRAYS. NOTE: THIS IS A PATCH
% CURRENTLY USED TO BE ABLE TO RUN THE CODE WITHOUT CRASHING WHEN A HEX IS
% OPERATING WITH A TWO-PHASE FLOW STREAM. HOWEVER, THIS IS CURRENTLY
% PHYSICALLY INCORRECT (Cp and Pr are not well defined inside the
% saturation curve) AND MUST BE ADDRESSED PROPERLY.
cond1 = S.mu == 0;
if any(cond1)
    S.mu(cond1) = NaN;
    S.mu = fillmissing(S.mu,'nearest');
end
cond2 = S.k == 0;
if any(cond2)
    S.k(cond2) = NaN;
    S.k  = fillmissing(S.k,'nearest');
end
Cp_mean = median(S.Cp);
Pr_mean = median(S.Pr);
cond3 = S.Cp <= 0 | S.Cp > 100*Cp_mean | S.Cp < Cp_mean/100 | isnan(S.Cp);
cond4 = S.Pr <= 0 | S.Pr > 100*Pr_mean | S.Pr < Pr_mean/100 | isnan(S.Pr);
if any([cond3;cond4])
    %{
    figure(31)
    plot(S.T)
    figure(32)
    semilogy(S.Pr)
    figure(33)
    semilogy(S.Cp)
    %}
    Cp_mean = median(S.Cp(~cond3));
    S.Cp(cond3) = Cp_mean;
    Pr_mean = median(S.Pr(~cond4));
    S.Pr(cond4) = Pr_mean;
    %{
    figure(34)
    semilogy(S.Pr)
    figure(35)
    semilogy(S.Cp)
    keyboard
    %}
end
%%% <--------------

end

