% Script that calculates the cost of each component and subsequently the
% total cost of the system, as well as other economic metrics (possibly
% ...)


% Compressors and expanders
for ii = 1 : Nc_ch
    CCMP(ii) = compexp_econ(CCMP(ii))  ;
    DEXP(ii) = compexp_econ(DEXP(ii))  ;
end

for ii = 1 : Ne_ch
    CEXP(ii) = compexp_econ(CEXP(ii))  ;
    DCMP(ii) = compexp_econ(DCMP(ii))  ;
end