function [ Cf, St ] = developed_flow( Re, Pr, shape )
%Compute friction factor and Stanton number from Reynolds and Prandtl
%numbers

if strcmp(shape,'circular')
    % Set coefficients for circular pipe
    aL  = 16.0;
    bL  = 4.36;
else
    error('not implemented')
end

% Set laminar to turbulent transition points
Re_lim1 = 2000; % first limit (laminar to transition)
Re_lim2 = 3000; % second limit (transition to turbulence)

% Obtain array indices for laminar, transition and turbulent regimes
lam = find(Re<=Re_lim1);
tra = find(Re>Re_lim1 & Re<=Re_lim2);
tur = find(Re>Re_lim2);

% Compute the transition "weight" parameter, which is used to compute
% parameters in the transition regime. 0<W<1 (0=laminar,1=turbulent)
W = (Re(tra) - Re_lim1)/(Re_lim2 - Re_lim1);

% Allocate matrices
f  = zeros(size(Re));
Nu = zeros(size(Re));

% Compute friction factor and Stanton numbers. Fully developed flow model
% with constant heat flux (Incropera):
% Friction factor
f(lam)  = aL*4./Re(lam);
f(tur)  = (0.790 * log(Re(tur)) - 1.64).^(-2);
f_lim1  = aL*4./Re_lim1;
f_lim2  = (0.790 * log(Re_lim2) - 1.64).^(-2);
f(tra)  = (1-W)*f_lim1 + W*f_lim2;
Cf = f/4;
% Nusselt number and Stanton number
Nu(lam) = bL + 0.*Re(lam);
Nu(tur) = (f(tur)/8).*(Re(tur) - 1000).*Pr(tur)./(1 + 12.7*(f(tur)/8).^(1/2).*(Pr(tur).^(2/3)-1));
Nu_lim1 = bL + 0.*Re_lim1;
Nu_lim2 = (f_lim2/8).*(Re_lim2 - 1000).*Pr(tra)./(1 + 12.7*(f_lim2/8).^(1/2).*(Pr(tra).^(2/3)-1));
Nu(tra) = (1-W)*Nu_lim1 + W.*Nu_lim2;
St = Nu./(Re.*Pr);

if any(isnan([Cf;St]))
    warning('error in developed_flow function')
    keyboard
end

% figure(3)
% plot(Re,f)
% xlabel('Reynolds number')
% ylabel('Friction factor')
% 
% figure(4)
% plot(Re,Nu)
% xlabel('Reynolds number')
% ylabel('Nusselt number')

end

