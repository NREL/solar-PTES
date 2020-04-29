function [ S ] = developed_flow( S, mode )
% Compute heat transfer coefficient and friction factor

% Compute Reynolds number
S.Re = S.D*S.G./S.mu;

% First, compute the friction factor, Stanton number and heat transfer
% coefficient assuming single phase flow along the whole stream
[ S.Cf, S.St, S.ht ] = single_phase_flow( S.Re, S.Pr, S.G, S.Cp, S.shape );

% Compute pressure gradient due to flow friction. First, create average
% arrays of Cf and v
S.dpdL = - 2*S.G^2*S.Cf.*S.v./S.D;

% Check if any values fall within the two-phase region
itp = 0<=S.x & S.x<=1;
if any(itp)
    % Store array of vapour quality along two-phase region and compute
    % average
    x = S.x(itp);
    
    % Compute the single-phase friction factors and heat
    % transfer coefficients. For saturated liquid conditions:
    ReL  = S.D * S.G ./ S.muL;
    CpL  = S.PrL .* S.kL ./ S.muL;
    [ CfL, ~, htL ] = single_phase_flow( ReL, S.PrL, S.G, CpL, S.shape );
    % And saturated gas conditions:
    ReG  = S.D * S.G ./ S.muG;
    CpG  = S.PrG .* S.kG ./ S.muG;
    [ CfG, ~, ~   ] = single_phase_flow( ReG, S.PrG, S.G, CpG, S.shape );
    
    % Compute pressure gradients for single-phase conditions. Saturated
    % liquid:
    vL     = 1./S.rhoL;
    dpdL_L = - 2*S.G^2*CfL.*vL./S.D;
    % And saturated gas:
    vG     = 1./S.rhoG;
    dpdL_G = - 2*S.G^2*CfG.*vG./S.D;
    
    % Apply correlation from Müller-Steinhagen and Heck (1986):
    A = dpdL_L;
    B = dpdL_G;
    C = A + 2*(B-A).*x;
    dpdL_TP = C.*(1-x).^(1/3) + B.*x.^3;
    S.dpdL(itp) = dpdL_TP;
    
    if any(isnan(S.dpdL))
        keyboard
    end
    
    switch mode
        case 'heating'
            % Use Kandlikar's boiling heat transfer correlation
            
            % Set the fluid depending parameter
            switch S.name
                case 'Water'
                    F = 1.0;
                otherwise
                    error('not implemented')
            end
            
            % Compute the Froude number and the Froude depending function
            Fr   = S.G^2.0 ./ ( S.rhoL.^2.0 * 9.8 * S.D);
            f_Fr = max([ones(size(Fr)),2.63*Fr.^0.3],[],2);
            
            % Compute the Boiling number
            Bo   = S.q(itp) ./ (S.G.*S.hLG);
            
            % Compute the Convection number
            Co   = ((1 - x)./x).^0.8 .* (S.rhoG./S.rhoL).^0.5;
            
            % Compute the heat transfer coefficient in the nucleate boiling region
            % and in the convective region
            hNBR = htL.*( 0.6683*Co.^(-0.2).*f_Fr + 1058.0*Bo.^0.7*F ).*(1-x).^0.8;
            hCR  = htL.*( 1.1360*Co.^(-0.9).*f_Fr +  667.2*Bo.^0.7*F ).*(1-x).^0.8;
            
            % The actual heat transfer coefficient is the maximum of the two
            S.ht(itp) = max([hNBR,hCR],[],2);
            
        case 'cooling'
            % Use Shah's heat transfer correlation
            
            % Compute reduced pressure
            Pcrit = RP1(0,0,0,'Pcrit',S);
            pr    = S.p/Pcrit; %reduced pressure
            
            % Apply correlation
            htTP = htL.*( (1-x).^0.8 + (3.8*x.^0.76.*(1-x).^0.04)./(pr.^0.38));
            S.ht(itp) = htTP;
            
        otherwise
            error('not implemented')
    end
    
    if any(isnan(S.ht)) || any(isnan(S.dpdL))
        warning('NaN value found inside developed_flow function')
        keyboard
    end
    
    % Set Stanton number to NaN for two-phase region (where it is not well
    % defined)
    S.St(itp) = NaN;
end

end

%%% SUPPORT FUNCTIONS %%%
function [ Cf, St, ht ] = single_phase_flow( Re, Pr, G, Cp, shape )

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
% Heat transfer coefficient
ht = G*Cp.*St;

end
