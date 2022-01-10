function [ S ] = developed_flow( S, mode )
% Compute heat transfer coefficient and friction factor

% Compute Reynolds number
S.Re = S.D*S.G./S.mu;

% Compute Graetz number
if ~isempty(S.LS)
    S.Gz = (S.D./S.LS).*S.Re.*S.Pr;
else
    S.Gz = zeros(size(S.Re));
end

% First, compute the friction factor, Stanton number and heat transfer
% coefficient assuming single phase flow along the whole stream
[ S.Cf, S.St, S.ht ] = single_phase_flow( S.Re, S.Pr, S.G, S.Cp, S.shape, S.Gz);

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
    if ~isempty(S.LS)
        GzL = (S.D./S.LS(itp)).*ReL.*S.PrL;
    else
        GzL = zeros(size(ReL));
    end
    [ CfL, ~, htL ] = single_phase_flow( ReL, S.PrL, S.G, CpL, S.shape, GzL );
    % And saturated gas conditions:
    ReG  = S.D * S.G ./ S.muG;
    CpG  = S.PrG .* S.kG ./ S.muG;
    if ~isempty(S.LS)
        GzG = (S.D./S.LS(itp)).*ReG.*S.PrG;
    else
        GzG = zeros(size(ReG));
    end
    [ CfG, ~, ~   ] = single_phase_flow( ReG, S.PrG, S.G, CpG, S.shape, GzG );
    
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
            Pcrit = RPN(0,0,0,'Pcrit',S);
            p     = S.p(itp);
            pr    = 0.5*p/Pcrit; %reduced pressure
            
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
function [ Cf, St, ht ] = single_phase_flow( Re, Pr, G, Cp, shape, Gz )

if strcmp(shape,'cross-flow')
    % Finned tubes (external transversal flow). Cf and St are only a
    % function of Re and Pr. There is no transition points
    
    % Compute Colburn factor and heat transfer coefficient
    j  = 0.1733 * Re.^(-0.4071);
    St = j./Pr.^(2/3);
    ht = G*Cp.*St;
    
    % Fanning friction factor (friction coefficient)
    Cf = 0.1275 * Re.^(-0.2128);
    
    return
end

% Set regime transition points. 'lim1' is the first limit (laminar to
% transition). 'lim2' is the second limit (transition to turbulence).
switch shape
    case 'circular' % Macroscopic circular pipe
        Re_lim1 = 2000;
        Re_lim2 = 3000;
        
    case 'PCHE' % Semi-circular PCHE (straight channels)
        Re_lim1 = 2000;
        Re_lim2 = 3000;
        
    case 'PCHE32' % Semi-circular PCHE (32deg zigzag channels for Hoopes2016 comparisson)
        Re_lim1 = 1000;
        Re_lim2 = 2000;
        
    otherwise
        error('not implemented')
end

% Obtain array indices for laminar, transition and turbulent regimes
lam = find(Re<=Re_lim1);
tra = find(Re>Re_lim1 & Re<=Re_lim2);
tur = find(Re>Re_lim2);

% Compute the transition "weight" parameter, which is used to compute
% parameters in the transition regime. 0<W<1 (0=laminar,1=turbulent)
W = (Re(tra) - Re_lim1)/(Re_lim2 - Re_lim1);

% Allocate matrices
Cf = zeros(size(Re));
Nu = zeros(size(Re));

% Compute friction coefficient and Nusselt numbers. Use fully developed
% flow model with constant heat flux.
% Sources:
% --> Incropera&DeWitt for circular channels

% Obtain friction coefficient
Cf(lam) = friction_coefficient( Re(lam), 'lam', shape );
Cf_lim1 = friction_coefficient( Re_lim1, 'lam', shape );
Cf(tur) = friction_coefficient( Re(tur), 'tur', shape );
Cf_lim2 = friction_coefficient( Re_lim2, 'tur', shape );


% Nusselt number
Nu(lam) = Nusselt_number( Re(lam),      [],      [], 'lam', shape);
Nu_lim1 = Nusselt_number( Re_lim1,      [],      [], 'lam', shape);
Nu(tur) = Nusselt_number( Re(tur), Pr(tur), Cf(tur), 'tur', shape);
Nu_lim2 = Nusselt_number( Re_lim2, Pr(tra), Cf_lim2, 'tur', shape);

% Use Graetz number to account for effects due to entry region
entry_effects = 0;
switch entry_effects
    case 0 % No entry effects
        Nu = Nu + 0.0.*Gz;
        
    case 1 % Consider entry effects
        %{
        Nu(lam) = 3.66 + 0.0668*Gz(lam)./(1 + 0.04.*Gz(lam).^(2/3));
        Gz_lim1 = Gz(tra).*Re_lim1./Re(tra);
        Nu_lim1 = 3.66 + 0.0668*Gz_lim1./(1 + 0.04.*Gz_lim1.^(2/3));
        %}
        
        %%{
        % Thermally developing entry
        % Data from Hong and Bergles as displayed by Rohsenow
        TAB_xt = [
            0.000458, 0.000954, 0.00149, 0.00208, 0.00271,...
            0.00375,  0.00493,  0.00627, 0.00777, 0.00946,...
            0.0128,   0.0168,   0.0217,  0.0279,  0.0351,...
            0.0442,   0.0552,   0.0686,  0.0849,  0.105,...
            0.130,    0.159,    0.500,   1.00,   10.0 ];
        
        TAB_Nu = [
            17.71, 13.72, 11.80, 10.55, 9.605,...
            8.475, 7.723, 7.137, 6.556, 6.300,...
            5.821, 5.396, 5.077, 4.767, 4.562,...
            4.429, 4.276, 4.217, 4.156, 4.124,...
            4.118, 4.108, 4.089, 4.089, 4.089];
        
        % Create fit using interpolation
        x = TAB_xt'; y = TAB_Nu';
        Nu(lam) = 10.^interp1(log10(x),log10(y),log10(1./Gz(lam)),'makima');
        Gz_lim1 = Gz(tra).*Re_lim1./Re(tra);
        Nu_lim1 = 10.^interp1(log10(x),log10(y),log10(1./Gz_lim1),'makima');
        Nu_lim1(isnan(Nu_lim1)) = 4.089;
        Nu(isnan(Nu)) = 4.089;
        if any(isnan(Nu))
            keyboard
        end
        %}
        
        %{
        % Hydrodynamically developing entry
        xp = Pr./Gz;
        Cf(lam) = ( 13.74.*xp(lam).^(1/2) + (1.25 + 64.*xp(lam) - 13.74.*xp(lam).^(1/2))./(1 + 0.00021.*xp(lam).^(-2)) )./(4.*xp(lam).*Re(lam));
        Cf(isnan(Cf)) = 16./Re(isnan(Cf));
        if any(isnan(Cf))
            keyboard
        end
        %}
end

% Compute transition points
Cf(tra) = (1-W).*Cf_lim1 + W.*Cf_lim2;
Nu(tra) = (1-W).*Nu_lim1 + W.*Nu_lim2;

%Nu = 3.5*Nu;
%{
if strcmp(shape,'PCHE')
    Cf = Cf*1.1;
    Nu = Nu*1.2;
    %Cf = Cf*1;
    %Nu = Nu*1;
end
%}

% Compute Stanton number and heat transfer coefficient
St = Nu./(Re.*Pr);
ht = G*Cp.*St; %ht = Nu*k/D;

end

function [ Cf ] = friction_coefficient ( Re , regime, shape )
% Cf is the friction coefficient, also known as the "fanning friction
% factor". It is defined here as: dp = - 2*Cf*v*G^2*dx/D

switch shape
    case 'circular'
        
        switch regime
            case 'lam'
                aL = 16;
                Cf = aL./Re;
                
            case 'tur'
                Cf = (1/4)*(0.790 * log(Re) - 1.64).^(-2);
                
            otherwise
                error("incorrect regime input, please set 'lam' or 'tur'")
        end
        
    case 'PCHE'
        
        switch regime
            case 'lam'
                %Figley2013
                aL = 15.78;
                Cf = aL./Re;
                
            case 'tur'
                %Figley2013
                Cf = (1/4)*(1./ (1.8*log10(Re) - 1.5)).^2;
                
                %Ngo
                %Cf = 0.1942*Re.^(-0.091)
                %keyboard
                
            %case 'tur' %Marchionni2019
            %    A  = -2.0.*log10(12./Re);
            %    B  = -2.0.*log10(2.51.*A./Re);
            %    Cf = (1/4)*( 4.781 - (A - 4.781).^2./(B - 2*A + 4.781) ).^(-2);
                
            otherwise
                error("incorrect regime input, please set 'lam' or 'tur'")
        end
        
        case 'PCHE32'
        
        switch regime
            case 'lam'
                %Figley2013
                aL = 15.78;
                Cf = aL./Re;
                
            case 'tur'
                %Kim2016
                Cf = 0.2515*Re.^(-0.2031);
                
            otherwise
                error("incorrect regime input, please set 'lam' or 'tur'")
        end
        
    otherwise
        error('not implemented')
end

end

function [ Nu ] = Nusselt_number ( Re , Pr, Cf, regime, shape )

switch regime
    case 'lam'
        
        switch shape
            case 'circular'
                Nu_lam = 4.36;
                
            case 'PCHE' %Figley2013
                Nu_lam = 4.089;
                
            case 'PCHE32'
                Nu_lam = 4.089;
                
            otherwise
                error('not implemented')
        end
        Nu = Nu_lam + 0.*Re;
        
    case 'tur'
        switch shape
            case 'circular'
                Nu = (Cf/2).*(Re - 1000).*Pr./...
                    (1 + 12.7*(Cf/2).^(1/2).*(Pr.^(2/3)-1));
                
            case 'PCHE'
                %j  = 0.125*Re.^(-0.36);
                %Nu = j.*Pr.^(1/3).*Re;
                Nu = (Cf/2).*(Re - 1000).*Pr./...
                    (1 + 12.7*(Cf/2).^(1/2).*(Pr.^(2/3)-1));
                %Nu = 3.5*Nu;
                %Nu = 0.1696.*Re.^(0.629).*Pr.^(0.317)
                %keyboard
                %Nu = 0.023.*Re.^(0.8).*Pr.^(0.33);
                
            case 'PCHE32' %Kim2016
                Nu = 0.0292*Re.^0.8138;
                
            otherwise
                error('not implemented')
        end
        
    otherwise
        error("incorrect regime input, please set 'lam' or 'tur'")
end

end