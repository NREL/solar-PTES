% Script that starts up the load cycle
if strcmp(Load(Iload),'c')
    fprintf('\nLOAD CYCLE #%i. CHARGING.\n\n',Iload);
elseif strcmp(Load(Iload),'d')
    fprintf('\nLOAD CYCLE #%i. DISCHARGING.\n\n',Iload);
elseif strcmp(Load(Iload),'s')
    fprintf('\nLOAD CYCLE #%i. STORAGE.\n\n',Iload);
end

Lend = false ; % Should the load cycle be ended?

% Set steam properties as constant along pipe
% Inlet properties of steam are known
if strcmp(Load(Iload),'c')
    
    % Find steam input conditions. Assume conditions are constant along
    % pipe for the starting point
    if Lreadload
        mdot(Iload) = (dsg_mdot(Iload) - dsg_mdot0) / Npipe ;
        xs    = dsg_q(Iload) * ones(Nx,2) ; % Dryness fraction
        Gs    = mdot(Iload)  * ones(Nx,2) / Ap ; % Mass flux
        Ts    = dsg_T(Iload) * ones(Nx,2) ; % Temperature
        Ps    = PsC   * ones(Nx,2) ; % Pressure
        hs    = RP1('PQ_INPUTS',PsC(1),xs(1),'H',steam) * ones(Nx,2) ;
        ss    = RP1('PQ_INPUTS',PsC(1),xs(1),'S',steam) * ones(Nx,2) ;
        rhos  = RP1('PQ_INPUTS',PsC(1),xs(1),'D',steam) * ones(Nx,2) ;
        mus   = RP1('PQ_INPUTS',PsC(1),xs(1),'VISCOSITY',steam) * ones(Nx,2) ;
        us    = Gs ./ rhos ;
    else
        xs    = xsC   * ones(Nx,2) ; % Dryness fraction
        Gs    = mdot(Iload) * ones(Nx,2) / Ap ; % Mass flux
        Ts    = TsC   * ones(Nx,2) ; % Temperature
        Ps    = PsC   * ones(Nx,2) ; % Pressure
        hs    = hsC   * ones(Nx,2) ; % Enthalpy
        ss    = ssC   * ones(Nx,2) ; % Entropy
        rhos  = rhosC * ones(Nx,2) ; % 's' stands for steam
        mus   = musC  * ones(Nx,2) ; % Viscosity
        us    = Gs ./ rhos ;
    end
    
    % Check whether fully charged
    % Refine this later
    if mean(xp(:,1)) > XPend_chg % if average melted fraction is > 90% consider it to be fully charged
        Lend = true ;
        fprintf('Charging cycle not started as storage is already charged.\n');
    end
    
    hl = hlC ;
    hv = hvC ;
    vl = vlC ;
    vv = vvC ;
    Tsat = TsC ;
    
    Tcoef = TcoefC ;
    Vcoef = VcoefC ;
elseif strcmp(Load(Iload),'d')
    if Lreadload
        mdot(Iload) = (dsg_mdot0 - dsg_mdot(Iload)) / Npipe ; % Probably need something more sophisticated than this based on power
        xs    = zeros(Nx,2) ; % Dryness fraction. Comes in as water, but aim to leave as dsg_q0
        Gs    = mdot(Iload)  * ones(Nx,2) / Ap ; % Mass flux
        Ts    = TsD  * ones(Nx,2) ; % Temperature
        Ps    = PsD   * ones(Nx,2) ; % Pressure
        hs    = RP1('PQ_INPUTS',PsD(1),xs(1),'H',steam) * ones(Nx,2) ;
        ss    = RP1('PQ_INPUTS',PsD(1),xs(1),'S',steam) * ones(Nx,2) ;
        rhos  = RP1('PQ_INPUTS',PsD(1),xs(1),'D',steam) * ones(Nx,2) ;
        mus   = RP1('PQ_INPUTS',PsD(1),xs(1),'VISCOSITY',steam) * ones(Nx,2) ;
        us    = Gs ./ rhos ;
    else
        rhos  = rhosD * ones(Nx,2) ; % 's' stands for steam
        Gs    = mdot(Iload) * ones(Nx,2) / Ap ; % Mass flux
        Ts    = TsD   * ones(Nx,2) ; % Temperature
        Ps    = PsD   * ones(Nx,2) ; % Pressure
        hs    = hsD   * ones(Nx,2) ; % Enthalpy
        ss    = ssD   * ones(Nx,2) ; % Entropy
        xs    = xsD   * ones(Nx,2) ; % Dryness fraction
        mus   = musD  * ones(Nx,2) ; % Viscosity
        us    = Gs ./ rhos ;
    end
    
    
    % Check whether fully discharged
    % Refine this later
    if mean(xp(:,1)) < XPend_dis % if average melted fraction is < 10% consider it to be fully discharged
        Lend = true ;
        fprintf('Discharging cycle not started as storage is already discharged.\n');
    end
         
    hl = hlD ;
    hv = hvD ;
    vl = vlD ;
    vv = vvD ;
    Tsat = TsD ;
    
    Tcoef = TcoefD ;
    Vcoef = VcoefD ;
else
    error('You`ve mucked up')
end
