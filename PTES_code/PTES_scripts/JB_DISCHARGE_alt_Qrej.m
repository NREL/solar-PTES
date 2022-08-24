% In this version of JB_DISCHARGE heat is rejected in an alternative
% location - specifically, heat is rejected after the compressor. This
% means that heat is rejected in the same point of the cycle as the charge
% cycle and that the same heat rejection equipment can be used for charge
% and discharge.


% Set stage indices
iG = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the Environment (heat rejection) stream number
iA = 1;  % keeps track of the Air (heat rejection) stream number
iPMP = 1 ; % Keeps track of which pump is being used

% Regenerator cold inlet
iReg1 = 1; % index regenerator hot inlet
switch Load.mode
    case 0 % PTES
        iReg2 = iReg1 + 1 + 3*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + cooling + compression)    
    case 2 % Heat engine only
        iReg2 = iReg1 + 1 + 2*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + compression)    
end

%fprintf('Heat exchanger summary\n');
%print_hexs(HX,1,'Charge:\n');
%print_hexs(HX,2,'Discharge:\n');
%keyboard

% Compute PR_dis based on charge pressure ratio and PRr
PRdis = PRr*PRch;

if design_mode == 1
    % Initial guess of discharge conditions
    % Expander outlet (regenerator hot inlet)
    gas.state(iL,1).p    = pbot; gas.state(iL,1).T = T1;
    gas.state(iL,1).mdot = Load.mdot(iL);
    [gas] = update(gas,[iL,1],1);
    gas.state(iL,iReg2).T = T0;
    gas.state(iL,iReg2).p = gas.state(iL,1).p*PRdis;
    gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iReg2],1);
    TOLconv = 1e-3 ;
    environ.T0 = T0 ;
else
    for ii = 1 : numel(D(D~=0))/2
        gas.state(iL,ii).T    = D(1,ii) ;
        gas.state(iL,ii).p    = D(2,ii) ;
        gas.state(iL,ii).mdot = Load.mdot(iL) ;
        
        % For inventory control, assume that the pressure scales with the off-design mass flow rate
        %gas.state(iL,ii).p = (DEXP.Pin/DEXP.pr0) * Load.mdot(iL) / DEXP.mdot0 ;
        gas.state(iL,ii).p = gas0.state(2,ii).p * (Load.mdot(iL) / DEXP.mdot0) *  sqrt(T0_off(iL) / T0) ; % Second ever run is discharging
        [gas] = update(gas,[iL,ii],1);
        
    end   
    p1guess = gas.state(iL,1).p ;
    switcher = 0;

    p1prev   = 0 ;
    erprevP   = 0 ;
    gradPP     = 0 ;
    gradPT = 0;

    T1prev   = 0 ;
    erprevT   = 0 ;
    gradTP     = 0 ;
    gradTT = 0;
    environ.T0 = T0_off(iL) ;
    TOLconv = 1e-4 ;
end

switch HX_model
    case 'eff'
        HX_model_temp = 'eff';
    case 'geom'
        % Set the heat exchanger models to 'eff' temporarily, to obtain
        % approximated cycle points before using the 'geom' model 
        HX_model_temp = 'eff';
        for ihx=1:length(HX)
            HX(ihx).model = HX_model_temp;
        end
    otherwise
        error('not implemented')
end

% Set matrix of temperature and pressure points to test convergence
D_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
max_iter = 150;
for counter = 1:max_iter
    fprintf(1,['Discharging JB PTES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    
    % REGENERATE (gas-gas)
    [HX(ihx_reg),gas,iG,~,~] = hex_func(HX(ihx_reg),iL,gas,iReg1,gas,iReg2,0,0);
    
    for iN = 1:Nc_dis

        if Nc_dis>1
            warning('Have not tested Nc_dis > 1.');
        end
        switch Load.mode
            case 0 % PTES
                % COOL (gas-liquid)
                fluidC.state(iL,iC).T = CT.B(iL).T; fluidC.state(iL,iC).p = CT.B(iL).p; %#ok<*SAGROW>
                [fluidC] = update(fluidC,[iL,iC],1);

                %CT.B(iL).T = 205;%fluidC.state(iL,iC).T ;
                %CT.B(iL).h = fluidC.state(iL,iC).h ;
                %CT.B(iL).H = CT.B(iL).h * CT.B(iL).M ;
                %CT         = update_tank_state(CT,CT.B(iL),T0,2);

                [HX(ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ihx_cld(iN)),iL,gas,iG,fluidC,iC,1,1.0);
                %[HX(ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ihx_cld(iN)),iL,gas,iG,fluidC,iC,3,fluidC.state(1,1).T);
                [DPMP(iPMP),fluidC,iC] = compexp_func (DPMP(iPMP),iL,fluidC,iC,'Paim',fluidC.state(iL,1).p,1);
                iC=iC+1; iPMP=iPMP+1;
            case 1 % Heat engine only
        end
        
        % COMPRESS
        PRc_dis = (PRdis)^(1/Nc_dis)/(1-ploss)^2; % stage compression pressure ratio
        p_aim = gas.state(iL,iG).p*PRc_dis;
        [DCMP(iN),gas,iG] = compexp_func (DCMP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ; 

        % REJECT HEAT (external HEX)
        air.state(iL,iA).T = environ.T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
        
        if design_mode == 0 && T0_off(iL) < T0
            % If off-design and T0 is colder than the design value
            % Reduce air mass flow rate so that expander inlet temperature
            % remains at the design point.
            [HX(ihx_rejd(iN)), gas, iG, air, iA] = hex_func(HX(ihx_rejd(iN)),iL,gas,iG,air,iA,5,gas0.state(iL,iG+1).T);
        else
            [HX(ihx_rejd(iN)), gas, iG, air, iA] = hex_func(HX(ihx_rejd(iN)),iL,gas,iG,air,iA,1,0.5);
        end
        [DFAN(1),air,iA] = compexp_func (DFAN(1),iL,air,iA,'Paim',p0,1);
    end
    
    % REGENERATE (gas-gas)
    [HX(ihx_reg),~,~,gas,iG] = hex_func(HX(ihx_reg),iL,gas,iReg1,gas,iReg2,0,0);
    
    for iN = 1:Ne_dis
        % HEAT (gas-fluid)
        fluidH.state(iL,iH).T = HT.B(iL).T; fluidH.state(iL,iH).p = HT.B(iL).p; THmin = HT.A(1).T;
        [fluidH] = update(fluidH,[iL,iH],1);
        Taim = THmin;

        [HX(ihx_hot(iN)),fluidH,iH,gas,iG] = hex_func(HX(ihx_hot(iN)),iL,fluidH,iH,gas,iG,4,THmin);
        [DPMP(iPMP),fluidH,iH] = compexp_func (DPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
        iH=iH+1; iPMP=iPMP+1;
        
        % EXPAND
        PRe_dis = (gas.state(iL,iG).p/pbot)^(1/(Ne_dis+1-iN));  % stage expansion pressure ratio
        p_aim = gas.state(iL,iG).p/PRe_dis;
        [DEXP(iN),gas,iG] = compexp_func (DEXP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ; 
    end
    
    % Determine convergence and proceed
    D = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    convergence = all(abs((D(D~=0) - D_0(D~=0))./D(D~=0))*100 < TOLconv);
    
    if (convergence && strcmp(HX_model_temp,HX_model)) || counter >= max_iter % is discharge cycle converged?
        % Close working fluid streams
        gas.stage(iL,iG).type = 'end';
        gas = count_Nstg(gas);
        
        % Close air (heat rejection) streams
        iA_out = 1:3:(iA-1); iA_in  = iA_out + 2;
        for i=iA_in, air.stage(iL,i).type = 'end'; end
        air = count_Nstg(air);
        
        % Close storage fluid streams
        iH_out = 1:3:(iH-1); iH_in  = iH_out + 2;
        iC_out = 1:3:(iC-1); iC_in  = iC_out + 2;
        for i=iH_in, fluidH.stage(iL,i).type = 'end'; end
        for i=iC_in, fluidC.stage(iL,i).type = 'end'; end
        fluidH = count_Nstg(fluidH);
        fluidC = count_Nstg(fluidC);
        
        % Uncomment these lines to print states
        
        print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
        print_states(fluidC,iL,1:fluidC.Nstg(iL)+1,Load);
        print_states(air,iL,1:air.Nstg(iL)+1,Load);
        %}
        
        % Exit loop
        break
        
    elseif convergence
        % If convergence has been reach but HX_model_temp~=HX_model, set
        % the heat exchanger models back to the original using the new
        % converged state, and resume iteration
        HX_model_temp = HX_model;
        for ihx=1:length(HX)
            HX(ihx).model = HX_model_temp;
        end
        gas.state(iL,1) = gas.state(iL,iG);        
        D_0 = D;
        iG=1;iH=1;iHc=1;iC=1;iE=1;iA=1;iPMP=1;
        
    else
        % Set new initial conditions
        if ~design_mode
            
            %print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
%              fprintf('p1:   %13.8f \n',gas.state(iL,1).p/1e5)
%              fprintf('pEND: %13.8f \n\n',gas.state(iL,iG).p/1e5)
%              fprintf('T1:   %13.8f \n',gas.state(iL,1).T)
%              fprintf('TEND: %13.8f \n\n',gas.state(iL,iG).T)

%{
            % Adjust inlet pressure to try to reach convergence. The
            % 'smoothing' factor has to be quite small (<0.1, say) for this to be stable
            smooth = 0.025;
            if switcher == 0 
                gas.state(iL,1).p = gas.state(iL,1).p - smooth * (gas.state(iL,iG).p - gas.state(iL,1).p) ;
                switcher = 1;
            elseif switcher == 1
            %gas.state(iL,1).p = gas.state(iL,1).p / (gas.state(iL,iReg2+2).p * DEXP.mdot0 / (DEXP.Pin * gas.state(iL,iReg2+2).mdot)) ;
                gas.state(iL,1).T = gas.state(iL,1).T + smooth * (gas.state(iL,iG).T - gas.state(iL,1).T) ;
                switcher = 0;
            elseif switcher == 2
                gas.state(iL,1).mdot = gas.state(iL,1).mdot + smooth * (Load.mdot(iL) /( gas.state(iL,1).p / p1guess * sqrt(gas0.state(iL,1).T / gas.state(iL,1).T)) - gas.state(iL,1).mdot);
                switcher = 0;
            end
            
            %gas.state(iL,1).mdot = Load.mdot(iL);
            %gas.state(iL,1).mdot = Load.mdot(iL) * gas.state(iL,1).p / p1guess * sqrt(gas0.state(iL,1).T / gas.state(iL,1).T);
%}

            if counter == 1
                
                p1prev = gas.state(iL,1).p ;
                T1prev = gas.state(iL,1).T ;
                erprevP = gas.state(iL,1).p  - gas.state(iL,iG).p ;
                erprevT = gas.state(iL,1).T  - gas.state(iL,iG).T ;
                smooth = 0.025 ;
                
                gas.state(iL,1).p = gas.state(iL,1).p - smooth * (gas.state(iL,iG).p - gas.state(iL,1).p) ;
                gas.state(iL,1).T = gas.state(iL,1).T + smooth * (gas.state(iL,iG).T - gas.state(iL,1).T) ;
            else
                ernewP = gas.state(iL,1).p  - gas.state(iL,iG).p ;
                gradPP  = (ernewP - erprevP) / (gas.state(iL,1).p - p1prev) ;
                gradPT  = (ernewP - erprevP) / (gas.state(iL,1).T - T1prev) ;
                p1prev = gas.state(iL,1).p ;
                
                ernewT = gas.state(iL,1).T  - gas.state(iL,iG).T ;
                gradTP  = (ernewT - erprevT) / (gas.state(iL,1).p - p1prev) ;
                gradTT  = (ernewT - erprevT) / (gas.state(iL,1).T - T1prev) ;
                T1prev = gas.state(iL,1).T ;

                erprevP = ernewP ;
                erprevT = ernewT ;

                gas.state(iL,1).p = gas.state(iL,1).p - 0.2 * ernewP / gradPP;% - 0.2 * ernewP / gradPT;
                gas.state(iL,1).T = gas.state(iL,1).T - 0.2 * ernewT / gradTT;% - 0.2 * ernewT / gradTP;
            end

            [gas] = update(gas,[iL,1],1);
            
        else
            gas.state(iL,1) = gas.state(iL,iG);
        end
        
        D_0 = D;
        iG=1;iH=1;iHc=1;iC=1;iE=1;iA=1;iPMP=1;
    end
end
if counter==max_iter
    warning('Exiting JB_DISCHARGE cycle without having reached convergence');
end

% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of (at least) one of the source tanks to stablish
% discharge time.

% Find t_dis
MdotH = total_mdot(fluidH,iL,iH_out);
t_dis  = HT.B(iL).M/MdotH;
if Load.mode == 0
    MdotC = total_mdot(fluidC,iL,iC_out);
    tC_dis  = CT.B(iL).M/MdotC;
    t_dis   = min([t_dis,tC_dis]);
end
Load.time(iL) = min([Load.time(iL),t_dis])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
HT = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
if Load.mode == 0
    % Cold tanks
    CT = run_tanks(CT,iL,fluidC,iC_out,iC_in,Load,T0);
end
% Atmospheric tanks
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);