% Set stage indices
iG = 1;  % keeps track of the gas stage number
iH = ones(1,Nhot);  % keeps track of the Hot fluid stream number
iC = ones(1,Ncld); % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the Environment (heat rejection) stream number
iA = 1 ;
iPMP = 1 ;

% Compute PR_dis based on charge pressure ratio and PRr
PRdis = PRr*PRch;

% Initial guess of discharge conditions
% Expander outlet (regenerator hot inlet)
gas.state(iL,1).p    = pbot; gas.state(iL,1).T = gas.state(1,1).T ;%fluidC(Ncld).state(1,1).T+1. ;
gas.state(iL,1).mdot = Load.mdot(iL);
[gas] = update(gas,[iL,1],1);

% Regenerator inlet indeces
if Nrcp == 0
    iT1 = 1 + Nc_ch*2 + 1 + Ne_ch*2 ;
    switch Load.mode
        case 4 % sCO2-PTES
            iT3 = 1 + 3*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + cooling + compression)    
        case 5 % Heat engine only
            iT3 = 1 + 2*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + compression)    
            error('Not implemented')
    end
    gas.state(iL,iT3).T = T0;
    gas.state(iL,iT3).p = gas.state(iL,1).p*PRdis;
    gas.state(iL,iT3).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iT3],1);
elseif Nrcp == 1
    % Regenerator cold inlet
    iReg1 = 1; % index regenerator hot inlet
    switch Load.mode
        case 4 % sCO2-PTES
            iReg2 = iReg1 + 1 + 3*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + cooling + compression)    
        case 5 % Heat engine only
            iReg2 = iReg1 + 1 + 2*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + compression)    
            error('Not implemented')
    end
    gas.state(iL,iReg2).T = T0;
    gas.state(iL,iReg2).p = gas.state(iL,1).p*PRdis;
    gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iReg2],1);
      
elseif Nrcp == 2
    % Regenerator cold inlet
    iReg1 = 1; % index regenerator hot inlet
    switch Load.mode
        case 4 % sCO2-PTES
            iReg2 = iReg1 + 3 + 3*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + cooling + compression)    
            iReg3 = iReg2 - 1 ;
            if Lrcmp
                iReg2 = iReg2 + 1 ;
            end
        case 5 % Heat engine only
            error('Not implemented')
            iReg2 = iReg1 + 2 + 3*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + compression)    
            iReg3 = iReg2 - 1 ;
            if Lrcmp
                iReg2 = iReg2 + 1 ;
            end
    end
    gas.state(iL,iReg2).T = TthreshD;
    gas.state(iL,iReg2).p = gas.state(iL,1).p*PRdis;
    gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iReg2],1);
    
    gas.state(iL,iReg3).T = T0;
    gas.state(iL,iReg3).p = gas.state(iL,1).p*PRdis;
    gas.state(iL,iReg3).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iReg3],1);
    
else
    error('Not implemented')
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
A_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
max_iter=100;
for counter=1:max_iter
    fprintf(1,['Discharging sCO2 PTES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    % REGENERATE (gas-gas)
    if Nrcp == 1
        [HX(Nhot+Ncld+1),gas,iG,~,~] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg2,0,0);
    elseif Nrcp == 2
        [HX(Nhot+Ncld+1),gas,iG,~,~] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg2,0,0);
        if Lrcmp
            [HX(Nhot+Ncld+2),gas,iG,~,~] = hex_func(HX(Nhot+Ncld+2),iL,gas,iG,gas,iReg3,3,min(TthreshD,gas.state(iL,iG).T-5)); % Low-temp regenerator
        else
            [HX(Nhot+Ncld+2),gas,iG,~,~] = hex_func(HX(Nhot+Ncld+2),iL,gas,iG,gas,iReg3,0,0); % Low-temp regenerator
        end
    end
    
    % RECOMPRESSION
    if Lrcmp
        % Set mass flows
        iRCMP = iReg2 + 4 ;
        gas.state(iL,iRCMP) = gas.state(iL,iG); 
        gas.state(iL,iRCMP).mdot = gas.state(iL,iG).mdot - gas.state(iL,iReg3).mdot ; % Mass flow into recompressor
        gas.state(iL,iG).mdot = gas.state(iL,iReg3).mdot ; % Set mass flow that goes on through heat rejection and cold store
        
        switch Load.mode
        case 4 % PTES
            ind = 1.0 ;
        case 2 % Heat engine only
            error('Not implemented')
            %ind = 0.0 ;
        end
        %p_aim = gas.state(iL,iRCMP).p*PRdis*(1-ploss)^ind;
        p_aim = gas.state(iL,iReg3).p ;
        [RCMP,gas,~] = compexp_func (RCMP,iL,gas,iRCMP,'Paim',p_aim,design_mode) ; 
        gas.stage(iL,iRCMP+1).type = 'comp'; % This seems to be necessary to get recompressor written and plotted
        gas.stage(iL,iRCMP+2).type = 'comp';% This seems to be necessary to get recompressor written and plotted
    end
    
    PRc_dis = (PRdis)^(1/Nc_dis)/(1-ploss)^2; % expansion pressure ratio
    for iN = 1:Nc_dis
        % REJECT HEAT (external HEX) *IF* cold store is below ambient (plus a bit)
        if fluidC(1).state(1,1).T <= T0+10
            T_aim = environ.T0 + Trej;    
            %[gas,environ,iG,iE] = hex_set(gas,[iL,iG],environ,[iL,iE],T_aim,eff,ploss);
            air.state(iL,iA).T = T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
            if T_aim >= gas.state(iL,iG).T 
                gas.state(iL,iG+1) = gas.state(iL,iG) ;
                iG = iG + 1 ;
            else
                [HX(end-1), gas, iG, air, iA] = hex_func(HX(end-1),iL,gas,iG,air,iA,5,T_aim);
                [DFAN(1),air,iA] = compexp_func (DFAN(1),iL,air,iA,'Paim',p0,1);
            end
        end    

        switch Load.mode
            case 4 % sCO2-PTES
                % COOL (gas-liquid)
                for ii = Ncld : -1 : 1
                    fluidC(ii).state(iL,iC(ii)).T = CT(ii).B(iL).T; fluidC(ii).state(iL,iC(ii)).p = CT(ii).B(iL).p;
                    TCoutMIN = min(gas.state(iL,iG).T-1,CT(ii).A(1).T) ;
                    [fluidC(ii)] = update(fluidC(ii),[iL,iC(ii)],1);
                    %[HX(Nhot+ii),gas,iG,fluidC(ii),iC(ii)] = hex_func(HX(Nhot+ii),iL,gas,iG,fluidC(ii),iC(ii),3,TCoutMIN);
                    %[HX(Nhot+ii),gas,iG,fluidC(ii),iC(ii)] = hex_func(HX(Nhot+ii),iL,gas,iG,fluidC(ii),iC(ii),1,0.75);
                    [HX(Nhot+ii),gas,iG,fluidC(ii),iC(ii)] = hex_func(HX(Nhot+ii),iL,gas,iG,fluidC(ii),iC(ii),1,1);
                    [DPMP(iPMP),fluidC(ii),iC(ii)] = compexp_func (DPMP(iPMP),iL,fluidC(ii),iC(ii),'Paim',fluidC(ii).state(iL,1).p,1);
                    iC(ii) = iC(ii) + 1 ; iPMP=iPMP+1;
                end
                
                % If cold store outlet is above ambient then reject heat
                if gas.state(iL,iG).T > T0
                    T_aim = environ.T0 + 5.;    
                    %[gas,environ,iG,iE] = hex_set(gas,[iL,iG],environ,[iL,iE],T_aim,eff,ploss);
                    air.state(iL,iA).T = T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
                    [HX(end), gas, iG, air, iA] = hex_func(HX(end),iL,gas,iG,air,iA,5,T_aim);
                    [DFAN(2),air,iA] = compexp_func (DFAN(2),iL,air,iA,'Paim',p0,1);
                end 
                
            case 5 % Heat engine only
        end
                
        % COMPRESS
        p_aim = gas.state(iL,iG).p*PRc_dis;
        [DCMP(iN),gas,iG] = compexp_func (DCMP(iN),iL,gas,iG,'Paim',p_aim,design_mode) ; 
    end
    
    % REGENERATE (gas-gas)
    if Nrcp == 1
        [HX(Nhot+Ncld+1),~,~,gas,iG] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg2,0,0);
    elseif Nrcp == 2
        if Lrcmp
            % If there is a recompression, adjust mass flows and add a mixer between the recomp. outlet and LTR cold outlet
            [HX(Nhot+Ncld+2),~,~,gas,iG] = hex_func(HX(Nhot+Ncld+2),iL,gas,iReg1+1,gas,iReg3,3,min(TthreshD,gas.state(iL,iReg1+1).T-5)); % Require cold side to reach a certain temp
            gas.state(iL,iRCMP).mdot   = gas.state(iL,iReg1).mdot - gas.state(iL,iReg3).mdot ;
            gas.state(iL,iRCMP+1).mdot = gas.state(iL,iRCMP).mdot ;
            [MIX(1),gas,~,~] = mix_streams(MIX(1),gas,[iL,iG],[iL,iRCMP+1]) ;
        else
            [HX(Nhot+Ncld+2),~,~,gas,~] = hex_func(HX(Nhot+Ncld+2),iL,gas,iReg1+1,gas,iReg3,0,0);
        end
        [HX(Nhot+Ncld+1),~,~,gas,iG] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg2,0,0);
    end
    
    for iN = 1:Ne_dis
        % HEAT (gas-fluid)
        for ii = Nhot : -1 : 1
            fluidH(ii).state(iL,iH(ii)).T = HT(ii).B(iL).T; fluidH(ii).state(iL,iH(ii)).p = HT(ii).B(iL).p; 
            THoutMAX = max(gas.state(iL,iG).T+1,HT(ii).A(1).T);
            [fluidH(ii)] = update(fluidH(ii),[iL,iH(ii)],1);
            [HX(ii),fluidH(ii),iH(ii),gas,iG] = hex_func(HX(ii),iL,fluidH(ii),iH(ii),gas,iG,4,THoutMAX); % Mode 4.
            [DPMP(iPMP),fluidH(ii),iH(ii)] = compexp_func (DPMP(iPMP),iL,fluidH(ii),iH(ii),'Paim',fluidH(ii).state(iL,1).p,1);
            iH(ii) = iH(ii) + 1 ; iPMP=iPMP+1;
        end
        
        % EXPAND
        PRe_dis = (gas.state(iL,iG).p/pbot)^(1/(Ne_dis+1-iN));  % expansion pressure ratio
        p_aim = gas.state(iL,iG).p/PRe_dis;
        [DEXP(iN),gas,iG] = compexp_func (DEXP(iN),iL,gas,iG,'Paim',p_aim,design_mode) ; 
    end
    
    % Determine convergence and proceed
    A = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    convergence = all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) ;
    if (convergence && strcmp(HX_model_temp,HX_model)) || counter==max_iter % is discharge cycle converged?
        % Close working fluid streams
        gas.stage(iL,iG).type = 'end';
        gas = count_Nstg(gas);
        
        % Close air (heat rejection) streams
        iA_out = 1:3:(iA-1); iA_in  = iA_out + 2;
        for i=iA_in, air.stage(iL,i).type = 'end'; end
        air = count_Nstg(air);
        
        % Close storage fluid streams
        for ii = 1 : Nhot
            iH_out = 1:3:(iH(ii)-1); iH_in  = iH_out + 2;
            for i=iH_in, fluidH(ii).stage(iL,i).type = 'end'; end
            fluidH(ii) = count_Nstg(fluidH(ii));
        end
        for ii = 1 : Ncld
            iC_out = 1:3:(iC(ii)-1); iC_in  = iC_out + 2;
            for i=iC_in, fluidC(ii).stage(iL,i).type = 'end'; end
            fluidC(ii) = count_Nstg(fluidC(ii));
        end
        
        % Uncomment these lines to print states
        %print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        %for ii=1:Nhot, print_states(fluidH(ii),iL,1:fluidH(ii).Nstg(iL)+1,Load); end
        %for ii=1:Ncld, print_states(fluidC(ii),iL,1:fluidC(ii).Nstg(iL)+1,Load); end
        
        % Exit loop
        break
    elseif convergence
        % If convergence has been reach but HX_model_temp~=HX_model, set
        % the heat exchanger models back to the original using the new
        % converged state, and resume iteration
        HX_model_temp = HX_model;
        for ihx=1:length(HX)
            HX(ihx).model = HX_model_temp;
            if strcmp(HX(ihx).name,'rej')
                HX(ihx).model = 'eff';
            end
        end
        gas.state(iL,1) = gas.state(iL,iG);
        A_0 = A;
        iG=1; iH=1; iC=1; iA=1; iPMP=1; iE=1;

    else
        % Set new initial conditions
        gas.state(iL,1) = gas.state(iL,iG);
        A_0 = A;
        iG=1; iH=ones(1,Nhot); iC=ones(1,Ncld); iE=1; iA=1; iPMP=1;
    end
end
if counter==max_iter
    warning('Exiting sCO2_DISCHARGE cycle without having reached convergence');
end


% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of (at least) one of the source tanks to stablish
% discharge time.

% Find t_dis (minimum for both sides to avoid depletion)
t_dis  = HT(1).B(iL).M/fluidH(1).state(iL,1).mdot;
for ii = 1:Nhot
    iH_out = 1:3:(iH(ii)-1);
    MdotH = total_mdot(fluidH(ii),iL,iH_out);
    t_dis  = min([t_dis,HT(ii).B(iL).M/MdotH]);
end
if Load.mode == 4
    for ii = 1:Ncld
        iC_out = 1:3:(iC(ii)-1);
        MdotC = total_mdot(fluidC(ii),iL,iC_out);
        t_dis  = min([t_dis,CT(ii).B(iL).M/MdotC]);
    end
elseif Load.mode == 5
end
Load.time(iL) = min([Load.time(iL),t_dis])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
for ii = 1 : Nhot
    iH_out = 1:3:(iH(ii)-1); iH_in  = iH_out + 2;
    [HT(ii)] = run_tanks(HT(ii),iL,fluidH(ii),iH_out,iH_in,Load,T0);
end
if Load.mode==4
    % Cold tanks
    for ii = 1 : Ncld
        iC_out = 1:3:(iC(ii)-1); iC_in  = iC_out + 2;
        [CT(ii)] = run_tanks(CT(ii),iL,fluidC(ii),iC_out,iC_in,Load,T0);
    end
end
% Atmospheric tanks
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);