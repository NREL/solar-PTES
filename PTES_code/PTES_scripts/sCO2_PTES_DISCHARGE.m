% Set stage indices
i  = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the Environment (heat rejection) stream number

% Compute PR_dis based on charge pressure ratio and PRr
PRdis = PRr*PRch;

% Initial guess of discharge conditions
% Expander outlet (regenerator hot inlet)
gas.state(iL,1).p    = pbot; gas.state(iL,1).T = fluidC(Ncld).state(1,1).T+1. ;
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


% Set matrix of temperature and pressure points to test convergence
A_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
while 1
    fprintf(1,'Hello discharge sCO2-PTES\n')
    % REGENERATE (gas-gas)
    if Nrcp == 1
        [gas,~,i,~] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0); %#ok<*SAGROW>
    elseif Nrcp == 2
        [gas,~,i,~] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0);   % High-temp regenerator
        if Lrcmp
            [gas,~,i,~] = hex_TQ(gas,[iL,i],gas,[iL,iReg3],eff,ploss,'regen',3,min(TthreshD,gas.state(iL,i).T-5));% Low-temp regenerator
        else
            [gas,~,i,~] = hex_TQ(gas,[iL,i],gas,[iL,iReg3],eff,ploss,'regen',0,0);       % Low-temp regenerator
        end
    end
    
    % RECOMPRESSION
    if Lrcmp
        % Set mass flows
        iRCMP = iReg2 + 4 ;
        gas.state(iL,iRCMP) = gas.state(iL,i); 
        gas.state(iL,iRCMP).mdot = gas.state(iL,i).mdot - gas.state(iL,iReg3).mdot ; % Mass flow into recompressor
        gas.state(iL,i).mdot = gas.state(iL,iReg3).mdot ; % Set mass flow that goes on through heat rejection and cold store
        
        switch Load.mode
        case 4 % PTES
            ind = 1.0 ;
        case 2 % Heat engine only
            error('Not implemented')
            %ind = 0.0 ;
        end
        p_aim = gas.state(iL,iRCMP).p*PRdis*(1-ploss)^ind;
        %[gas,~] = compexp(gas,[iL,iRCMP],eta,p_aim,3);
        [RCMP,gas,~] = compexp_func (RCMP,gas,[iL,iRCMP],'Paim',p_aim) ; 
        gas.stage(iL,iRCMP+1).type = 'comp'; % This seems to be necessary to get recompressor written and plotted
        gas.stage(iL,iRCMP+2).type = 'comp';% This seems to be necessary to get recompressor written and plotted
    end
    
    PRc_dis = (PRdis)^(1/Nc_dis)/(1-ploss)^2; % expansion pressure ratio
    for iN = 1:Nc_dis
        % REJECT HEAT (external HEX) *IF* cold store is below ambient (plus a bit)
        if fluidC(1).state(1,1).T <= T0+10
            T_aim = environ.T0;    
            [gas,environ,i,iE] = hex_set(gas,[iL,i],environ,[iL,iE],T_aim,eff,ploss);
        end    

        switch Load.mode
            case 4 % sCO2-PTES
                % COOL (gas-liquid)
                for ii = Ncld : -1 : 1
                    fluidC(ii).state(iL,1).T = CT(ii).B(iL).T; fluidC(ii).state(iL,1).p = CT(ii).B(iL).p;
                    TCoutMIN = min(gas.state(iL,i).T-5,CT(ii).A(1).T) ;
                    [fluidC(ii)] = update(fluidC(ii),[iL,1],1);
                    [gas,fluidC(ii),i,~] = hex_TQ(gas,[iL,i],fluidC(ii),[iL,1],eff,ploss,'hex',3,TCoutMIN); %Mode 1. Hot gas mdot known, cold fluid mdot not known
                    iC = ii ;
                end
                iC = iC + 1 ;
                
                % If cold store outlet is above ambient then reject heat
                if gas.state(iL,i).T > T0
                    T_aim = environ.T0;    
                    [gas,environ,i,iE] = hex_set(gas,[iL,i],environ,[iL,iE],T_aim,eff,ploss);
                end 
                
            case 5 % Heat engine only
        end
                
        % COMPRESS
        p_aim = gas.state(iL,i).p*PRc_dis;
        %[gas,i] = compexp(gas,[iL,i],eta,p_aim,3); %#ok<*SAGROW>
        [DCMP(iN),gas,i] = compexp_func (DCMP(iN),gas,[iL,i],'Paim',p_aim) ; 
        
        %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(iL,i0).T,gas.state(iL,i0).p/1e5,gas.stage(iL,i0).type,i0); end; fprintf(1,'\n');
    end
    
    % REGENERATE (gas-gas)
    if Nrcp == 1
        [~,gas,~,i] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0); %#ok<*SAGROW>
    elseif Nrcp == 2
        if Lrcmp
            % If there is a recompression, adjust mass flows and add a mixer between the recomp. outlet and LTR cold outlet
            [~,gas,~,i] = hex_TQ(gas,[iL,iReg1+1],gas,[iL,iReg3],eff,ploss,'regen',3,min(TthreshD,gas.state(iL,iReg1+1).T-5)); %% Require cold side to reach a certain temp
            gas.state(iL,iRCMP).mdot   = gas.state(iL,iReg1).mdot - gas.state(iL,iReg3).mdot ;
            gas.state(iL,iRCMP+1).mdot = gas.state(iL,iRCMP).mdot ;
            [gas,i,~] = mix_streams(gas,[iL,i],[iL,iRCMP+1]) ;
        else
            [~,gas,~,~] = hex_TQ(gas,[iL,iReg1+1],gas,[iL,iReg3],eff,ploss,'regen',0,0); 
        end
        [~,gas,~,i] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0); 
    end
    
    
    for iN = 1:Ne_dis
        % HEAT (gas-fluid)
        for ii = Nhot : -1 : 1
            fluidH(ii).state(iL,1).T = HT(ii).B(iL).T; fluidH(ii).state(iL,1).p = HT(ii).B(iL).p; 
            THoutMAX = max(gas.state(iL,i).T+1,HT(ii).A(1).T);
            [fluidH(ii)] = update(fluidH(ii),[iL,1],1);
            [fluidH(ii),gas,~,i] = hex_TQ(fluidH(ii),[iL,1],gas,[iL,i],eff,ploss,'hex',4, THoutMAX); % Mode 4.       
            iH = ii ;
        end
        iH = iH + 1 ;
               
        % EXPAND
        PRe_dis = (gas.state(iL,i).p/pbot)^(1/(Ne_dis+1-iN));  % expansion pressure ratio
        p_aim = gas.state(iL,i).p/PRe_dis;
        %[gas,i] = compexp(gas,[iL,i],eta,p_aim,1);
        [DEXP(iN),gas,i] = compexp_func (DEXP(iN),gas,[iL,i],'Paim',p_aim) ; 
    end
    
    % Close cycle
    gas.stage(iL,i).type = gas.stage(iL,1).type;
    gas = count_Nstg(gas);
    
    % Determine convergence and proceed
    A = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(iL,i0).T,gas.state(iL,i0).p/1e5,gas.stage(iL,i0).type,i0); end; fprintf(1,'\n');
    %disp((A(A~=0) - A_0(A~=0))./A(A~=0)*100);
    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is discharge cycle converged?
        % Exit discharge cycle
        %T1d = gas.state(iL,i).T;
        %pbotd = gas.state(iL,i).p;
        gas_min_rho_dis = gas.state(iL,i); %take data for power density calculation
        %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(iL,i0).T,gas.state(iL,i0).p/1e5,gas.stage(iL,i0).type,i0); end; fprintf(1,'\n');
        break
    else
        % Set new initial conditions
        gas.state(iL,1) = gas.state(iL,i);
        A_0 = A;
        i=1;iH=1;iHc=1;iC=1;iE=1;
    end
end


% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of (at least) one of the source tanks to stablish
% discharge time.

% Find t_dis (minimum for both cycles to avoid depletion)
% Think this needs to be modified to account for numerous tanks

if Nhot > 1 || Ncld > 1
    % Hot tanks
    tH_dis = 1.e11;
    for ii = 1:Nhot
       tH_dis = min(tH_dis, HT(ii).B(iL).M / fluidH(ii).state(iL,1).mdot);
    end
    
    % Cold tanks
    tC_dis = 1.e11;
    for ii = 1:Ncld
       tC_dis = min(tC_dis, CT(ii).B(iL).M / fluidC(ii).state(iL,1).mdot);
    end
    t_dis   = min([tH_dis,tC_dis]);
    Load.time(iL) = min([Load.time(iL),t_dis])*(1-1e-6);
else
    % This was the original method, which works for multiple
    % compressions/expansions
    [MdotH] = total_Mdot(fluidH,[iL,1]);
    t_dis  = HT(1).B(iL).M/MdotH;
    if Load.mode == 4
        [MdotC] = total_Mdot(fluidC,[iL,1]);
        tC_dis  = CT(1).B(3).M/MdotC;
        t_dis   = min([t_dis,tC_dis]);
    elseif Load.mode == 5
        t_dis = Load.time(iL) ;
    end
    Load.time(iL) = min([Load.time(iL),t_dis])*(1-1e-6);
end
% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
for ii = 1 : Nhot
    [HT(ii)] = run_tanks(HT(ii),fluidH(ii),iL,Load,T0);
end

if Load.mode == 4
    % Cold tanks
    for ii = 1 : Ncld
        [CT(ii)] = run_tanks(CT(ii),fluidC(ii),iL,Load,T0);
    end
end