% This script runs the discharging phase of a "Time-Shifted" Recompression
% sCO2 Power cycle. This power cycle does not use a recompressor. Instead
% the required heat comes from a 'low-temperature' hot storage system that
% was charged earlier.

% Set stage indices
iG = 1;  % keeps track of the gas stage number
iH = ones(1,Nhot);  % keeps track of the Hot fluid stream number
iC = ones(1,Ncld); % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the Environment (heat rejection) stream number

% Compute PR_dis based on charge pressure ratio and PRr
PRdis = PRr*PRch;

% Initial guess of discharge conditions
% Expander outlet (regenerator hot inlet)
gas.state(iL,1).p    = pbot; gas.state(iL,1).T = HT(1).B(1).T ;
gas.state(iL,1).mdot = Load.mdot(iL);
[gas] = update(gas,[iL,1],1);

% Regenerator inlet indeces
iReg1 = 1; % index HOT regenerator hot inlet
iReg2 = iReg1 + 4 + 3*Nc_dis; % index HOT regenerator cold inlet (after regeneration + heat rejection + compression)
iReg3 = iReg2 - 2 ;

gas.state(iL,iReg2).T = TthreshD;
gas.state(iL,iReg2).p = gas.state(iL,1).p*PRdis;
gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
[gas] = update(gas,[iL,iReg2],1);

gas.state(iL,iReg3).T = T0;
gas.state(iL,iReg3).p = gas.state(iL,1).p*PRdis;
gas.state(iL,iReg3).mdot = gas.state(iL,1).mdot;
[gas] = update(gas,[iL,iReg3],1);

% Set matrix of temperature and pressure points to test convergence
A_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
while 1
    fprintf(1,'Hello discharge TSRC-sCO2 PTES cycle\n')
    % REGENERATE (gas-gas)
    if new_hex_calls
        [HX(Nhot+Ncld+1),gas,iG,~,~] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg2,0,0);
        [HX(Nhot+Ncld+2),gas,iG,~,~] = hex_func(HX(Nhot+Ncld+2),iL,gas,iG,gas,iReg3,3,min(TthreshD,gas.state(iL,iG).T-5)); % Low-temp regenerator
    else
        [gas,~,iG,~] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0); %#ok<*SAGROW>
        [gas,~,iG,~] = hex_TQ(gas,[iL,iG],gas,[iL,iReg3],eff,ploss,'regen',3,min(TthreshD,gas.state(iL,iG).T-5));% Low-temp regenerator
    end
    
    % Split flows. Part goes through cold store. Part goes through heat rejection.
    % Set mass flows
    iCLD = iReg2 + 4 ; % Cold store is tabulated after the rest of the cycle
    gas.state(iL,iCLD)      = gas.state(iL,iG);
    gas.state(iL,iCLD).mdot = gas.state(iL,iG).mdot - gas.state(iL,iReg3).mdot ; % Mass flow into recompressor
    gas.state(iL,iG).mdot   = gas.state(iL,iReg3).mdot ; % Set mass flow that goes on through heat rejection and cold store
    
    % Check that isn't zero
    if gas.state(iL,iCLD).mdot == 0
        gas.state(iL,iCLD).mdot = 0.1 * gas.state(iL,iG).mdot ;
    end
    
    fluidC.state(iL,iC).T = CT.B(iL).T; fluidC.state(iL,iC).p = CT.B(iL).p;
    TCoutMIN = min(gas.state(iL,iG).T-5,CT.A(1).T) ;
    [fluidC] = update(fluidC,[iL,iC],1);
    if new_hex_calls
        [HX(Nhot+1),gas,~,fluidC,iC] = hex_func(HX(Nhot+1),iL,gas,iCLD,fluidC,iC,3,TCoutMIN);
    else
        [gas,fluidC,~,iC] = hex_TQ(gas,[iL,iCLD],fluidC,[iL,iC],eff,ploss,'hex',3,TCoutMIN);
    end
    iC = iC + 1 ;
                
    % Reject heat to environment (external HEX)
    T_aim = environ.T0 + Trej;
    [gas,environ,iG,iE] = hex_set(gas,[iL,iG],environ,[iL,iE],T_aim,eff,ploss);
    
    % Mixer
    [gas,iG,~] = mix_streams(gas,[iL,iG],[iL,iCLD+1]) ;
    
    PRc_dis = (PRdis)^(1/Nc_dis)/(1-ploss)^1; % expansion pressure ratio
    
    % COMPRESS
    p_aim = gas.state(iL,iG).p*PRc_dis;
    [DCMP(1),gas,iG] = compexp_func (DCMP(1),iL,gas,iG,'Paim',p_aim) ;
    
    % REGENERATE (gas-gas)
    if new_hex_calls
        
        % Split flow. Part goes through low-temp recuperator. Part goes
        % through low-temp storage.
        
        % Mass flow rates get adjusted
        % Low-temp recuperator:
        [HX(Nhot+Ncld+2),~,~,gas,iG] = hex_func(HX(Nhot+Ncld+2),iL,gas,iReg1+1,gas,iReg3,3,min(TthreshD,gas.state(iL,iReg1+1).T-5)); % Require cold side to reach a certain temp
        gas.state(iL,iCLD).mdot      = gas.state(iL,iReg1).mdot - gas.state(iL,iReg3).mdot ;
        gas.state(iL,iCLD+1).mdot    = gas.state(iL,iCLD).mdot ;
        
        gas.state(iL,iCLD+2)      = gas.state(iL,iReg3) ;
        gas.state(iL,iCLD+2).mdot = gas.state(iL,iCLD).mdot ;
        gas.state(iL,iCLD+3).mdot = gas.state(iL,iCLD).mdot ;
        
        % Now low-temp storage
        fluidH(2).state(iL,iH(2)).T = HT(2).B(iL).T; fluidH(2).state(iL,iH(2)).p = HT(2).B(iL).p; THmin = HT(2).A(1).T;
        [fluidH(2)] = update(fluidH(2),[iL,iH(2)],1);
        [HX(Nhot),fluidH(2),iH(2),gas,~] = hex_func(HX(Nhot),iL,fluidH(2),iH(2),gas,iCLD+2,4,THmin); % Require cold side to reach a certain temp
        iH(2) = iH(2) + 1; 
        
        [gas,~,~] = mix_streams(gas,[iL,iG],[iL,iCLD+3]) ;
        
        % High-temp recuperator
        [HX(Nhot+Ncld+1),~,~,gas,iG] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg2,0,0);
               
    else
        error('Not implemented')
        % If there is a recompression, adjust mass flows and add a mixer between the recomp. outlet and LTR cold outlet
%         [~,gas,~,iG] = hex_TQ(gas,[iL,iReg1+1],gas,[iL,iReg3],eff,ploss,'regen',3,min(TthreshD,gas.state(iL,iReg1+1).T-5)); %% Require cold side to reach a certain temp
%         gas.state(iL,iRCMP).mdot   = gas.state(iL,iReg1).mdot - gas.state(iL,iReg3).mdot ;
%         gas.state(iL,iRCMP+1).mdot = gas.state(iL,iRCMP).mdot ;
%         [gas,~,~] = mix_streams(gas,[iL,iG],[iL,iRCMP+1]) ;
%         
%         [~,gas,~,iG] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0);
               
    end
    
    for iN = 1:Ne_dis
        % HEAT (gas-fluid)
        fluidH(1).state(iL,iH(1)).T = HT(1).B(iL).T; fluidH(1).state(iL,iH(1)).p = HT(1).B(iL).p;
        THoutMAX = max(gas.state(iL,iG).T+1,HT(1).A(1).T);
        [fluidH(1)] = update(fluidH(1),[iL,iH(1)],1);
        if new_hex_calls
            [HX(1),fluidH(1),iH(1),gas,iG] = hex_func(HX(1),iL,fluidH(1),iH(1),gas,iG,4,THoutMAX); % Mode 4.
        else
            [fluidH(1),gas,iH(1),iG] = hex_TQ(fluidH(1),[iL,iH(1)],gas,[iL,iG],eff,ploss,'hex',4, THoutMAX); % Mode 4.
        end
        iH(1) = iH(1) + 1 ;
        
        % EXPAND
        PRe_dis = (gas.state(iL,iG).p/pbot)^(1/(Ne_dis+1-iN));  % expansion pressure ratio
        p_aim = gas.state(iL,iG).p/PRe_dis;
        [DEXP(iN),gas,iG] = compexp_func (DEXP(iN),iL,gas,iG,'Paim',p_aim) ; 
    end
    
    % Determine convergence and proceed
    A = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is discharge cycle converged?
        % Close working fluid streams
        gas.stage(iL,iG).type = 'end';
        gas = count_Nstg(gas);
        
        % Close storage fluid streams
        for ii = 1 : Nhot
            iH_out = 1:2:(iH(ii)-1); iH_in  = iH_out + 1;
            for i=iH_in, fluidH(ii).stage(iL,i).type = 'end'; end
            fluidH(ii) = count_Nstg(fluidH(ii));
        end
        for ii = 1 : Ncld
            iC_out = 1:2:(iC(ii)-1); iC_in  = iC_out + 1;
            for i=iC_in, fluidC(ii).stage(iL,i).type = 'end'; end
            fluidC(ii) = count_Nstg(fluidC(ii));
        end
        
        % Uncomment these lines to print states
        %print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        %for ii=1:Nhot, print_states(fluidH(ii),iL,1:fluidH(ii).Nstg(iL)+1,Load); end
        %for ii=1:Ncld, print_states(fluidC(ii),iL,1:fluidC(ii).Nstg(iL)+1,Load); end
        
        % Exit loop
        break
    else
        % Set new initial conditions
        gas.state(iL,1) = gas.state(iL,iG);
        A_0 = A;
        iG=1; iH=ones(1,Nhot); iC=ones(1,Ncld); iE=1;
    end
end


% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of (at least) one of the source tanks to stablish
% discharge time.

% Find t_dis (minimum for both sides to avoid depletion)
t_dis  = HT(1).B(iL).M/fluidH(1).state(iL,1).mdot;
for ii = 1:Nhot
    iH_out = 1:2:(iH(ii)-1);
    MdotH = total_mdot(fluidH(ii),iL,iH_out);
    t_dis  = min([t_dis,HT(ii).B(iL).M/MdotH]);
end
for ii = 1:Ncld
    iC_out = 1:2:(iC(ii)-1);
    MdotC = total_mdot(fluidC(ii),iL,iC_out);
    t_dis  = min([t_dis,CT(ii).B(iL).M/MdotC]);
end
Load.time(iL) = min([Load.time(iL),t_dis])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
for ii = 1 : Nhot
    iH_out = 1:2:(iH(ii)-1); iH_in  = iH_out + 1;
    [HT(ii)] = run_tanks(HT(ii),iL,fluidH(ii),iH_out,iH_in,Load,T0);
end

% Cold tanks
for ii = 1 : Ncld
    iC_out = 1:2:(iC(ii)-1); iC_in  = iC_out + 1;
    [CT(ii)] = run_tanks(CT(ii),iL,fluidC(ii),iC_out,iC_in,Load,T0);
end
