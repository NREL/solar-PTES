% Set stage indices
iG = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the Environment (heat rejection) stream number
iA = 1;  % keeps track of the Air (heat rejection) stream number

% Compute PR_dis based on charge pressure ratio and PRr
PRdis = PRr*PRch;

% Initial guess of discharge conditions
% Expander outlet (regenerator hot inlet)
gas.state(iL,1).p    = pbot; gas.state(iL,1).T = T1;
gas.state(iL,1).mdot = Load.mdot(iL);
[gas] = update(gas,[iL,1],1);
% Regenerator cold inlet
iReg1 = 1; % index regenerator hot inlet
switch Load.mode
    case 0 % PTES
        iReg2 = iReg1 + 1 + 3*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + cooling + compression)    
    case 2 % Heat engine only
        iReg2 = iReg1 + 1 + 2*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + compression)    
end
gas.state(iL,iReg2).T = T0;
gas.state(iL,iReg2).p = gas.state(iL,1).p*PRdis;
gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
[gas] = update(gas,[iL,iReg2],1);

% Set matrix of temperature and pressure points to test convergence
A_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
while 1
    fprintf(1,'Hello discharge PTES\n')
    
    % REGENERATE (gas-gas)
    if new_hex_calls
        %[REGEN,gas,iG,~,~] = hex_func(REGEN,iL,gas,iReg1,gas,iReg2,0,0);
        [RCP,gas,iG,~,~] = hex_func(RCP,iL,gas,iReg1,gas,iReg2,0,0); % New call using hx_class
    else
        [gas,~,iG,~] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0);
    end
        
    for iN = 1:Nc_dis
        % REJECT HEAT (external HEX)
        % REPLACE THIS WITH hex_func call?
        T_aim = environ.T0;
        [gas,environ,iG,iE] = hex_set(gas,[iL,iG],environ,[iL,iE],T_aim,eff,ploss);
        
        switch Load.mode
            case 0 % PTES
                % COOL (gas-liquid)
                fluidC.state(iL,iC).T = CT.B(iL).T; fluidC.state(iL,iC).p = CT.B(iL).p; %#ok<*SAGROW>
                [fluidC] = update(fluidC,[iL,iC],1);
                if new_hex_calls
                    %[HX,gas,iG,fluidC,iC] = hex_func(HX,iL,gas,iG,fluidC,iC,1,1.0);
                    [HXc,gas,iG,fluidC,iC] = hex_func(HXc,iL,gas,iG,fluidC,iC,1,1.0);
                else
                    [gas,fluidC,iG,iC] = hex_TQ(gas,[iL,iG],fluidC,[iL,iC],eff,ploss,'hex',1,1.0);
                end
                iC=iC+1;
            case 1 % Heat engine only
        end
        
        % COMPRESS
        PRc_dis = (PRdis)^(1/Nc_dis)/(1-ploss)^2; % stage compression pressure ratio
        p_aim = gas.state(iL,iG).p*PRc_dis;
        [DCMP(iN),gas,iG] = compexp_func (DCMP(iN),iL,gas,iG,'Paim',p_aim) ; 
    end
    
    % REGENERATE (gas-gas)
    if new_hex_calls
        %[REGEN,~,~,gas,iG] = hex_func(REGEN,iL,gas,iReg1,gas,iReg2,0,0);
        [RCP,~,~,gas,iG] = hex_func(RCP,iL,gas,iReg1,gas,iReg2,0,0); % New call using hx_class
    else
        [~,gas,~,iG] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0);
    end
        
    for iN = 1:Ne_dis
        % HEAT (gas-fluid)
        fluidH.state(iL,iH).T = HT.B(iL).T; fluidH.state(iL,iH).p = HT.B(iL).p; THmin = HT.A(1).T;
        [fluidH] = update(fluidH,[iL,iH],1);
        Taim = THmin;
        if new_hex_calls
            %[HX,fluidH,iH,gas,iG] = hex_func(HX,iL,fluidH,iH,gas,iG,2,1.0);
            [HXh,fluidH,iH,gas,iG] = hex_func(HXh,iL,fluidH,iH,gas,iG,2,1.0); % New call using hx_class
        else
            [fluidH,gas,iH,iG] = hex_TQ(fluidH,[iL,iH],gas,[iL,iG],eff,ploss,'hex',2,1.0);
        end            
        iH=iH+1;
        
        % EXPAND
        PRe_dis = (gas.state(iL,iG).p/pbot)^(1/(Ne_dis+1-iN));  % stage expansion pressure ratio
        p_aim = gas.state(iL,iG).p/PRe_dis;
        [DEXP(iN),gas,iG] = compexp_func (DEXP(iN),iL,gas,iG,'Paim',p_aim) ; 
    end
    
    % Determine convergence and proceed
    A = [[gas.state(iL,:).T];[gas.state(iL,:).p]];

    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is discharge cycle converged?
        % Close working fluid streams
        gas.stage(iL,iG).type = 'end';
        gas = count_Nstg(gas);
        
        % Close air (heat rejection) streams
        iA_out = 1:2:(iA-1); iA_in  = iA_out + 1;
        for i=iA_in, air.stage(iL,i).type = 'end'; end
        air = count_Nstg(air);
        
        % Close storage fluid streams
        iH_out = 1:2:(iH-1); iH_in  = iH_out + 1;
        iC_out = 1:2:(iC-1); iC_in  = iC_out + 1;
        for i=iH_in, fluidH.stage(iL,i).type = 'end'; end
        for i=iC_in, fluidC.stage(iL,i).type = 'end'; end
        fluidH = count_Nstg(fluidH);
        fluidC = count_Nstg(fluidC);
        
        % Uncomment these lines to print states
        %print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        %print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
        %print_states(fluidC,iL,1:fluidC.Nstg(iL)+1,Load);
        
        % Exit loop
        break
    else
        % Set new initial conditions
        gas.state(iL,1) = gas.state(iL,iG);
        A_0 = A;
        iG=1;iH=1;iHc=1;iC=1;iE=1;
    end
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
iA_out = 0; iA_in = 0;
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);