% Set stage indices
i  = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the Environment (heat rejection) stream number

% Compute PR_dis based on charge pressure ratio and PRr
PRdis = PRr*PRch;

% Initial guess of discharge conditions
% Expander outlet (regenerator hot inlet)
gas.state(iL,1).p    = pbot; gas.state(iL,1).T = T1;
gas.state(iL,1).mdot = Load.mdot(iL);
[gas] = update(gas,[iL,1],1);
% Regenerator cold inlet
iReg1 = 1; % index regenerator hot inlet
switch mode
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
    [gas,~,i,~] = hex_TQ_cond(gas,[iL,iReg1],gas,[iL,iReg2],eff,0,ploss,'regen',0,0); %#ok<*SAGROW>
    
    PRc_dis = (PRdis)^(1/Nc_dis)/(1-ploss)^2; % expansion pressure ratio
    for iN = 1:Nc_dis
        % REJECT HEAT (external HEX)
        T_aim = environ.T0;
        [gas,environ,i,iE] = hex_set(gas,[iL,i],environ,[iL,iE],T_aim,eff,ploss);
        
        switch mode
            case 0 % PTES
                % COOL (gas-liquid)
                fluidC(iC).state(iL,1).T = CT.B(3).T; fluidC(iC).state(iL,1).p = CT.B(3).p;% fluidC(iC).state(iL,1).mdot = fluidC(Ne_ch+1-iC).state(1,1).mdot;
                [fluidC(iC)] = update(fluidC(iC),[iL,1],1);
                [gas,fluidC(iC),i,~] = hex_TQ_cond(gas,[iL,i],fluidC(iC),[iL,1],eff,1.0,ploss,'hex',0,0);
                iC=iC+1;
            case 1 % Heat engine only
        end
        
        % COMPRESS
        p_aim = gas.state(iL,i).p*PRc_dis;
        [gas,i] = compexp(gas,[iL,i],eta,p_aim,3); %#ok<*SAGROW>
        
        %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(iL,i0).T,gas.state(iL,i0).p/1e5,gas.stage(iL,i0).type,i0); end; fprintf(1,'\n');
    end
    
    % REGENERATE (gas-gas)
    [~,gas,~,i] = hex_TQ_cond(gas,[iL,iReg1],gas,[iL,iReg2],eff,0,ploss,'regen',0,0); %#ok<*SAGROW>
    
    for iN = 1:Ne_dis
        % HEAT (gas-fluid)
        fluidH(iH).state(iL,1).T = HT.B(3).T; fluidH(iH).state(iL,1).p = HT.B(3).p; THmin = HT.A(1).T;
        [fluidH(iH)] = update(fluidH(iH),[iL,1],1);
        [fluidH(iH),gas,~,i] = hex_TQ_cond(fluidH(iH),[iL,1],gas,[iL,i],eff,1.0,ploss,'hex',2, THmin);
        iH=iH+1;
        
        % EXPAND
        PRe_dis = (gas.state(iL,i).p/pbot)^(1/(Ne_dis+1-iN));  % expansion pressure ratio
        p_aim = gas.state(iL,i).p/PRe_dis;
        [gas,i] = compexp(gas,[iL,i],eta,p_aim,1);
    end
    
    % Close cycle
    gas.stage(iL,i).type = gas.stage(iL,1).type;
    
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
[MdotH] = total_Mdot(fluidH,[iL,1]);
t_dis  = HT.B(3).M/MdotH;
if mode == 0
    [MdotC] = total_Mdot(fluidC,[iL,1]);
    tC_dis  = CT.B(3).M/MdotC;
    t_dis   = min([t_dis,tC_dis]);
end
Load.time(iL) = min([Load.time(iL),t_dis])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
[HT] = run_tanks(HT,fluidH,iL,Load,T0);
if mode == 0
    % Cold tanks
    [CT] = run_tanks(CT,fluidC,iL,Load,T0);
end