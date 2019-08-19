% Set stage indices
i  = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the Environment (heat rejection) stream number

% Compute PR_dis based on charge pressure ratio and PRr
PRdis = PRr*PRch;

% Initial guess of discharge conditions
% Expander outlet (regenerator hot inlet)
gas.state(2,1).p    = pbot; gas.state(2,1).T = T1;
gas.state(2,1).mdot = mdot;
[gas] = update(gas,[2,1],1);
% Regenerator cold inlet
iReg1 = 1; % index regenerator hot inlet
switch mode
    case 0 % PTES
        iReg2 = iReg1 + 1 + 3*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + cooling + compression)    
    case 2 % Heat engine only
        iReg2 = iReg1 + 1 + 2*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + compression)    
end
gas.state(2,iReg2).T = T0;
gas.state(2,iReg2).p = gas.state(2,1).p*PRdis;
gas.state(2,iReg2).mdot = gas.state(2,1).mdot;
[gas] = update(gas,[2,iReg2],1);

% Set matrix of temperature and pressure points to test convergence
A_0 = [[gas.state(2,:).T];[gas.state(2,:).p]];
while 1
    %fprintf(1,'Hello discharge PTES\n')
    
    % REGENERATE (gas-gas)    
    [gas,~,i,~] = hex_TQ_cond(gas,[2,iReg1],gas,[2,iReg2],eff,0,ploss,'regen',0,0); %#ok<*SAGROW>
    
    PRc_dis = (PRdis)^(1/Nc_dis)/(1-ploss)^2; % expansion pressure ratio
    for iN = 1:Nc_dis
        % REJECT HEAT (external HEX)
        T_aim = environ.T0;
        [gas,environ,i,iE] = hex_set(gas,[2,i],environ,[2,iE],T_aim,eff,ploss);
        
        switch mode
            case 0 % PTES
                % COOL (gas-liquid)
                fluidC(iC).state(2,1).T = CT.B(3).T; fluidC(iC).state(2,1).p = CT.B(3).p;% fluidC(iC).state(2,1).mdot = fluidC(Ne_ch+1-iC).state(1,1).mdot;
                [fluidC(iC)] = update(fluidC(iC),[2,1],1);
                [gas,fluidC(iC),i,~] = hex_TQ_cond(gas,[2,i],fluidC(iC),[2,1],eff,1.0,ploss,'hex',0,0);
                iC=iC+1;
            case 1 % Heat engine only
        end
        
        % COMPRESS
        p_aim = gas.state(2,i).p*PRc_dis;
        [gas,i] = compexp(gas,[2,i],eta,p_aim,3); %#ok<*SAGROW>
        
        %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(2,i0).T,gas.state(2,i0).p/1e5,gas.stage(2,i0).type,i0); end; fprintf(1,'\n');
    end
    
    % REGENERATE (gas-gas)
    [~,gas,~,i] = hex_TQ_cond(gas,[2,iReg1],gas,[2,iReg2],eff,0,ploss,'regen',0,0); %#ok<*SAGROW>
    
    for iN = 1:Ne_dis
        % HEAT (gas-fluid)
        fluidH(iH).state(2,1).T = HT.B(3).T; fluidH(iH).state(2,1).p = HT.B(3).p; THmin = HT.A(1).T;
        [fluidH(iH)] = update(fluidH(iH),[2,1],1);
        [fluidH(iH),gas,~,i] = hex_TQ_cond(fluidH(iH),[2,1],gas,[2,i],eff,1.0,ploss,'hex',2, THmin);
        iH=iH+1;
        
        % EXPAND
        PRe_dis = (gas.state(2,i).p/pbot)^(1/(Ne_dis+1-iN));  % expansion pressure ratio
        p_aim = gas.state(2,i).p/PRe_dis;
        [gas,i] = compexp(gas,[2,i],eta,p_aim,1);
    end
    
    % Close cycle
    gas.stage(2,i).type = gas.stage(2,1).type;
    
    % Determine convergence and proceed
    A = [[gas.state(2,:).T];[gas.state(2,:).p]];
    %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(2,i0).T,gas.state(2,i0).p/1e5,gas.stage(2,i0).type,i0); end; fprintf(1,'\n');
    %disp((A(A~=0) - A_0(A~=0))./A(A~=0)*100);
    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is discharge cycle converged?
        % Exit discharge cycle
        %T1d = gas.state(2,i).T;
        %pbotd = gas.state(2,i).p;
        gas_min_rho_dis = gas.state(2,i); %take data for power density calculation
        %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(2,i0).T,gas.state(2,i0).p/1e5,gas.stage(2,i0).type,i0); end; fprintf(1,'\n');
        break
    else
        % Set new initial conditions
        gas.state(2,1) = gas.state(2,i);
        A_0 = A;
        i=1;iH=1;iHc=1;iC=1;iE=1;
    end
end


% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of (at least) one of the source tanks to stablish
% discharge time.

% Find t_dis (minimum for both cycles to avoid depletion)
[MdotH] = total_Mdot(fluidH,[2,1]);
t_dis  = HT.B(3).M/MdotH;
if mode == 0
    [MdotC] = total_Mdot(fluidC,[2,1]);
    tC_dis  = CT.B(3).M/MdotC;
    t_dis   = min([t_dis,tC_dis]);
end

% During discharge, B is source and A is sink
[HT.B(4),HT.A(3),HT.A(4),WL_mixH_dis] = liquid_tanks_compute(fluidH,2,2,HT.B(3),t_dis,T0);
if mode == 0
    [CT.B(4),CT.A(3),CT.A(4),WL_mixC_dis] = liquid_tanks_compute(fluidC,2,2,CT.B(3),t_dis,T0);
end