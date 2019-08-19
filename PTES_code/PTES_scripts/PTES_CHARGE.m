% Set stage indices
i  = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the heat rejection stream number

% Initial guess of charge conditions
% Compressor inlet (regenerator hot outlet)
gas.state(1,1).p    = pbot; gas.state(1,1).T = T1;
gas.state(1,1).mdot = mdot;
[gas] = update(gas,[1,1],1);
% Regenerator inlet indeces
iReg1 = 1 + Nc_ch*2;         % index hot regenerator inlet (after compression + cooling)
iReg2 = iReg1 + 2 + Ne_ch*2; % index cold regenerator inlet (after regeneration + heat_reject + expansion + heating)
gas.state(1,iReg2).T = T0;
gas.state(1,iReg2).p = pbot;
gas.state(1,iReg2).mdot = gas.state(1,1).mdot;
[gas] = update(gas,[1,iReg2],1);

% Set matrix of temperature and pressure points to test convergence
A_0 = [[gas.state(1,:).T];[gas.state(1,:).p]];
while 1
    
    %fprintf(1,'Hello charge PTES\n')
    
    for iN = 1:Nc_ch
        % COMPRESS
        if setTmax
            T_aim = Tmax;
            [gas,i] = compexp(gas,[1,i],eta,T_aim,0);
        else
            PRc = (pmax/gas.state(1,i).p)^(1/(Nc_ch+1-iN));
            %PRc = (PRch)^(1/Nc_ch)/(1-ploss) % pressure ratio for each compression stage
            p_aim = gas.state(1,i).p*PRc;
            [gas,i] = compexp(gas,[1,i],eta,p_aim,3);
        end
        ptop  = gas.state(1,i).p;
        
        % COOL (gas-liquid)
        fluidH(iH).state(1,1).T = HT.A(1).T; fluidH(iH).state(1,1).p = HT.A(1).p; %#ok<*SAGROW>
        [fluidH(iH)] = update(fluidH(iH),[1,1],1);
        [gas,fluidH(iH),i,~] = hex_TQ_cond(gas,[1,i],fluidH(iH),[1,1],eff,1.0,ploss,'hex',0,0);
        iH=iH+1;
    end    
    if setTmax, PRch = ptop/pbot; end
    
    % REGENERATE (gas-gas)
    [gas,~,i,~] = hex_TQ_cond(gas,[1,iReg1],gas,[1,iReg2],eff,0,ploss,'regen',0,0);
    
    %     % REJECT HEAT (external HEX)
    T_aim = environ.T0;
    [gas,environ,i,iE] = hex_set(gas,[1,i],environ,[1,iE],T_aim,eff,ploss);
    
    for iN = 1:Ne_ch
        % EXPAND
        PRe = (gas.state(1,i).p/pbot)^(1/(Ne_ch+1-iN)); %expansion pressure ratio
        p_aim = gas.state(1,i).p/PRe;
        [gas,i] = compexp(gas,[1,i],eta,p_aim,1);        
        
        % HEAT (gas-liquid)
        fluidC(iC).state(1,1).T = CT.A(1).T; fluidC(iC).state(1,1).p = CT.A(1).p;
        [fluidC(iC)] = update(fluidC(iC),[1,1],1);
        [fluidC(iC),gas,~,i] = hex_TQ_cond(fluidC(iC),[1,1],gas,[1,i],eff,1.0,ploss,'hex', 0, 0);
        iC=iC+1;
    end
    
    % REGENERATE (gas-gas)
    [~,gas,~,i] = hex_TQ_cond(gas,[1,iReg1],gas,[1,iReg2],eff,0,ploss,'regen',0,0); 
    
    % Close cycle
    gas.stage(1,i).type = gas.stage(1,1).type;
    
    % Determine convergence and proceed
    A = [[gas.state(1,:).T];[gas.state(1,:).p]];
    %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(1,i0).T,gas.state(1,i0).p/1e5,gas.stage(1,i0).type,i0); end; fprintf(1,'\n');
    %disp((A(A~=0) - A_0(A~=0))./A(A~=0)*100);
    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is charge cycle converged?
        % Exit charge cycle
        %T1 = gas.state(1,1).T;
        %pbot = gas.state(1,1).p;
        gas_min_rho_ch = gas.state(1,1); %take data for power density calculation
        %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(1,i0).T,gas.state(1,i0).p/1e5,gas.stage(1,i0).type,i0); end; fprintf(1,'\n');
        break
    else
        % Set new initial conditions
        gas.state(1,1) = gas.state(1,i);        
        A_0 = A;
        i=1; iH=1; iC=1; iE=1;
    end
end

% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of source tanks to stablish initial mass
% Hot tanks
[HT.A(1)] = liquid_tank_start(fluidH,[1,1],HT.A(1),t_ch,T0);
[HT.A(2),HT.B(1),HT.B(2),WL_mixH_ch] = liquid_tanks_compute(fluidH,1,2,HT.A(1),t_ch,T0);

% Cold tanks
[CT.A(1)] = liquid_tank_start(fluidC,[1,1],CT.A(1),t_ch,T0);
[CT.A(2),CT.B(1),CT.B(2),WL_mixC_ch] = liquid_tanks_compute(fluidC,1,2,CT.A(1),t_ch,T0);