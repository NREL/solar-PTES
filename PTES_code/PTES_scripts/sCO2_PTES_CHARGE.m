% Set stage indices
i  = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the heat rejection stream number

% Initial guess of charge conditions
% Compressor inlet (regenerator hot outlet)
gas.state(iL,1).p    = pbot; gas.state(iL,1).T = T1;
gas.state(iL,1).mdot = Load.mdot(iL);
[gas] = update(gas,[iL,1],1);
% Regenerator inlet indeces
if Nrcp == 0
    iT1 = Ncld + Nc_ch*2 + Nhot + Ne_ch*2 ;
    gas.state(iL,iT1).T = T0; % Assume T1 is at ambient conditions
    gas.state(iL,iT1).p = pbot;
    gas.state(iL,iT1).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iT1],1);
elseif Nrcp == 1
    iReg1 = Ncld + Nc_ch*2;         % index hot regenerator inlet (after compression + cooling)
    iReg2 = iReg1 + 1 + Nhot + Ne_ch*2; % index cold regenerator inlet (after regeneration + heat_reject + expansion + heating)
    gas.state(iL,iReg2).T = T0;
    gas.state(iL,iReg2).p = pbot;
    gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iReg2],1);
else
    error('Not implemented')
end

% Set matrix of temperature and pressure points to test convergence
A_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
while 1
    
    fprintf(1,'Hello charge PTES\n')
    
    for iN = 1:Nc_ch
        % COMPRESS
        if setTmax
            T_aim = Tmax;
            [gas,i] = compexp(gas,[iL,i],eta,T_aim,0);
        else
            PRc = (pmax/gas.state(iL,i).p)^(1/(Nc_ch+1-iN));
            %PRc = (PRch)^(1/Nc_ch)/(1-ploss) % pressure ratio for each compression stage
            p_aim = gas.state(iL,i).p*PRc;
            [gas,i] = compexp(gas,[iL,i],eta,p_aim,3);
        end
        ptop  = gas.state(iL,i).p;
        
        % COOL (gas-liquid)
        fluidH(iH).state(iL,1).T = HT.A(iL).T; fluidH(iH).state(iL,1).p = HT.A(iL).p; %#ok<*SAGROW>
        [fluidH(iH)] = update(fluidH(iH),[iL,1],1);
        [gas,fluidH(iH),i,~] = hex_TQ_cond(gas,[iL,i],fluidH(iH),[iL,1],eff,1.0,ploss,'hex',0,0);
        plot_hex(gas,[iL,i-1],fluidH(iH),[iL,1],100,10); % Hot storage
        
        if Nhot == 2
            fluidH2(iH).state(iL,1).T = HT2.A(iL).T; fluidH2(iH).state(iL,1).p = HT2.A(iL).p;
            [fluidH2(iH)] = update(fluidH2(iH),[iL,1],1);
            [gas,fluidH2(iH),i,~] = hex_TQ_cond(gas,[iL,i],fluidH2(iH),[iL,1],eff,1.0,ploss,'hex',0,0);
            plot_hex(gas,[iL,i-1],fluidH2(iH),[iL,1],100,11); % Cold storage
        end
                
        iH=iH+1;
    end    
    if setTmax, PRch = ptop/pbot; end
    
    % REGENERATE (gas-gas)
    if Nrcp == 1
        iIN = i;
        [gas,~,i,~] = hex_TQ_cond(gas,[iL,iReg1],gas,[iL,iReg2],eff,0,ploss,'regen',0,0);
    end
    
    % May wish to make cold store as cold as possible, or avoid rejecting
    % heat here to avoid the worst of the CO2 c_p variation
    if Lcld
        T_aim = environ.T0 ;
    else
        T_aim = gas.state(iL,i).T - 1.0 ;
    end
    %     % REJECT HEAT (external HEX)
    [gas,environ,i,iE] = hex_set(gas,[iL,i],environ,[iL,iE],T_aim,eff,ploss);
    
    for iN = 1:Ne_ch
        % EXPAND
        PRe = (gas.state(iL,i).p/pbot)^(1/(Ne_ch+1-iN)); %expansion pressure ratio
        p_aim = gas.state(iL,i).p/PRe;
        [gas,i] = compexp(gas,[iL,i],eta,p_aim,1);        
        
        % HEAT (gas-liquid)
        fluidC(iC).state(iL,1).T = CT.A(iL).T; fluidC(iC).state(iL,1).p = CT.A(iL).p;
        [fluidC(iC)] = update(fluidC(iC),[iL,1],1);
        [fluidC(iC),gas,~,i] = hex_TQ_cond(fluidC(iC),[iL,1],gas,[iL,i],eff,1.6,ploss,'hex', 0, 0);
        plot_hex(gas,[iL,i-1],fluidC(iC),[iL,1],100,12); % Cold storage
        
        if Ncld == 2
            fluidC2(iC).state(iL,1).T = CT2.A(iL).T; fluidC2(iC).state(iL,1).p = CT2.A(iL).p;
            [fluidC2(iC)] = update(fluidC2(iC),[iL,1],1);
            [fluidC2(iC),gas,~,i] = hex_TQ_cond(fluidC2(iC),[iL,1],gas,[iL,i],eff,1.0,ploss,'hex', 0, 0);
            plot_hex(gas,[iL,i-1],fluidC2(iC),[iL,1],100,13); % Cold storage
        end
        
        iC=iC+1;
    end
    
    % REGENERATE (gas-gas)
    if Nrcp == 1
        [~,gas,~,i] = hex_TQ_cond(gas,[iL,iReg1],gas,[iL,iReg2],eff,0,ploss,'regen',0,0); 
        plot_hex(gas,[iL,i-1],gas,[iL,iIN],100,14); % Recuperator
    end
    
    % Close cycle
    gas.stage(iL,i).type = gas.stage(iL,1).type;
    gas = count_Nstg(gas);
    
    % Determine convergence and proceed
    A = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(iL,i0).T,gas.state(iL,i0).p/1e5,gas.stage(iL,i0).type,i0); end; fprintf(1,'\n');
    %disp((A(A~=0) - A_0(A~=0))./A(A~=0)*100);
    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is charge cycle converged?
        % Exit charge cycle
        %T1 = gas.state(iL,1).T;
        %pbot = gas.state(iL,1).p;
        gas_min_rho_ch = gas.state(iL,1); %take data for power density calculation
        %for i0=1:i, fprintf(1,'\n %f\t%f\t%10s\t%d',gas.state(iL,i0).T,gas.state(iL,i0).p/1e5,gas.stage(iL,i0).type,i0); end; fprintf(1,'\n');
        break
    else
        % Set new initial conditions
        gas.state(iL,1) = gas.state(iL,i);        
        A_0 = A;
        i=1; iH=1; iC=1; iE=1;
    end
end

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
[HT] = run_tanks(HT,fluidH,iL,Load,T0);
if Nhot == 2; [HT2] = run_tanks(HT2,fluidH2,iL,Load,T0);end
% Cold tanks
[CT] = run_tanks(CT,fluidC,iL,Load,T0);
if Ncld == 2; [CT2] = run_tanks(CT2,fluidC2,iL,Load,T0);end