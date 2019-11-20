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
    iReg1 = Nhot + Nc_ch*2;         % index hot regenerator inlet (after compression + cooling)
    iReg2 = iReg1 + 1 + Ncld + Ne_ch*2; % index cold regenerator inlet (after regeneration + heat_reject + expansion + heating)
    gas.state(iL,iReg2).T = T0;
    gas.state(iL,iReg2).p = pbot;
    gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iReg2],1);
elseif Nrcp == 2
    iReg1 = Nhot + Nc_ch*2;             % index hot regenerator inlet (after compression + cooling)
    iReg2 = iReg1 + 2 + Ncld + Ne_ch*2; % index cold regenerator inlet (after 1xregeneration + heat_reject + expansion + heating)
    iReg3 = iReg2 + 1 ;                 % index cold regenerator inlet (after 2xregeneration + heat_reject + expansion + heating)
    gas.state(iL,iReg2).T = T0;         % Update cold inlet to cold temp regenerator
    gas.state(iL,iReg2).p = pbot;
    gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iReg2],1);
    
    gas.state(iL,iReg3).T = TthreshC;         % Update cold inlet to hot temp regenerator
    gas.state(iL,iReg3).p = pbot;
    gas.state(iL,iReg3).mdot = gas.state(iL,3).mdot;
    [gas] = update(gas,[iL,iReg3],1);
else
    error('Not implemented')
end

% Set matrix of temperature and pressure points to test convergence
A_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
while 1
    
    fprintf(1,'Hello charge sCO2-PTES\n')
    
    for iN = 1:Nc_ch
        % COMPRESS
        if setTmax
            T_aim = Tmax;
            %[gas,i] = compexp(gas,[iL,i],eta,T_aim,0);
            [CCMP(iN),gas,i] = compexp_func (CCMP(iN),gas,[iL,i],'Taim',T_aim) ; 
        else
            PRc = (pmax/gas.state(iL,i).p)^(1/(Nc_ch+1-iN));
            %PRc = (PRch)^(1/Nc_ch)/(1-ploss) % pressure ratio for each compression stage
            p_aim = gas.state(iL,i).p*PRc;
            %[gas,i] = compexp(gas,[iL,i],eta,p_aim,3);
            [CCMP(iN),gas,i] = compexp_func (CCMP(iN),gas,[iL,i],'Paim',p_aim) ; 
        end
        ptop  = gas.state(iL,i).p;
        
        % COOL (gas-liquid)
        for ii = 1 : Nhot
            fluidH(ii).state(iL,1).T = HT(ii).A(iL).T; fluidH(ii).state(iL,1).p = HT(ii).A(iL).p; %#ok<*SAGROW>
            [fluidH(ii)] = update(fluidH(ii),[iL,1],1);
            if ii == 1
                %[gas,fluidH(ii),i,~] = hex_TQ(gas,[iL,i],fluidH(ii),[iL,1],eff,ploss,'hex',1,1); % Mode 1: don't know hot fluid mdot. Crat = 1
                [gas,fluidH(ii),i,zzz,HX1] = hex_TQ(gas,[iL,i],fluidH(ii),[iL,1],eff,ploss,'hex',1,1);
            else
                TCoutMIN = fluidH(ii-1).state(iL,1).T ;
                [gas,fluidH(ii),i,~] = hex_TQ(gas,[iL,i],fluidH(ii),[iL,1],eff,ploss,'hex',3,TCoutMIN); % Mode 3. 
            end
            iH = ii ;
        end
        iH = iH + 1 ;        
    end    
    if setTmax, PRch = ptop/pbot; end
    
    % REGENERATE (gas-gas)
    if Nrcp == 1
        [gas,~,i,~] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0);
    elseif Nrcp == 2
        gas.state(iL,iReg3).mdot = gas.state(iL,iReg1).mdot ;
        [gas,~,i,~] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg3],eff,ploss,'regen',0,0); % High-temp regenerator
        [gas,~,i,~] = hex_TQ(gas,[iL,i],gas,[iL,iReg2],eff,ploss,'regen',0,0); % Low-temp regenerator
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
        %[gas,i] = compexp(gas,[iL,i],eta,p_aim,1);   
        [CEXP(iN),gas,i] = compexp_func (CEXP(iN),gas,[iL,i],'Paim',p_aim) ; 
        
        % HEAT (gas-liquid)
        for ii = 1 : Ncld
            fluidC(ii).state(iL,1).T = CT(ii).A(iL).T; fluidC(ii).state(iL,1).p = CT(ii).A(iL).p;
            [fluidC(ii)] = update(fluidC(ii),[iL,1],1);
            if ii == 1
                [fluidC(ii),gas,~,i] = hex_TQ(fluidC(ii),[iL,1],gas,[iL,i],eff,ploss,'hex', 2, 1.0); % Mode 2: sCO2 mdot is known, cold fluid mdot is not
            else
                THoutMAX = fluidC(ii-1).state(iL,1).T ;
                [fluidC(ii),gas,~,i] = hex_TQ(fluidC(ii),[iL,1],gas,[iL,i],eff,ploss,'hex', 4, THoutMAX); % Mode 
            end
            iC = ii ;
        end        
        iC=iC+1;
    end
    
    % REGENERATE (gas-gas)
    if Nrcp == 1
        [~,gas,~,i] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0); 
    elseif Nrcp == 2
        [~,gas,~,~] = hex_TQ(gas,[iL,iReg1+1],gas,[iL,iReg2],eff,ploss,'regen',0,0); % Low-temp regenerator
        [~,gas,~,i] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg3],eff,ploss,'regen',0,0); 
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
for ii = 1 : Nhot
    [HT(ii)] = run_tanks(HT(ii),fluidH(ii),iL,Load,T0);
end

% Cold tanks
for ii = 1 : Ncld
    [CT(ii)] = run_tanks(CT(ii),fluidC(ii),iL,Load,T0);
end