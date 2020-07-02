% Set stage indices
iG = 1;  % keeps track of the gas stage number
iH = ones(1,Nhot);  % keeps track of the Hot fluid stream number
iC = ones(1,Ncld); % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the heat rejection stream number
iA = 1; % Keeps track of air stream

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
            [CCMP(iN),gas,iG] = compexp_func (CCMP(iN),iL,gas,iG,'Taim',T_aim,design_mode) ; 
        else
            PRc = (pmax/gas.state(iL,iG).p)^(1/(Nc_ch+1-iN));
            p_aim = gas.state(iL,iG).p*PRc;
            [CCMP(iN),gas,iG] = compexp_func (CCMP(iN),iL,gas,iG,'Paim',p_aim,design_mode) ; 
        end
        ptop  = gas.state(iL,iG).p;
        
        % COOL (gas-liquid)
        for ii = 1 : Nhot
            fluidH(ii).state(iL,iH(ii)).T = HT(ii).A(iL).T; fluidH(ii).state(iL,iH(ii)).p = HT(ii).A(iL).p; %#ok<*SAGROW>
            [fluidH(ii)] = update(fluidH(ii),[iL,iH(ii)],1);
            if ii == 1
                [HX(ii),gas,iG,fluidH(ii),iH(ii)] = hex_func(HX(ii),iL,gas,iG,fluidH(ii),iH(ii),1,1.0); % Mode 1: don't know hot fluid mdot. Crat = 1
            else
                TCoutMIN = fluidH(ii-1).state(iL,1).T ;
                [HX(ii),gas,iG,fluidH(ii),iH(ii)] = hex_func(HX(ii),iL,gas,iG,fluidH(ii),iH(ii),3,TCoutMIN); % Mode 3.
            end
            iH(ii) = iH(ii) + 1;
        end
    end    
    if setTmax, PRch = ptop/pbot; end
    
    % REGENERATE (gas-gas)
    if Nrcp == 1
        [HX(Nhot+Ncld+1),gas,iG,~,~] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg2,0,0);
    elseif Nrcp == 2
        gas.state(iL,iReg3).mdot = gas.state(iL,iReg1).mdot ;
        [HX(Nhot+Ncld+1),gas,iG,~,~] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg3,0,0); % High-temp regenerator
        [HX(Nhot+Ncld+2),gas,iG,~,~] = hex_func(HX(Nhot+Ncld+2),iL,gas,iG,   gas,iReg2,0,0); % Low-temp regenerator
    end
    
    % May wish to make cold store as cold as possible, or avoid rejecting
    % heat here to avoid the worst of the CO2 c_p variation
    if Lcld
        T_aim = environ.T0 + T0_inc;
    else
        T_aim = gas.state(iL,iG).T - 0.1 ;
    end
    % REJECT HEAT (external HEX)
    [gas,environ,iG,iE] = hex_set(gas,[iL,iG],environ,[iL,iE],T_aim,eff,ploss);
    
    for iN = 1:Ne_ch
        % EXPAND
        PRe = (gas.state(iL,iG).p/pbot)^(1/(Ne_ch+1-iN)); %expansion pressure ratio
        p_aim = gas.state(iL,iG).p/PRe;
        [CEXP(iN),gas,iG] = compexp_func (CEXP(iN),iL,gas,iG,'Paim',p_aim,design_mode) ; 
        
        % HEAT (gas-liquid)
        for ii = 1 : Ncld
            fluidC(ii).state(iL,iC(ii)).T = CT(ii).A(iL).T; fluidC(ii).state(iL,iC(ii)).p = CT(ii).A(iL).p;
            [fluidC(ii)] = update(fluidC(ii),[iL,iC(ii)],1);
            if ii == 1
                [HX(Nhot+ii),fluidC(ii),iC(ii),gas,iG] = hex_func(HX(Nhot+ii),iL,fluidC(ii),iC(ii),gas,iG,2,1.0); % Mode 2: sCO2 mdot is known, cold fluid mdot is not
            else
                THoutMAX = fluidC(ii-1).state(iL,1).T ;
                [HX(Nhot+ii),fluidC(ii),iC(ii),gas,iG] = hex_func(HX(Nhot+ii),iL,fluidC(ii),iC(ii),gas,iG,4,THoutMAX);
            end
            iC(ii) = iC(ii) + 1;
        end
    end
    
    % REGENERATE (gas-gas)
    if Nrcp == 1
        [HX(Nhot+Ncld+1),~,~,gas,iG] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg2,0,0);
    elseif Nrcp == 2
        [HX(Nhot+Ncld+2),~,~,gas,~] = hex_func(HX(Nhot+Ncld+2),iL,gas,iReg1+1,gas,iReg2,0,0);
        [HX(Nhot+Ncld+1),~,~,gas,iG] = hex_func(HX(Nhot+Ncld+1),iL,gas,iReg1,gas,iReg3,0,0);
    end
    
    % Determine convergence and proceed
    A = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is charge cycle converged?
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


% Atmospheric tanks
iA_out = 0; iA_in = 0;
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);