% Set stage indices
iG = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the heat rejection stream number
iA = 1;  % keeps track of the Air (heat rejection) stream number

% Initial guess of charge conditions
% Compressor inlet (regenerator hot outlet)
gas.state(iL,1).p    = pbot; gas.state(iL,1).T = T1;
gas.state(iL,1).mdot = Load.mdot(iL);
[gas] = update(gas,[iL,1],1);
% Regenerator inlet indeces
iReg1 = 1 + Nc_ch*2;         % index hot regenerator inlet (after compression + cooling)
iReg2 = iReg1 + 2 + Ne_ch*2; % index cold regenerator inlet (after regeneration + heat_reject + expansion + heating)
gas.state(iL,iReg2).T = T0;
gas.state(iL,iReg2).p = pbot;
gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
[gas] = update(gas,[iL,iReg2],1);

% Set matrix of temperature and pressure points to test convergence
A_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
while 1
    
    fprintf(1,'Hello charge PTES\n')
    
    for iN = 1:Nc_ch
        % COMPRESS
        if setTmax
            T_aim = Tmax;
            [CCMP(iN),gas,iG] = compexp_func (CCMP(iN),iL,gas,iG,'Taim',T_aim) ; %#ok<*SAGROW>
        else
            PRc = (pmax/gas.state(iL,iG).p)^(1/(Nc_ch+1-iN)); % stage compression pressure ratio
            p_aim = gas.state(iL,iG).p*PRc;
            [CCMP(iN),gas,iG] = compexp_func (CCMP(iN),iL,gas,iG,'Paim',p_aim) ; 
        end
        ptop  = gas.state(iL,iG).p;
        
        % COOL (gas-liquid)
        fluidH.state(iL,iH).T = HT.A(iL).T; fluidH.state(iL,iH).p = HT.A(iL).p;
        [fluidH] = update(fluidH,[iL,iH],1);
        if new_hex_calls
            %[HX1,gas,iG,fluidH,iH] =
            %hex_func(HX1,iL,gas,iG,fluidH,iH,1,1.0); % Original call
            [HXh,gas,iG,fluidH,iH] = hex_func(HXh,iL,gas,iG,fluidH,iH,1,1.0); % New call using hx_class
        else
            [gas,fluidH,iG,iH,HX] = hex_TQ(gas,[iL,iG],fluidH,[iL,iH],eff,ploss,'hex',1,1.0);
        end
        iH=iH+1;
    end    
    if setTmax, PRch = ptop/pbot; end
    
    % REGENERATE (gas-gas)
    if new_hex_calls
        %[REGEN,gas,iG,~,~] = hex_func(REGEN,iL,gas,iReg1,gas,iReg2,0,0);
        [RCP,gas,iG,~,~] = hex_func(RCP,iL,gas,iReg1,gas,iReg2,0,0);% New call using hx_class
    else
        [gas,~,iG,~] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0);
    end
        
    % REJECT HEAT (external HEX)
    T_aim = environ.T0;
    [gas,environ,iG,iE] = hex_set(gas,[iL,iG],environ,[iL,iE],T_aim,eff,ploss);
    
    for iN = 1:Ne_ch
        % EXPAND
        PRe = (gas.state(iL,iG).p/pbot)^(1/(Ne_ch+1-iN)); % stage expansion pressure ratio
        p_aim = gas.state(iL,iG).p/PRe;
        [CEXP(iN),gas,iG] = compexp_func (CEXP(iN),iL,gas,iG,'Paim',p_aim) ; 
        
        % HEAT (gas-liquid)
        fluidC.state(iL,iC).T = CT.A(iL).T; fluidC.state(iL,iC).p = CT.A(iL).p;
        [fluidC] = update(fluidC,[iL,iC],1);
        if new_hex_calls
            %[HX,fluidC,iC,gas,iG] = hex_func(HX,iL,fluidC,iC,gas,iG,2,1.0);
            [HXc,fluidC,iC,gas,iG] = hex_func(HXc,iL,fluidC,iC,gas,iG,2,1.0); % New call using hx_class
        else
            [fluidC,gas,iC,iG] = hex_TQ(fluidC,[iL,iC],gas,[iL,iG],eff,ploss,'hex',2,1.0);
        end
        iC=iC+1;
    end
    
    % REGENERATE (gas-gas)
    if new_hex_calls
        %[REGEN,~,~,gas,iG] = hex_func(REGEN,iL,gas,iReg1,gas,iReg2,0,0);
        [RCP,~,~,gas,iG] = hex_func(RCP,iL,gas,iReg1,gas,iReg2,0,0); % New call using hx_class
    else
        [~,gas,~,iG] = hex_TQ(gas,[iL,iReg1],gas,[iL,iReg2],eff,ploss,'regen',0,0);
    end
    
    % Determine convergence and proceed
    A = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    
    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is charge cycle converged?
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
        iG=1; iH=1; iC=1; iE=1;
    end
end

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
HT = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
% Cold tanks
CT = run_tanks(CT,iL,fluidC,iC_out,iC_in,Load,T0);
% Atmospheric tanks
iA_out = 0; iA_in = 0;
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);