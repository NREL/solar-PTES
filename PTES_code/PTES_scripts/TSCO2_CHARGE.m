% This script runs the charging cycle of a "Time-Shifted Recompression sCO2
% power cycle.

% Before the charging cycle can be run, have to have an idea of what the
% normal recompression-sCO2 cycle looks like. Therefore run a solar cycle
% to charge the hot tanks, then a recompression cycle. Then extract the key
% information, such as temperatures and mass flow rates. Use these to
% design the TSRC-sCO2 charging cycle.

% Only want to do this all once - very first time charge cycle is called

if iL == 1

    % First reset the solar tanks
    HT(1) = reset_tanks(HT(1),TH_dis0(1),p0,0,TH_chg0(1),p0,MH_dis0(1),T0);
    
    HT(1).A(2) = HT(1).A(1); % Not sure these are necessary
    HT(1).B(2) = HT(1).B(1);
    
    % Now run the sCO2-recompression cycle
    Load.mode(iL) = 5 ;
    Load.type(iL) = 'rcmpCO2' ;
    Nhot = 1 ;
    sCO2_RECOMP
    
    % Extract key design points from the sCO2-recompression cycle
    TH_dis0(2) = gas.state(1,5).T ;
    TH_chg0(2) = gas.state(1,6).T ;
    HT(2)  = double_tank_class(fluidH(2),TH_dis0(2),p0,MH_dis0(2),TH_chg0(2),p0,MH_chg0(2),T0,Load.num+1); %hot double tank
    
    TC_dis0 =  gas.state(1,3).T ;
    TC_chg0 =  gas.state(1,4).T ;
    CT  = double_tank_class(fluidC,TC_dis0,p0,MC_dis0,TC_chg0,p0,MC_chg0,T0,Load.num+1); %cold double tank
    
    Load.mdot(iL) = gas.state(1,11).mdot ;
    
    % Reset stage and load counters to over-write the temporary data that was
    % just generated
    iL = 1 ;
    Load.mode(iL) = 6 ;
    Load.type(iL) = 'chgTSCO2' ;
    Nhot = 2 ;
    
    % Reset cold tanks
    for ir = 1 : Ncld
        CT(ir) = reset_tanks(CT(ir),TC_dis0(ir),p0,MC_dis0(ir),TC_chg0(ir),p0,MC_chg0(ir),T0);
    end
    
    % Reset hot tanks
    for ir = 1 : Nhot
        HT(ir) = reset_tanks(HT(ir),TH_dis0(ir),p0,MH_dis0(ir),TH_chg0(ir),p0,MH_chg0(ir),T0);
    end
    
    for ir = 1:length(fluidH)
        fluidH(ir) = reset_fluid(fluidH(ir)); %#ok<*SAGROW>
    end

end

% Set stage indices
iG = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the heat rejection stream number
iA = 1;  % keeps track of the Air (heat rejection) stream number

% Initial guess of charge conditions
% Compressor inlet (cold storage outlet)
gas = reset_fluid(gas);
gas.state(iL,1).p    = pbot; gas.state(iL,1).T = TC_dis0;
gas.state(iL,1).mdot = Load.mdot(iL);
[gas] = update(gas,[iL,1],1);

% Cold storage inlet
gas.state(iL,4).T = TC_chg0;
gas.state(iL,4).p = pbot;
gas.state(iL,4).mdot = gas.state(iL,1).mdot;
[gas] = update(gas,[iL,4],1);

% Set matrix of temperature and pressure points to test convergence
A_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
while 1
    
    fprintf(1,'Hello charge TSRC-sCO2 PTES\n')
    
    % COMPRESS
    T_aim = TH_chg0(2);
    [CCMP(1),gas,iG] = compexp_func (CCMP(1),iL,gas,iG,'Taim',T_aim) ;
    ptop  = gas.state(iL,iG).p;
    
    % COOL (gas-liquid)
    fluidH(2).state(iL,iH).T = HT(2).A(iL).T; fluidH(2).state(iL,iH).p = HT(2).A(iL).p;
    [fluidH(2)] = update(fluidH(2),[iL,iH],1);
    [HX(1),gas,iG,fluidH(2),iH] = hex_func(HX(1),iL,gas,iG,fluidH(2),iH,1,1.0);
    iH=iH+1;
    
    PRch = ptop/pbot ;
        
    % REJECT HEAT (external HEX)
    if Lcld
        T_aim = environ.T0;
    else
        T_aim = gas.state(iL,iG).T - 0.1 ; % Need to reject a little heat for the energy balance
    end
    [gas,environ,iG,iE] = hex_set(gas,[iL,iG],environ,[iL,iE],T_aim,eff,ploss);
    
    % EXPAND
    PRe = (gas.state(iL,iG).p/pbot)^(1/(Ne_ch+1-iN)); % stage expansion pressure ratio
    p_aim = gas.state(iL,iG).p/PRe;
    [CEXP(1),gas,iG] = compexp_func (CEXP(1),iL,gas,iG,'Paim',p_aim) ;
    
    % HEAT (gas-liquid)
    fluidC.state(iL,iC).T = CT.A(iL).T; fluidC.state(iL,iC).p = CT.A(iL).p;
    [fluidC] = update(fluidC,[iL,iC],1);
    [HX(2),fluidC,iC,gas,iG] = hex_func(HX(2),iL,fluidC,iC,gas,iG,2,1.0);
    iC=iC+1;
    
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
        for i=iH_in, fluidH(2).stage(iL,i).type = 'end'; end
        for i=iC_in, fluidC.stage(iL,i).type = 'end'; end
        fluidH(2) = count_Nstg(fluidH(2));
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
% Hot tanks. Have to do something about hot tank 1 - e.g. the one charged
% by solar???

% Solar tank behaviour
HT(1) = reset_tanks(HT(1),TH_dis0(1),p0,0,TH_chg0(1),p0,MH_dis0(1),T0);
HT(1).A(2) = HT(1).A(1);
HT(1).B(2) = HT(1).B(1);

HT(2) = run_tanks(HT(2),iL,fluidH(2),iH_out,iH_in,Load,T0);
% Cold tanks
CT = run_tanks(CT,iL,fluidC,iC_out,iC_in,Load,T0);


