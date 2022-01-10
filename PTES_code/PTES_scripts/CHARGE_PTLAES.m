%%% INTEGRATED PTES-LAES PLANT LAYOUT DURING CHARGE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SHARED
%  ------
%  A1: Compressor inlet (atmosphere)        - xa   - CCMP(1)
%  A2: Hot HEX inlet                        - xa   - HX(1)
%  A3: Rejection unit inlet                 - xa   - HX(2)
%  A4: Separator inlet (makes b1 & c1)      - xa   - ...
%
%  LAES SIDE
%  ---------
%  B1: Compressor inlet                     - xb   - CCMP(2)
%  B2: Hot HEX inlet                        - xb   - HX(3)
%  B3: Rejection unit inlet                 - xb   - HX(4)
%  B4: Coupler A inlet (vs. d2)             - xb   - HX(5)
%  B5: Coupler B inlet (vs. e3)             - xb   - HX(6)
%  B6: Expander inlet                       - xb   - CEXP(1)
%  B7: Liquid air tank inlet                - xb   - ...
%
%  PTES SIDE
%  ---------
%  C1: Separator inlet (makes d1 & e1)      - xc   - ...
%
%  D1: Expander inlet                       - xd   - CEXP(2)
%  D2: Coupler A inlet (vs. b4)             - xd   - HX(5)
%  D3: End of stream xd                     - xd   - end
%
%  E1: Regen inlet (vs. e4)                 - xe   - HX(7)
%  E2: Expander inlet                       - xe   - CEXP(3)
%  E3: Coupler B inlet (vs. b5)             - xe   - HX(6)
%  E4: Regen inlet (vs. e1)                 - xe   - HX(7)
%  E5: End of stream xe                     - xe   - end
%
% where:
% A4:   xa = xb + xc
% C1:   xc = xd + xe
% so:   xa = xb + xd + xe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set stage indices
iA1 = 1+0;
iB1 = iA1+4;
iC1 = iB1+7;
iD1 = iC1+1;
iE1 = iD1+3;
iA = iA1; iB = iB1; iC = iC1; iD = iD1; iE = iE1;

iH   = 1;  % keeps track of the Hot fluid stream number
iM   = 1;  % keeps track of the Medium fluid stream number
iR0  = iE1+5;
iR   = iR0;% keeps track of the heat Rejection air stream number
iPMP = 1 ; % Keeps track of which pump is being used
iFAN = 1 ; % Keeps track of which fan is being used

% Set charge pressure ratio
gamma    = 1.4;
phi_exp  = eta*(gamma-1)/gamma;
phi_comp = (gamma-1)/(gamma*eta);
%TR_exp   = Tpnts(2)/Tpnts(1);
%PRC_chg  = (TR_exp/eff)^(1/phi_exp); % pressure ratio PTES 'cold' line
%PR_PTES_chg  = min([(Tmax/T0)^(1/phi_comp),PRC_chg]); % pressure ratio PTES side

PR_PTES_chg  = sqrt(pmax_LA/p0)/(1-ploss)^2/0.98; % pressure ratio PTES side
TR_exp = PR_PTES_chg^(phi_exp)*eff;
TB4chg = T0;
TB5chg = TB4chg/TR_exp;
TB6chg = TB5chg/TR_exp;

% Obtain mass flow rates at PTES side. This requires computation of
% specific heat capacity of supercritical air between points b6 and b4.
% Assume 'logspace distribution' of temperatures around heat exchanger
% couplers (2 couplers, 3 temperature points).
%Tpnts = logspace(log10(TLA_obj),log10(T0),3);
%TB6chg = Tpnts(1);
%TB5chg = Tpnts(2);
%TB4chg = Tpnts(3);
nT    = 100;
TvA   = linspace(TB5chg,TB4chg,nT)';
TvB   = linspace(TB6chg,TB5chg,nT)';
pv    = pmax_LA.*ones(size(TvA));
CpvA  = RPN('PT_INPUTS',pv,TvA,'C',gas);
CpvB  = RPN('PT_INPUTS',pv,TvB,'C',gas);
%figure(10); plot([TvA;TvB],[CpvA;CpvB]); xlabel('Temperature, K'); ylabel('Cp, J/kg/K');
%legend(sprintf('%.1f bar',pmax_LA/1e5));
CpA = mean(CpvA);
CpB = mean(CpvB);
CpN = RPN('PT_INPUTS',p0,300,'C',gas);
mLA = Load.mdot(iL);
xb  = 1.0;
xd  = xb*CpA/CpN;
xe  = xb*CpB/CpN;
xc  = xd + xe;
xa  = xb + xc;

% Initial guess of charge conditions
% Compressor inlet
gas.state(iL,iA1).p = p0;
gas.state(iL,iA1).T = T0;
gas.state(iL,iA1).mdot = xa*mLA;
[gas] = update(gas,[iL,iA1],1);
% PTES Regenerator inlet (cold, low p side - e4)
iE4 = iE1+3;
gas.state(iL,iE4).p    = p0/(1-ploss);
gas.state(iL,iE4).T    = TB5chg;
gas.state(iL,iE4).mdot = xe*mLA;
[gas] = update(gas,[iL,iE4],1);


% Set matrix of temperature and pressure points to test convergence
C_0 = [[gas.state(iL,:).T],[gas.state(iL,:).p]];
max_iter=100;
for counter=1:max_iter
    
    fprintf(1,['Charging CCES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PTES & LAES %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    % FIRST COMPRESSION (A1-->A2)
    p_aim = p0*PR_PTES_chg;
    [CCMP(1),gas,iA] = compexp_func(CCMP(1),iL,gas,iA,'Paim',p_aim,1);
    
    % COOL (A2-->A3)
    % Set fluidH
    fluidH.state(iL,iH).T = HT.A(iL).T; fluidH.state(iL,iH).p = HT.A(iL).p;
    [fluidH] = update(fluidH,[iL,iH],1);
    % Run HEX
    [HX(1),gas,iA,fluidH,iH] = hex_func(HX(1),iL,gas,iA,fluidH,iH,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidH,iH] = compexp_func(CPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
    iH=iH+1;iPMP=iPMP+1;
    
    % REJECT HEAT (external HEX) (A3-->A4)
    gas.state(iL,iR).T = T0; gas.state(iL,iR).p = p0; gas = update(gas,[iL,iR],1);
    [HX(2), gas, iA, gasX, iR] = hex_func(HX(2),iL,gas,iA,gas,iR,1,0.5);
    gas.state(iL,iR)   = gasX.state(iL,iR);
    gas.state(iL,iR-1) = gasX.state(iL,iR-1);
    gas.stage(iL,iR-1) = gasX.stage(iL,iR-1);
    [CFAN(iFAN),gas,iR]   = compexp_func (CFAN(iFAN),iL,gas,iR,'Paim',p0,1);
    iFAN=iFAN+1; iR=iR+1;
    
    % SEPARATE (A4-->B1 & C1)
    x_aim = xc/xa; % mass flow rate fraction of secondary stream
    [gas,iB] = split_stream(gas,iL,iA,iC1,x_aim);
    
    %%%%%%%%%%%%%%%%
    %%%%% LAES %%%%%
    %%%%%%%%%%%%%%%%
    % SECOND COMPRESSION (B1-->B2)
    p_aim = pmax_LA;
    [CCMP(2),gas,iB] = compexp_func(CCMP(2),iL,gas,iB,'Paim',p_aim,1);
    
    % COOL (B2-->B3)
    % Set fluidM
    fluidM.state(iL,iM).T = MT.A(iL).T; fluidM.state(iL,iM).p = MT.A(iL).p;
    [fluidM] = update(fluidM,[iL,iM],1);
    % Run HEX
    [HX(3),gas,iB,fluidM,iM] = hex_func(HX(3),iL,gas,iB,fluidM,iM,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidM,iM] = compexp_func(CPMP(iPMP),iL,fluidM,iM,'Paim',fluidM.state(iL,1).p,1); %#ok<*SAGROW>
    iM=iM+1;iPMP=iPMP+1;
    
    % REJECT HEAT (external HEX) (B3-->B4)
    gas.state(iL,iR).T = T0; gas.state(iL,iR).p = p0; gas = update(gas,[iL,iR],1);
    [HX(4), gas, iB, gasX, iR] = hex_func(HX(4),iL,gas,iB,gas,iR,1,0.5);
    gas.state(iL,iR)   = gasX.state(iL,iR);
    gas.state(iL,iR-1) = gasX.state(iL,iR-1);
    gas.stage(iL,iR-1) = gasX.stage(iL,iR-1);
    [CFAN(iFAN),gas,iR]   = compexp_func (CFAN(iFAN),iL,gas,iR,'Paim',p0,1);
    iFAN=iFAN+1; iR=iR+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PTES & COUPLERS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SEPARATE (C1-->D1 & E1)
    x_aim = xe/xc; % mass flow rate fraction of secondary stream
    [gas,iD] = split_stream(gas,iL,iC1,iE1,x_aim);
    
    % EXPANSION (D1-->D2)
    p_aim = p0/(1-ploss);
    [CEXP(2),gas,iD] = compexp_func (CEXP(2),iL,gas,iD,'Paim',p_aim,1);
    
    % COUPLER A (B4-->B5) (D2->D3)
    [HX(5),gas,iB,gasX,iD] = hex_func(HX(5),iL,gas,iB,gas,iD,0,0);
    gas.state(iL,iD)   = gasX.state(iL,iD);
    gas.stage(iL,iD-1) = gasX.stage(iL,iD-1);
    
    % REGENERATOR (E1--E2) & (E4-->E5)
    [HX(7),gas,iE,gasX,iE5] = hex_func(HX(7),iL,gas,iE1,gas,iE4,0,0);
    gas.state(iL,iE5) = gasX.state(iL,iE5);
    gas.stage(iL,iE4) = gasX.stage(iL,iE4);
    
    % EXPANSION (E2-->E3)
    p_aim = p0/(1-ploss)^2;
    [CEXP(3),gas,iE] = compexp_func (CEXP(3),iL,gas,iE,'Paim',p_aim,1);
    Te3 = gas.state(iL,iE).T;
    
    % COUPLER B (B5-->B6) (E3->E4)
    [HX(6),gas,iB,gasX,iE4] = hex_func(HX(6),iL,gas,iB,gas,iE,0,0);
    gas.state(iL,iE4)   = gasX.state(iL,iE4);
    gas.stage(iL,iE4-1) = gasX.stage(iL,iE4-1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%% CLOSE LAES %%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    % LIQUID AIR EXPANSION (B6-->B7)
    pLA = max([RPN('QT_INPUTS',0.0,gas.state(iL,iB).T,'P',gas)*1.2,p0]);
    [CEXP(1),gas,iB] = compexp_func (CEXP(1),iL,gas,iB,'Paim',pLA, 1) ;
    TLA = gas.state(iL,iB).T;
    hLA = gas.state(iL,iB).h;
    iLA_in  = iB;
    iLA_out = [];
    
    % Determine convergence and proceed
    C = [[gas.state(iL,:).T],[gas.state(iL,:).p]];
    err = abs((C(C~=0) - C_0(C~=0))./C(C~=0))*100; % percentage error
    fprintf('Convergence error = %f %%\n',max(err))
    convergence = all(err < 1e-3);
    %print_states(gas,iL,1:28,Load);
    %keyboard
    
    if convergence || counter==max_iter % is charge cycle converged?

        % Close working fluid streams and air (heat rejection) streams
        % (they are the same in this cycle!)
        gas.stage(iL,iB).type  = 'end';
        gas.stage(iL,iD).type  = 'end';
        gas.stage(iL,iE5).type = 'end';
        iR_out = [1,iR0:3:(iR-1)]; iR_in  = [iD,iE5,(iR0:3:(iR-1)) + 2];
        for i=iR_in, gas.stage(iL,i).type = 'end'; end
        gas = count_Nstg(gas);
        
        % Close storage fluid streams
        iH_out = 1:3:(iH-1); iH_in  = iH_out + 2;
        iM_out = 1:3:(iM-1); iM_in  = iM_out + 2;
        for i=iH_in, fluidH.stage(iL,i).type = 'end'; end
        for i=iM_in, fluidM.stage(iL,i).type = 'end'; end
        fluidH = count_Nstg(fluidH);
        fluidM = count_Nstg(fluidM);
                
        % Uncomment these lines to print states
        %{
        print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
        print_states(fluidM,iL,1:fluidM.Nstg(iL)+1,Load);
        keyboard
        %}
        
        % Uncomment these lines to prints HEXs
        %{
        plot_hex(HX(1) ,1,41,'K',true,'Shared hot HEX')
        plot_hex(HX(2) ,1,42,'K',true,'Shared heat rejection')
        plot_hex(HX(3) ,1,43,'K',true,'LAES hot HEX')
        plot_hex(HX(4) ,1,44,'K',true,'LAES heat rejection')
        plot_hex(HX(5) ,1,45,'K',true,'Coupler A')
        plot_hex(HX(6), 1,46,'K',true,'Coupler B')
        plot_hex(HX(7) ,1,47,'K',true,'PTES Regen')
        %keyboard
        %}
        
        % Exit loop
        break
        
    else
        C_0 = C;
        iA = iA1; iB = iB1; iC = iC1; iD = iD1; iE = iE1;
        iH=1; iM=1; iR=iR0; iFAN=1; iPMP=1;
        
    end
end
if counter==max_iter
    warning('Exiting JB_CHARGE cycle without having reached convergence');
end

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
HT = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
% Cold tanks
MT = run_tanks(MT,iL,fluidM,iM_out,iM_in,Load,T0);
% Atmospheric tank
AT = run_single_tank(AT,iL,gas,iR_out,iR_in,Load,T0);
% Liquid air tank
LAT = run_single_tank(LAT,iL,gas,iLA_out,iLA_in,Load,T0);


% PLOT CHARGE CYCLE
%PLOT_CYCLE
%keyboard