%%% CCES PLANT LAYOUT DURING DISCHARGE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAES SIDE
%  ---------
%  B1: Pump inlet (liquid air tank)     - xb   - DCMP(1)
%  B2: Coupler B inlet (vs. E2)         - xb   - HX(6)
%  B3: Cold HEX                         - xb   - HX(5)
%  B4: Rejection unit inlet             - xb   - HX(4)
%  B5: Hot HEX inlet                    - xb   - HX(3)
%  B6: Turbine inlet                    - xb   - DEXP(1)
%  B7: Mixer inlet (vs. C3, makes A1)   - xb   - MIX(1)
%
%  SHARED
%  ------
%  A1: Rejection unit inlet             - xa   - HX(2)
%  A2: Hot HEX inlet                    - xa   - HX(1)
%  A3: Turbine inlet                    - xa   - DEXP(2)
%  A4: Atmosphere (turbine outlet)      - xa   - end
%
%  PTES SIDE
%  ---------
%  E1: Regen inlet (vs. E4)             - xe   - HX(7)
%  E2: Coupler B inlet (vs. B2)         - xe   - HX(6)
%  E3: Compressor inlet                 - xe   - DCMP(2)
%  E4: Regen inlet (vs. E1)             - xe   - HX(7)
%  E5: Auxiliary compressor inlet       - xc   - DCMP(4)
%  E5: Mixer inlet (vs. D3, makes C1)   - xe   - MIX(2)
%  E6: Mixer inlet (vs. B7, makes A1)   - xc   - MIX(1)
%  -->E7: END OF STREAM AFTER MIXING    - 0    - end
%
% where:
% A1:   xa = xb + xe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set stage indices
iB1 = 1+0;
iA1 = iB1+7;
iE1 = iA1+4;

iH   = 1;  % keeps track of the Hot fluid stream number
iM   = 1;  % keeps track of the Medium fluid stream number
iC   = 1;  % keeps track of the Cold fluid stream number
iR0  = iE1+7; 
iR   = iR0;% keeps track of the Air (heat rejection) stream number
iPMP = 1 ; % Keeps track of which pump is being used
iFAN = 1 ; % Keeps track of which fan is being used

% Read the pressure, temperature and enthalpy of the liquid air tank (LAT)
pLA = LAT.A(2).p;
TLA = LAT.A(2).T;
hLA = LAT.A(2).h;

% Set fixed conditions at LAES side
% Liquid air pump inlet
gas.state(iL,1).p = pLA;
gas.state(iL,1).h = hLA;
mLA = Load.mdot(iL);
gas.state(iL,1).mdot = mLA;
[gas] = update(gas,[iL,1],2);

% Perform pump compression. This is only done once to obtain conditions at
% the outlet
p_aim = pmax_LA*PRr;
[DCMP(1),gas,~] = compexp_func(DCMP(1),iL,gas,iB1,'Paim',p_aim,1);
TSA = gas.state(iL,2).T; % Temperature supercritical air at pump outlet

% Estimate temperature points between/around the Couplers
TB4dis = T0;                        % Coupler A air outlet
TB3dis = CT.A(1).T-5;               % between Coupler B and Coupler A
TB2dis = TSA + (CT.A(1).T-TSA)*(1-eff); % Coupler B air inlet
TR_dis = TB3dis/TB2dis;

% Set discharge pressure ratio for PTES side. PR_PTES_dis is the total after
% the auxiliary compressor ('hot' line), while PRC_dis is before the
% auxiliary compressor ('cold' line)
PR_PTES_dis = PR_PTES_chg*PRr;
PRC_dis = TR_dis^(1/phi_comp);
PRaux_dis = PR_PTES_dis/PRC_dis;

% Use temperature points to compute mass flow rates at the PTES side
nT      = 100;
TvA_dis = linspace(TB3dis,TB4dis,nT)';
TvB_dis = linspace(TB2dis,TB3dis,nT)';
pv      = pmax_LA.*ones(size(TvA_dis));
CpvA    = RPN('PT_INPUTS',pv,TvA_dis,'C',gas);
CpvB    = RPN('PT_INPUTS',pv,TvB_dis,'C',gas);
CpA     = mean(CpvA);
CpB     = mean(CpvB);
CpN     = RPN('PT_INPUTS',p0,300,'C',gas);
xb_dis  = 1.0;
xe_dis  = xb*CpB/CpN;
xa_dis  = xb_dis + xe_dis;

% Set fixed conditions at PTES side
% Regenerator inlet (warm, low p side - E1)
gas.state(iL,iE1).p    = p0;
gas.state(iL,iE1).T    = T0;
gas.state(iL,iE1).mdot = xe_dis*mLA;
[gas] = update(gas,[iL,iE1],1);

% Set initial guess conditions (PTES side)
% Regenerator inlet (cold, high p side - E4)
iE4 = iE1+3;
gas.state(iL,iE4).p    = p0*PRC_dis/(1-ploss);
gas.state(iL,iE4).T    = TB3dis;
gas.state(iL,iE4).mdot = xe_dis*mLA;
[gas] = update(gas,[iL,iE4],1);

% Set matrix of temperature and pressure points to test convergence
C_0 = [[gas.state(iL,:).T],[gas.state(iL,:).p]];
max_iter=100;
for counter=1:max_iter
    
    fprintf(1,['Discharging CCES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%% LAES - PUMP %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    % PUMP COMPRESSION (B1-->B2)
    p_aim = pmax_LA*PRr;
    [DCMP(1),gas,iB2] = compexp_func(DCMP(1),iL,gas,iB1,'Paim',p_aim,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PTES & COUPLERS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REGENERATOR (E1-->E2) & (E4-->E5)
    [HX(7),gas,iE2,gasX,iE5] = hex_func(HX(7),iL,gas,iE1,gas,iE4,0,0);
    gas.state(iL,iE5) = gasX.state(iL,iE5);
    gas.stage(iL,iE4) = gasX.stage(iL,iE4);
    
    % COUPLER B (E2->E3) (B2-->B3)
    [HX(6),gas,iE3,gasX,iB3] = hex_func(HX(6),iL,gas,iE2,gas,iB2,0,0);
    gas.state(iL,iB3) = gasX.state(iL,iB3);
    gas.stage(iL,iB2) = gasX.stage(iL,iB2);
    
    % COMPRESSION (E3-->E4)
    p_aim = p0*PRC_dis/(1-ploss);
    [DCMP(2),gas,iE4] = compexp_func (DCMP(2),iL,gas,iE3,'Paim',p_aim,1);
    
    % AUXILIARY COMPRESSION (E5-->E6)
    p_aim = p0*PR_PTES_dis/(1-ploss);
    [DCMP(4),gas,iE6] = compexp_func (DCMP(4),iL,gas,iE5,'Paim',p_aim,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% LAES - HEAT & TURBINE %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HEAT (B3-->B4)
    % Set fluidC
    fluidC.state(iL,iC).T = CT.B(iL).T; fluidC.state(iL,iC).p = CT.B(iL).p;
    fluidC.state(iL,iC).mdot = fluidC.state(1,iC).mdot;
    [fluidC] = update(fluidC,[iL,iC],1);
    % Run HEX
    [HX(5),fluidC,iC,gas,iB4] = hex_func(HX(5),iL,fluidC,iC,gas,iB3,0,0);
    % Run Pump
    [DPMP(iPMP),fluidC,iC] = compexp_func(DPMP(iPMP),iL,fluidC,iC,'Paim',fluidC.state(iL,1).p,1); %#ok<*SAGROW>
    iC=iC+1;iPMP=iPMP+1;
    
    % REJECT HEAT (external HEX) (B4-->B5)
    gas.state(iL,iR).T = T0; gas.state(iL,iR).p = p0; gas = update(gas,[iL,iR],1);
    [HX(4), gas, iR, gasX, iB5] = hex_func(HX(4),iL,gas,iR,gas,iB4,2,2.0);
    gas.state(iL,iB5) = gasX.state(iL,iB5);
    gas.stage(iL,iB4) = gasX.stage(iL,iB4);
    [DFAN(iFAN),gas,iR] = compexp_func (DFAN(iFAN),iL,gas,iR,'Paim',p0,1);
    iFAN=iFAN+1; iR=iR+1;
    
    % HEAT (B5-->B6)
    % Set fluidM
    fluidM.state(iL,iM).T = MT.B(iL).T; fluidM.state(iL,iM).p = MT.B(iL).p;
    fluidM.state(iL,iM).mdot = fluidM.state(1,iM).mdot;
    [fluidM] = update(fluidM,[iL,iM],1);
    % Run HEX
    [HX(3),fluidM,iM,gas,iB6] = hex_func(HX(3),iL,fluidM,iM,gas,iB5,0,0);
    % Run Pump
    [DPMP(iPMP),fluidM,iM] = compexp_func(DPMP(iPMP),iL,fluidM,iM,'Paim',fluidM.state(iL,1).p,1); %#ok<*SAGROW>
    iM=iM+1;iPMP=iPMP+1;
    
    % EXPANSION (B6-->B7)
    p_aim = gas.state(iL,iE6).p;
    [DEXP(1),gas,iB7] = compexp_func (DEXP(1),iL,gas,iB6,'Paim',p_aim, 1) ;
    
    % MIX (B7 & E6-->A1)
    [MIX(2),gas,iA1] = mix_streams2(MIX(2),gas,[iL,iB7],[iL,iE6]);
    iE7 = iE6+1;
    gas.state(iL,iE7) = gas.state(iL,iA1);
    gas.state(iL,iE7).mdot = 0;
    gas.stage(iL,iE7).type = 'end';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% CLOSE SHARED %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    % REJECT HEAT (external HEX) (A1-->A2)
    gas.state(iL,iR).T = T0; gas.state(iL,iR).p = p0; gas = update(gas,[iL,iR],1);
    [HX(2), gas, iA2, gasX, iR] = hex_func(HX(2),iL,gas,iA1,gas,iR,1,0.5);
    gas.state(iL,iR)   = gasX.state(iL,iR);
    gas.state(iL,iR-1) = gasX.state(iL,iR-1);
    gas.stage(iL,iR-1) = gasX.stage(iL,iR-1);
    [DFAN(iFAN),gas,iR] = compexp_func (DFAN(iFAN),iL,gas,iR,'Paim',p0,1);
    iFAN=iFAN+1; iR=iR+1;
    
    % HEAT (A2-->A3)
    % Set fluidH
    fluidH.state(iL,iH).T = HT.B(iL).T; fluidH.state(iL,iH).p = HT.B(iL).p;
    fluidH.state(iL,iH).mdot = fluidH.state(1,iH).mdot;
    [fluidH] = update(fluidH,[iL,iH],1);
    % Run HEX
    [HX(1),fluidH,iH,gas,iA3] = hex_func(HX(1),iL,fluidH,iH,gas,iA2,0,0);
    % Run Pump
    [DPMP(iPMP),fluidH,iH] = compexp_func(DPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1); %#ok<*SAGROW>
    iH=iH+1;iPMP=iPMP+1;
    
    % EXPANSION (A3-->A4)
    p_aim = p0;
    [DEXP(2),gas,iA4] = compexp_func (DEXP(2),iL,gas,iA3,'Paim',p_aim, 1) ;
    iLA_out = 1;
    iLA_in  = [];
    
    % Determine convergence and proceed
    C = [[gas.state(iL,:).T],[gas.state(iL,:).p]];
    err = abs((C(C~=0) - C_0(C~=0))./C(C~=0))*100; % percentage error
    fprintf('Convergence error = %f %%\n',max(err))
    convergence = all(err < 1e-3);
    %print_states(gas,iL,1:33,Load);
    %keyboard
    
    if convergence || counter==max_iter % is charge cycle converged?
        
        % Close working fluid streams and air (heat rejection) streams
        % (they are the same in this cycle!)
        gas.stage(iL,iA4).type   = 'end';
        iR_out = [iE1,iR0:3:(iR-1)]; iR_in  = [(iR0:3:(iR-1)) + 2,iA4];
        for i=iR_in, gas.stage(iL,i).type = 'end'; end
        gas = count_Nstg(gas);
        
        % Close storage fluid streams
        iH_out = 1:3:(iH-1); iH_in  = iH_out + 2;
        iM_out = 1:3:(iM-1); iM_in  = iM_out + 2;
        iC_out = 1:3:(iC-1); iC_in  = iC_out + 2;
        for i=iH_in, fluidH.stage(iL,i).type = 'end'; end
        for i=iM_in, fluidM.stage(iL,i).type = 'end'; end
        for i=iC_in, fluidC.stage(iL,i).type = 'end'; end
        fluidH = count_Nstg(fluidH);
        fluidM = count_Nstg(fluidM);
        fluidC = count_Nstg(fluidC);
                
        % Uncomment these lines to print states
        %{
        print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
        print_states(fluidM,iL,1:fluidM.Nstg(iL)+1,Load);
        print_states(fluidC,iL,1:fluidC.Nstg(iL)+1,Load);
        keyboard
        %}
        
        % Uncomment these lines to prints HEXs
        %{
        plot_hex(HX(1) ,iL,41,'K',true,'Shared hot HEX')
        plot_hex(HX(2) ,iL,42,'K',true,'Shared heat rejection')
        plot_hex(HX(3) ,iL,43,'K',true,'LAES hot HEX')
        plot_hex(HX(4) ,iL,44,'K',true,'LAES heat rejection')
        plot_hex(HX(5) ,iL,45,'K',true,'LAES cold HEX')
        plot_hex(HX(6), iL,46,'K',true,'Coupler B')
        plot_hex(HX(7) ,iL,47,'K',true,'PTES Regen')
        %keyboard
        %}
        
        % Exit loop
        break
        
    else
        C_0 = C;
        iH=1; iM=1; iC=1; iR=iR0; iFAN=1; iPMP=1;
        
    end
end
if counter==max_iter
    warning('Exiting JB_CHARGE cycle without having reached convergence');
end

% Compute total mass flow rates of the hot and cold storage fluids and
% compute the end of discharge time (stop when one tank becomes empty)
[MdotH]  = total_mdot(fluidH, iL, iH_out);
t_disH   = HT.B(iL).M/MdotH;
[MdotM]  = total_mdot(fluidM, iL, iM_out);
t_disM   = MT.B(iL).M/MdotM;
[MdotC]  = total_mdot(fluidC, iL, iC_out);
t_disC   = CT.B(iL).M/MdotC;
[MdotLA] = total_mdot(gas, iL, iLA_out);
t_disLA  = LAT.A(iL).M/MdotLA;
Load.time(iL) = min([Load.time(iL),t_disH,t_disM,t_disC,t_disLA])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
HT = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
% Medium tanks
MT = run_tanks(MT,iL,fluidM,iM_out,iM_in,Load,T0);
% Cold tanks
CT = run_tanks(CT,iL,fluidC,iC_out,iC_in,Load,T0);
% Atmospheric tanks
AT = run_single_tank(AT,iL,gas,iR_out,iR_in,Load,T0);
% Liquid air tanks
LAT = run_single_tank(LAT,iL,gas,iLA_out,iLA_in,Load,T0);

% PLOT DISCHARGE CYCLE
%PLOT_CYCLE
%keyboard