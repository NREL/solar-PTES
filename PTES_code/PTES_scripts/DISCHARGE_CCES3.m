%%% CCES PLANT LAYOUT DURING DISCHARGE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAES SIDE
%  ---------
%  L1: Pump inlet (liquid air tank)     - DCMP(1)
%  L2: Coupler C inlet (vs. P6)         - HX(6)
%  L3: Coupler B inlet (vs. P22)        - HX(5)
%  L4: Coupler A inlet (vs. P18)        - HX(4)
%  L5: Hot HEX inlet                    - HX(2)
%  L6: Turbine inlet                    - DEXP(1)
%  L7: Rejection unit inlet             - HX(3)
%  L8: Hot HEX inlet                    - HX(1) 
%  L9: Turbine inlet                    - DEXP(2)
%  L10: Atmosphere (turbine outlet)     - end
%
%  PTES SIDE
%  ---------
%  P1: Regen A inlet (vs. P14)                  - xa   - HX(8)
%  P2: Rejection unit inlet                     - xa   - HX(9)
%  P3: Separator inlet (makes P4 & P18)         - xa   - ...
%  P4: Regen B inlet (vs. P11)                  - xb   - HX(10)
%  P5: Separator inlet (makes P6 & P22)         - xb   - ...
%  P6: Regen C inlet (vs. P9)                   - xc   - HX(11)
%  P7: Coupler C inlet (vs. L2)                 - xc   - HX(6)
%  P8: Compressor inlet                         - xc   - DCMP(2)
%  P9: Regen C inlet (vs. P6)                   - xc   - HX(11)
% P10: Mixer inlet (mixes with P24, makes P10)  - xc   - MIX(3)
% P11: Regen B inlet (vs. P4)                   - xb   - HX(10)
% P12: Mixer inlet (mixes with P20, makes P13)  - xb   - MIX(4)
% P13: Auxiliary compressor inlet               - xa   - DCMP(5)
% P14: Regen A inlet (vs. P1)                   - xa   - HX(8)
% P15: Hot HEX inlet                            - xa   - HX(7)
% P16: Turbine inlet                            - xa   - DEXP(3)
% P17: Closure of main cycle stream = P1        - xa   - end
% P18: Coupler A inlet (vs. L4)                 - xe   - HX(4)
% P19: Compressor inlet                         - xe   - DCMP(4)
% P20: Mixer inlet (mixes with P12, makes P13)  - xe   - MIX(4)
% P21: End of stream xe                         - 0    - end
% P22: Coupler B inlet (vs. L3)                 - xd   - HX(5)
% P23: Compressor inlet                         - xd   - DCMP(3)
% P24: Mixer inlet (mixes with P9, makes P10)   - xd   - MIX(3)
% P25: End of stream xd                         - 0    - end
%
% where:
% P3:   xa = xb + xe
% P4:   xb = xc + xd
% P10:  xc + xd = xb
% P12:  xb + xe = xa
% so:   xa = xc + xd + xe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set stage indices
iG1  = 1;  % keeps track of the gas stage number (LAES)
iP1  = 1;  % keeps track of the gas stage number (PTES)
iH   = 1;  % keeps track of the Hot fluid stream number
iM   = 1;  % keeps track of the Medium fluid stream number
iA   = 1;  % keeps track of the Air (heat rejection) stream number
iPMP = 1 ; % Keeps track of which pump is being used
iFAN = 1 ; % Keeps track of which fan is being used

% Read the pressure, temperature and enthalpy of the liquid air tank (LAT)
pLA = LAT.B(2).p;
TLA = LAT.B(2).T;
hLA = LAT.B(2).h;

% Set discharge pressure ratio for PTES side. PRH_dis is the total after
% the auxiliary compressor ('hot' line), while PRC_dis is before the
% auxiliary compressor ('cold' line)
PRH_dis = PRH_chg^(1/eta^2);%PRH_chg*PRr;
Tpnts = logspace(log10(TLA),log10(T0),4);
TR_comp  = Tpnts(2)/Tpnts(1);
PRC_dis = (TR_comp*eff)^(1/phi_comp)/(1-ploss)^2; % pressure ratio PTES 'cold' line
PRC_dis = min([PRC_dis,PRH_dis/1.01]);
PRaux_dis = PRH_dis/PRC_dis;

% Set conditions at liquid air pump inlet
gas1.state(iL,1).p = pLA;
gas1.state(iL,1).h = hLA;
mLA = Load.mdot(iL);
gas1.state(iL,1).mdot = mLA;
[gas1] = update(gas1,[iL,1],2);

% Perform pump compression. This is only done once to obtain conditions at
% the outlet
p_aim = pmax_LA;
[DCMP(1),gas1,~] = compexp_func(DCMP(1),iL,gas1,iG1,'Paim',p_aim,1);

% Estimate temperature points between/around the Couplers
TR_dis = PRC_dis^(phi_comp)/eff;
TL2dis = gas1.state(iL,2).T;        % Coupler C air inlet
TL3dis = TL2dis*TR_dis;             % between Coupler C and Coupler B
TL4dis = TL3dis*TR_dis;             % between Coupler B and Coupler A
TL5dis = T0+(HT.A(1).T-T0)*(1-eff); % Coupler A air outlet

% Use temperature points to compute mass flow rates at the PTES side
nT  = 100;
TvC_dis = linspace(TL2dis,TL3dis,nT)';
TvB_dis = linspace(TL3dis,TL4dis,nT)';
TvA_dis = linspace(TL4dis,TL5dis,nT)';
pv      = pmax_LA.*ones(size(TvC_dis));
CpvC    = RPN('PT_INPUTS',pv,TvC_dis,'C',gas1);
CpvB    = RPN('PT_INPUTS',pv,TvB_dis,'C',gas1);
CpvA    = RPN('PT_INPUTS',pv,TvA_dis,'C',gas1);
%figure(10); plot([TvC;TvB;TvA],[CpvC;CpvB;CpvA]); xlabel('Temperature, K'); ylabel('Cp, J/kg/K');
%legend(sprintf('%.1f bar',pmax_LA/1e5));
CpC = mean(CpvC);
CpB = mean(CpvB);
CpA = mean(CpvA);
CpN = RPN('PT_INPUTS',1e5,300,'C',gas2);
xc_dis  = mLA*CpC/CpN;
xd_dis  = mLA*CpB/CpN;
xe_dis  = mLA*CpA/CpN;
xb_dis  = xc_dis + xd_dis;
xa_dis  = xb_dis + xe_dis;


% Initial guess of charge conditions (PTES side).
% Turbine outlet (P1)
gas2.state(iL,iP1).p    = pbot_PT;
gas2.state(iL,iP1).T    = HT.A(2).T;
gas2.state(iL,iP1).mdot = xa_dis;
[gas2] = update(gas2,[iL,iP1],1);
% Regenerator A inlet (cold, high p side - P14)
iP14 = 14;
gas2.state(iL,iP14).p    = pbot_PT*PRC_dis*(1-ploss);
gas2.state(iL,iP14).T    = T0;
gas2.state(iL,iP14).mdot = xa_dis;
[gas2] = update(gas2,[iL,iP14],1);
% Regenerator B inlet (cold, high p side - P11)
iP11 = 11;
gas2.state(iL,iP11).p    = pbot_PT*PRC_dis/(1-ploss);
gas2.state(iL,iP11).T    = TL4dis;
gas2.state(iL,iP11).mdot = xb_dis;
[gas2] = update(gas2,[iL,iP11],1);
% Regenerator C inlet (cold, low p side - P9)
iP9 = 9;
gas2.state(iL,iP9).p    = pbot_PT*PRC_dis/(1-ploss)^2;
gas2.state(iL,iP9).T    = TL3dis;
gas2.state(iL,iP9).mdot = xc_dis;
[gas2] = update(gas2,[iL,iP9],1);


% Set matrix of temperature and pressure points to test convergence
C_0 = [[gas1.state(iL,:).T],[gas1.state(iL,:).p],[gas2.state(iL,:).T],[gas2.state(iL,:).p]];
max_iter=100;
for counter=1:max_iter
    
    fprintf(1,['Discharging CCES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%% LAES - PUMP %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%    
    % PUMP COMPRESSION (L1-->L2)
    p_aim = pmax_LA;
    [DCMP(1),gas1,iG1] = compexp_func(DCMP(1),iL,gas1,iG1,'Paim',p_aim,1);
    
    
    %%%%%%%%%%%%%%%%
    %%%%% PTES %%%%%
    %%%%%%%%%%%%%%%%
    
    % REGENERATOR A (neon-neon) (P1-->P2) & (P14-->P15)
    [HX(8),gas2,iP2,gas2X,iP15] = hex_func(HX(8),iL,gas2,iP1,gas2,iP14,0,0);
    gas2.state(iL,iP15) = gas2X.state(iL,iP15);
    gas2.stage(iL,iP14) = gas2X.stage(iL,iP14);
    
    % REJECT HEAT (external HEX) (P2-->P3)
    air.state(iL,iA).T = T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
    if gas2.state(iL,iP2).T >= T0
        [HX(9), gas2, iP3, air, iA] = hex_func(HX(9),iL,gas2,iP2,air,iA,1,0.5);
    else
        [HX(9), air, iA, gas2, iP3] = hex_func(HX(9),iL,air,iA,gas2,iP2,2,2.0);
    end
    [DFAN(iFAN),air,iA] = compexp_func (DFAN(iFAN),iL,air,iA,'Paim',p0,1); %#ok<*SAGROW>
    iFAN=iFAN+1; iA=iA+1;
    
    % SEPARATE (P3-->P4 & P18)
    iP18  = 18;    % index of secondary stream born from split
    x_aim = xe_dis/xa_dis; % mass flow rate fraction of secondary stream
    [gas2,iP4] = split_stream(gas2,iL,iP3,iP18,x_aim);
    
    % REGENERATOR B (neon-neon) (P4-->P5) & (P11-->P12)
    [HX(10),gas2,iP5,gas2X,iP12] = hex_func(HX(10),iL,gas2,iP4,gas2,iP11,0,0);
    gas2.state(iL,iP12) = gas2X.state(iL,iP12);
    gas2.stage(iL,iP11) = gas2X.stage(iL,iP11);
    
    % SEPARATE (P5-->P6 & P22)
    iP22  = 22;    % index of secondary stream born from split
    x_aim = xd_dis/xb_dis; % mass flow rate fraction of secondary stream
    [gas2,iP6] = split_stream(gas2,iL,iP5,iP22,x_aim);
    
    % REGENERATOR C (neon-neon) (P6-->P7) & (P9-->P10)
    [HX(11),gas2,iP7,gas2X,iP10] = hex_func(HX(11),iL,gas2,iP6,gas2,iP9,0,0);
    gas2.state(iL,iP10) = gas2X.state(iL,iP10);
    gas2.stage(iL,iP9) = gas2X.stage(iL,iP9);
    
    % COUPLER C (neon-air) (P7->P8) (L2-->L3)
    [HX(6),gas2,iP8,gas1,iG1] = hex_func(HX(6),iL,gas2,iP7,gas1,iG1,0,0);
    
    % COMPRESSION (P8-->P9)
    p_aim = gas2.state(iL,iP8).p*PRC_dis;
    [DCMP(2),gas2,iP9] = compexp_func (DCMP(2),iL,gas2,iP8,'Paim',p_aim,1);
    
    % COUPLER B (neon-air) (P22->P23) (L3-->L4)
    [HX(5),gas2,iP23,gas1,iG1] = hex_func(HX(5),iL,gas2,iP22,gas1,iG1,0,0);
    
    % COMPRESSION (P23-->P24)
    p_aim = gas2.state(iL,iP10).p;
    [DCMP(3),gas2,iP24] = compexp_func (DCMP(3),iL,gas2,iP23,'Paim',p_aim,1);
    
    % MIX (P10 & P24-->P11)
    [MIX(3),gas2,iP11] = mix_streams2(MIX(3),gas2,[iL,iP10],[iL,iP24]);
    iP25 = iP24+1;
    gas2.state(iL,iP25) = gas2.state(iL,iP11);
    gas2.state(iL,iP25).mdot = 0;
    gas2.stage(iL,iP25).type = 'end';
    
    % COUPLER A (neon-air) (P18->P19) (L4-->L5)
    [HX(4),gas2,iP19,gas1,iG1] = hex_func(HX(4),iL,gas2,iP18,gas1,iG1,0,0);
    
    % COMPRESSION (P19-->P20)
    p_aim = gas2.state(iL,iP12).p;
    [DCMP(4),gas2,iP20] = compexp_func (DCMP(4),iL,gas2,iP19,'Paim',p_aim,1);
    
    % MIX (P12 & P20-->P13)
    [MIX(4),gas2,iP13] = mix_streams2(MIX(4),gas2,[iL,iP12],[iL,iP20]);
    iP21 = iP20+1;
    gas2.state(iL,iP21) = gas2.state(iL,iP13);
    gas2.state(iL,iP21).mdot = 0;
    gas2.stage(iL,iP21).type = 'end';
    
    % AUXILIARY COMPRESSION (P13-->P14)
    p_aim = pbot_PT*PRH_dis;
    %p_aim = gas2.state(iL,iP13).p*PRaux_dis;
    [DCMP(5),gas2,iP14] = compexp_func (DCMP(5),iL,gas2,iP13,'Paim',p_aim,1);
    
    % HEAT (P15-->P16)
    % Set fluidH
    fluidH.state(iL,iH).T = HT.B(iL).T; fluidH.state(iL,iH).p = HT.B(iL).p;
    [fluidH] = update(fluidH,[iL,iH],1);
    THmin = HT.A(1).T; Taim = THmin;
    % Run HEX
    %[HX(7),fluidH,iH,gas2,iP16] = hex_func(HX(7),iL,fluidH,iH,gas2,iP15,4,THmin);
    [HX(7),fluidH,iH,gas2,iP16] = hex_func(HX(7),iL,fluidH,iH,gas2,iP15,2,1.0);
    % Run Pump
    [DPMP(iPMP),fluidH,iH] = compexp_func (DPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
    iH=iH+1; iPMP=iPMP+1;
    
    % EXPANSION (P16-->P17=P1)
    p_aim = gas2.state(iL,iP1).p;
    [DEXP(3),gas2,iP17] = compexp_func(DEXP(3),iL,gas2,iP16,'Paim',p_aim,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%% CLOSE LAES %%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    % HEAT (L5-->L6)
    % Set fluidM
    fluidM.state(iL,iM).T = MT.B(iL).T; fluidM.state(iL,iM).p = MT.B(iL).p;
    [fluidM] = update(fluidM,[iL,iM],1);
    % Run HEX
    [HX(2),fluidM,iM,gas1,iG1] = hex_func(HX(2),iL,fluidM,iM,gas1,iG1,2,1.0);
    % Run Pump
    [DPMP(iPMP),fluidM,iM] = compexp_func(DPMP(iPMP),iL,fluidM,iM,'Paim',fluidM.state(iL,1).p,1); %#ok<*SAGROW>
    iM=iM+1;iPMP=iPMP+1;
    
    % EXPANSION (L6-->L7)
    p_aim = p0*(gas1.state(iL,iG1).p/p0)^(1/2);
    [DEXP(1),gas1,iG1] = compexp_func (DEXP(1),iL,gas1,iG1,'Paim',p_aim, 1) ;
    
    % REJECT HEAT (external HEX) (L7-->L8)
    air.state(iL,iA).T = T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
    [HX(3), gas1, iG1, air, iA] = hex_func(HX(3),iL,gas1,iG1,air,iA,1,0.5);
    [DFAN(iFAN),air,iA] = compexp_func (DFAN(iFAN),iL,air,iA,'Paim',p0,1);
    iFAN=iFAN+1; iA=iA+1;
    
    % HEAT (L8-->L9)
    % Set fluidM
    fluidM.state(iL,iM).T = MT.B(iL).T; fluidM.state(iL,iM).p = MT.B(iL).p;
    [fluidM] = update(fluidM,[iL,iM],1);
    % Run HEX
    [HX(1),fluidM,iM,gas1,iG1] = hex_func(HX(1),iL,fluidM,iM,gas1,iG1,2,1.0);
    % Run Pump
    [DPMP(iPMP),fluidM,iM] = compexp_func(DPMP(iPMP),iL,fluidM,iM,'Paim',fluidM.state(iL,1).p,1); %#ok<*SAGROW>
    iM=iM+1;iPMP=iPMP+1;
    
    % EXPANSION (L9-->L10)
    p_aim = p0;
    [DEXP(2),gas1,iG1] = compexp_func (DEXP(2),iL,gas1,iG1,'Paim',p_aim, 1) ;
    
    % Determine convergence and proceed
    C = [[gas1.state(iL,:).T],[gas1.state(iL,:).p],[gas2.state(iL,:).T],[gas2.state(iL,:).p]];
    err = abs((C(C~=0) - C_0(C~=0))./C(C~=0))*100; % percentage error
    fprintf('Convergence error = %f %%\n',max(err))
    convergence = all(err < 1e-3);
    
    if convergence || counter==max_iter % is charge cycle converged?

        % Close working fluid streams
        gas1.stage(iL,iG1).type = 'end';
        gas1 = count_Nstg(gas1);
        gas2.stage(iL,iP17).type = 'end';
        gas2 = count_Nstg(gas2);
        iLA_out = 1;
        iLA_in  = gas1.Nstg(iL)+1;
        
        % Close air (heat rejection) streams
        iA_out = 1:3:(iA-1); iA_in  = iA_out + 2;
        for i=iA_in, air.stage(iL,i).type = 'end'; end
        air = count_Nstg(air);
        
        % Close storage fluid streams
        iH_out = 1:3:(iH-1); iH_in  = iH_out + 2;
        iM_out = 1:3:(iM-1); iM_in  = iM_out + 2;
        for i=iH_in, fluidH.stage(iL,i).type = 'end'; end
        for i=iM_in, fluidM.stage(iL,i).type = 'end'; end
        fluidH = count_Nstg(fluidH);
        fluidM = count_Nstg(fluidM);
                
        % Uncomment these lines to print states
        %{
        print_states(gas1,iL,1:gas1.Nstg(iL)+1,Load);
        print_states(gas2,iL,1:gas2.Nstg(iL)+1,Load);
        print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
        print_states(fluidM,iL,1:fluidM.Nstg(iL)+1,Load);
        print_states(air,iL,1:air.Nstg(iL)+1,Load);
        %keyboard
        %}
        
        % Uncomment these lines to prints HEXs
        %{
        plot_hex(HX(1) ,iL,41,'K',true,'LAES hot HEX 1')
        plot_hex(HX(2) ,iL,42,'K',true,'LAES hot HEX 2')
        plot_hex(HX(3) ,iL,43,'K',true,'LAES heat rejection')
        plot_hex(HX(4) ,iL,44,'K',true,'Coupler A')
        plot_hex(HX(5) ,iL,45,'K',true,'Coupler B')
        plot_hex(HX(6), iL,46,'K',true,'Coupler C')
        plot_hex(HX(7) ,iL,47,'K',true,'PTES hot HEX')
        plot_hex(HX(8) ,iL,48,'K',true,'PTES Regen A')
        plot_hex(HX(9) ,iL,49,'K',true,'PTES heat rejection')
        plot_hex(HX(10),iL,50,'K',true,'PTES Regen B')
        plot_hex(HX(11),iL,51,'K',true,'PTES Regen C')
        %keyboard
        %}
        
        % Exit loop
        break
        
    else
        gas2.state(iL,1) = gas2.state(iL,iP17);
        C_0 = C;
        iG1=1; iH=1; iM=1; iA=1; iFAN=1; iPMP=1;
        
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
[MdotLA] = total_mdot(gas1, iL, iLA_out);
t_disLA  = LAT.B(iL).M/MdotLA;
Load.time(iL) = min([Load.time(iL),t_disH,t_disM,t_disLA])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
HT = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
% Cold tanks
MT = run_tanks(MT,iL,fluidM,iM_out,iM_in,Load,T0);
% Atmospheric tanks
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);
% Liquid air tanks
LAT = run_tanks(LAT,iL,gas1,iLA_out,iLA_in,Load,T0);


% PLOT DISCHARGE CYCLE
%PLOT_CYCLE
%keyboard