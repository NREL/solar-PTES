%%% CCES PLANT LAYOUT DURING CHARGE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAES SIDE
%  ---------
%  L1: Compressor inlet (atmosphere)    - CCMP(1)
%  L2: Hot HEX inlet                    - HX(1)
%  L3: Compressor inlet                 - CCMP(2)
%  L4: Hot HEX inlet                    - HX(2)
%  L5: Rejection unit inlet             - HX(3)
%  L6: Coupler A inlet (vs. P19)        - HX(4)
%  L7: Coupler B inlet (vs. P23)        - HX(5)
%  L8: Coupler C inlet (vs. P11)        - HX(6)
%  L9: Expander inlet                   - CEXP(1)
% L10: Liquid air tank inlet            - ...
%
%  PTES SIDE
%  ---------
%  P1: Compressor inlet                         - xa   - CCMP(3)
%  P2: Hot HEX inlet                            - xa   - HX(7)
%  P3: Regen A inlet (vs. P16)                  - xa   - HX(8)
%  P4: Auxiliary compressor inlet               - xa   - CCMP(4)
%  P5: Rejection unit inlet                     - xa   - HX(9)
%  P6: Separator inlet (makes P7 & P18)         - xa   - ...
%  P7: Regen B inlet (vs. P14)                  - xb   - HX(10)
%  P8: Separator inlet (makes P9 & P22)         - xb   - ...
%  P9: Regen C inlet (vs. P12)                  - xc   - HX(11)
% P10: Expander inlet                           - xc   - CEXP(2)
% P11: Coupler C inlet (vs. L8)                 - xc   - HX(6)
% P12: Regen C inlet (vs. P9)                   - xc   - HX(11)
% P13: Mixer inlet (mixes with P24, makes P14)  - xc   - MIX(1)
% P14: Regen B inlet (vs. P7)                   - xb   - HX(10)
% P15: Mixer inlet (mixes with P20, makes P16)  - xb   - MIX(2)
% P16: Regen A inlet (vs. P3)                   - xa   - HX(8)
% P17: Closure of main cycle stream = P1        - xa   - end
% P18: Expander inlet                           - xe   - CEXP(4)
% P19: Coupler A inlet (vs. L6)                 - xe   - HX(4)
% P20: Mixer inlet (mixes with P15, makes P16)  - xe   - MIX(2)
% P21: End of stream xe                         - 0    - end
% P22: Expander inlet                           - xd   - CEXP(3)
% P23: Coupler B inlet (vs. L7)                 - xd   - HX(5)
% P24: Mixer inlet (mixes with P13, makes P14)  - xd   - MIX(1)
% P25: End of stream xd                         - 0    - end
%
% where:
% P6:   xa = xb + xe
% P8:   xb = xc + xd
% P13:  xc + xd = xb
% P15:  xb + xe = xa
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

% Initial guess of charge conditions (LAES side)
% Compressor inlet
gas1.state(iL,1).p = p0;
gas1.state(iL,1).T = T0;
mLA = Load.mdot(iL);
gas1.state(iL,1).mdot = mLA;
[gas1] = update(gas1,[iL,1],1);

% Obtain mass flow rates at PTES side. This requires computation of
% specific heat capacity of supercritical air between points L6 and L9.
% Assume 'logspace distribution' of temperatures between liquid air
% temperature and ambient temperature (3 couplers, 4 temperature points).
Tpnts = logspace(log10(TLA_obj),log10(T0),4); 
TL9chg = Tpnts(1);
TL8chg = Tpnts(2);
TL7chg = Tpnts(3);
TL6chg = Tpnts(4);
nT  = 100;
TvC = linspace(TL9chg,TL8chg,nT)';
TvB = linspace(TL8chg,TL7chg,nT)';
TvA = linspace(TL7chg,TL6chg,nT)';
pv  = pmax_LA.*ones(size(TvC));
Cpv1 = RPN('PT_INPUTS',pv,TvC,'C',gas1);
Cpv2 = RPN('PT_INPUTS',pv,TvB,'C',gas1);
Cpv3 = RPN('PT_INPUTS',pv,TvA,'C',gas1);
%figure(10); plot([Tv1;Tv2;Tv3],[Cpv1;Cpv2;Cpv3]); xlabel('Temperature, K'); ylabel('Cp, J/kg/K');
%legend(sprintf('%.1f bar',pmax_LA/1e5));
Cp1 = mean(Cpv1);
Cp2 = mean(Cpv2);
Cp3 = mean(Cpv3);
CpN = RPN('PT_INPUTS',pbot_PT,300,'C',gas2);
xc  = mLA*Cp1/CpN;
xd  = mLA*Cp2/CpN;
xe  = mLA*Cp3/CpN;
xb  = xc + xd;
xa  = xb + xe;

% Set charge pressure ratio
if strcmp(gas2.name,'Neon')
    gamma = 1.667;
    phi_exp  = eta*(gamma-1)/gamma;
    phi_comp = (gamma-1)/(gamma*eta);
else
    error('not implemented')
end
TR_exp  = Tpnts(2)/Tpnts(1);
PRC_chg = (TR_exp/eff)^(1/phi_exp)/(1-ploss)^2; % pressure ratio PTES 'cold' line

% Initial guess of charge conditions (PTES side).
% Compressor inlet (P1)
gas2.state(iL,iP1).p    = pbot_PT;
gas2.state(iL,iP1).T    = HT.A(1).T;
gas2.state(iL,iP1).mdot = xa;
[gas2] = update(gas2,[iL,iP1],1);
% Regenerator A inlet (cold, low p side - P16)
iP16 = 16;
gas2.state(iL,iP16).p    = pbot_PT/(1-ploss);
gas2.state(iL,iP16).T    = T0;
gas2.state(iL,iP16).mdot = xa;
[gas2] = update(gas2,[iL,iP16],1);
% Regenerator B inlet (cold, low p side - P14)
iP14 = 14;
gas2.state(iL,iP14).p    = pbot_PT/(1-ploss)^2;
gas2.state(iL,iP14).T    = TL7chg;
gas2.state(iL,iP14).mdot = xb;
[gas2] = update(gas2,[iL,iP14],1);
% Regenerator C inlet (cold, low p side - P12)
iP12 = 12;
gas2.state(iL,iP12).p    = pbot_PT/(1-ploss)^3;
gas2.state(iL,iP12).T    = TL8chg;
gas2.state(iL,iP12).mdot = xc;
[gas2] = update(gas2,[iL,iP12],1);


% Set matrix of temperature and pressure points to test convergence
C_0 = [[gas1.state(iL,:).T],[gas1.state(iL,:).p],[gas2.state(iL,:).T],[gas2.state(iL,:).p]];
max_iter=100;
for counter=1:max_iter
    
    fprintf(1,['Charging CCES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    
    %%%%%%%%%%%%%%%%
    %%%%% LAES %%%%%
    %%%%%%%%%%%%%%%%    
    % FIRST COMPRESSION (L1-->L2)
    p_aim = (pmax_LA/p0)^(1/2)*p0*1.15;
    [CCMP(1),gas1,iG1] = compexp_func(CCMP(1),iL,gas1,iG1,'Paim',p_aim,1);
    
    % COOL (L2-->L3)
    % Set fluidM
    fluidM.state(iL,iM).T = MT.A(iL).T; fluidM.state(iL,iM).p = MT.A(iL).p;
    [fluidM] = update(fluidM,[iL,iM],1);
    % Run HEX
    [HX(1),gas1,iG1,fluidM,iM] = hex_func(HX(1),iL,gas1,iG1,fluidM,iM,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidM,iM] = compexp_func(CPMP(iPMP),iL,fluidM,iM,'Paim',fluidM.state(iL,1).p,1); %#ok<*SAGROW>
    iM=iM+1;iPMP=iPMP+1;
    
    % SECOND COMPRESSION (L3-->L4)
    p_aim = pmax_LA;
    [CCMP(2),gas1,iG1] = compexp_func(CCMP(2),iL,gas1,iG1,'Paim',p_aim,1);
    
    % COOL (L4-->L5)
    % Set fluidM
    fluidM.state(iL,iM).T = MT.A(iL).T; fluidM.state(iL,iM).p = MT.A(iL).p;
    [fluidM] = update(fluidM,[iL,iM],1);
    % Run HEX
    [HX(2),gas1,iG1,fluidM,iM] = hex_func(HX(2),iL,gas1,iG1,fluidM,iM,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidM,iM] = compexp_func(CPMP(iPMP),iL,fluidM,iM,'Paim',fluidM.state(iL,1).p,1); %#ok<*SAGROW>
    iM=iM+1;iPMP=iPMP+1;
    
    % REJECT HEAT (external HEX) (L5-->L6)
    air.state(iL,iA).T = T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
    [HX(3), gas1, iG1, air, iA] = hex_func(HX(3),iL,gas1,iG1,air,iA,1,0.5);
    [CFAN(iFAN),air,iA] = compexp_func (CFAN(iFAN),iL,air,iA,'Paim',p0,1);
    iFAN=iFAN+1; iA=iA+1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PTES - COMP & REGENS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % FIRST COMPRESSION (P1-->P2)
    T_aim = Tmax;
    [CCMP(3),gas2,iP2] = compexp_func(CCMP(3),iL,gas2,iP1,'Taim',T_aim,1);
    PRH_chg = gas2.state(iL,iP2).p/gas2.state(iL,iP1).p; % pressure ratio PTES 'hot' line
    
    % COOL (P2-->P3)
    % Set fluidH
    fluidH.state(iL,iH).T = HT.A(iL).T; fluidH.state(iL,iH).p = HT.A(iL).p;
    [fluidH] = update(fluidH,[iL,iH],1);
    % Run HEX
    [HX(7),gas2,iP3,fluidH,iH] = hex_func(HX(7),iL,gas2,iP2,fluidH,iH,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidH,iH] = compexp_func(CPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
    iH=iH+1;iPMP=iPMP+1;
    
    % REGENERATOR A (neon-neon) (P3-->P4) & (P16-->P17=P1)
    [HX(8),gas2,iP4,gas2X,iP17] = hex_func(HX(8),iL,gas2,iP3,gas2,iP16,0,0);
    gas2.state(iL,iP17)   = gas2X.state(iL,iP17);
    gas2.stage(iL,iP17-1) = gas2X.stage(iL,iP17-1);
    
    % SECOND COMPRESSION (P4-->P5)
    p_aim = PRC_chg*gas2.state(iL,iP1).p;
    [CCMP(4),gas2,iP5] = compexp_func (CCMP(4),iL,gas2,iP4,'Paim',p_aim, 1) ;
    
    % REJECT HEAT (external HEX) (P5-->P6)
    air.state(iL,iA).T = T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
    [HX(9), gas2, iP6, air, iA] = hex_func(HX(9),iL,gas2,iP5,air,iA,1,0.5);
    [CFAN(iFAN),air,iA] = compexp_func (CFAN(iFAN),iL,air,iA,'Paim',p0,1);
    iFAN=iFAN+1; iA=iA+1;
        
    % SEPARATE (P6-->P7 & P18)
    iP18  = 18;    % index of secondary stream born from split
    x_aim = xe/xa; % mass flow rate fraction of secondary stream
    [gas2,iP7] = split_stream(gas2,iL,iP6,iP18,x_aim);

    % REGENERATOR B (neon-neon) (P7-->P8) & (P14-->P15)
    [HX(10),gas2,iP8,gas2X,iP15] = hex_func(HX(10),iL,gas2,iP7,gas2,iP14,0,0);
    gas2.state(iL,iP15) = gas2X.state(iL,iP15);
    gas2.stage(iL,iP14) = gas2X.stage(iL,iP14);
    
    % SEPARATE (P8-->P9 & P22)
    iP22  = 22;    % index of secondary stream born from split
    x_aim = xd/xb; % mass flow rate fraction of secondary stream
    [gas2,iP9] = split_stream(gas2,iL,iP8,iP22,x_aim);
    
    % REGENERATOR C (neon-neon) (P9-->P10) & (P12-->P13)
    [HX(11),gas2,iP10,gas2X,iP13] = hex_func(HX(11),iL,gas2,iP9,gas2,iP12,0,0);
    gas2.state(iL,iP13) = gas2X.state(iL,iP13);
    gas2.stage(iL,iP12) = gas2X.stage(iL,iP12);
    
    %print_states(gas1,iL,1:12,Load);
    %print_states(gas2,iL,1:26,Load);
    %keyboard
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PTES - EXP & COUPLERS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % EXPANSION (P18-->P19)
    p_aim = pbot_PT/(1-ploss)^2;
    [CEXP(4),gas2,iP19] = compexp_func (CEXP(4),iL,gas2,iP18,'Paim',p_aim,1);
    
    % COUPLER A (air-neon) (L6-->L7) (P19->P20)
    [HX(4),gas1,iG1,gas2,iP20] = hex_func(HX(4),iL,gas1,iG1,gas2,iP19,0,0);
    
    % EXPANSION (P22-->P23)
    p_aim = pbot_PT/(1-ploss)^3;
    [CEXP(3),gas2,iP23] = compexp_func (CEXP(3),iL,gas2,iP22,'Paim',p_aim,1);
    
    % COUPLER B (air-neon) (L7-->L8) (P23->P24)
    [HX(5),gas1,iG1,gas2,iP24] = hex_func(HX(5),iL,gas1,iG1,gas2,iP23,0,0);
    
    % EXPANSION (P10-->P11)
    p_aim = pbot_PT/(1-ploss)^4;
    [CEXP(2),gas2,iP11] = compexp_func (CEXP(2),iL,gas2,iP10,'Paim',p_aim,1);
    
    % COUPLER C (air-neon) (L8-->L9) (P11->P12)
    [HX(6),gas1,iG1,gas2,iP12] = hex_func(HX(6),iL,gas1,iG1,gas2,iP11,0,0);
    
    %print_states(gas1,iL,1:12,Load);
    %print_states(gas2,iL,1:26,Load);
    %keyboard
   
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%% CLOSE LAES %%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    % LIQUID AIR EXPANSION (L7-->L8)
    pLA = max([RPN('QT_INPUTS',0.0,gas1.state(iL,iG1).T,'P',gas1)*1.2,p0]);
    [CEXP(1),gas1,iG1] = compexp_func (CEXP(1),iL,gas1,iG1,'Paim',pLA, 1) ;
    TLA = gas1.state(iL,iG1).T;
    hLA = gas1.state(iL,iG1).h;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% CLOSE PTES  - MIX %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % MIX (P13 & P24-->P14)
    [MIX(1),gas2,iP14] = mix_streams2(MIX(1),gas2,[iL,iP13],[iL,iP24]);
    iP25 = iP24+1;
    gas2.state(iL,iP25) = gas2.state(iL,iP14);
    gas2.state(iL,iP25).mdot = 0;
    gas2.stage(iL,iP25).type = 'end';
    
    % MIX (P15 & P20-->P16)
    [MIX(2),gas2,iP16] = mix_streams2(MIX(2),gas2,[iL,iP15],[iL,iP20]);
    iP21 = iP20+1;
    gas2.state(iL,iP21) = gas2.state(iL,iP16);
    gas2.state(iL,iP21).mdot = 0;
    gas2.stage(iL,iP21).type = 'end';
    
    % Determine convergence and proceed
    C = [[gas1.state(iL,:).T],[gas1.state(iL,:).p],[gas2.state(iL,:).T],[gas2.state(iL,:).p]];
    err = abs((C(C~=0) - C_0(C~=0))./C(C~=0))*100; % percentage error
    fprintf('Convergence error = %f %%\n',max(err))
    convergence = all(err < 1e-3);
    %print_states(gas1,iL,1:12,Load);
    %print_states(gas2,iL,1:26,Load);
    %keyboard
    
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
        plot_hex(HX(1) ,1,41,'K',true,'LAES hot HEX 1')
        plot_hex(HX(2) ,1,42,'K',true,'LAES hot HEX 2')
        plot_hex(HX(3) ,1,43,'K',true,'LAES heat rejection')
        plot_hex(HX(4) ,1,44,'K',true,'Coupler A')
        plot_hex(HX(5) ,1,45,'K',true,'Coupler B')
        plot_hex(HX(6), 1,46,'K',true,'Coupler C')
        plot_hex(HX(7) ,1,47,'K',true,'PTES hot HEX')
        plot_hex(HX(8) ,1,48,'K',true,'PTES Regen A')
        plot_hex(HX(9) ,1,49,'K',true,'PTES heat rejection')
        plot_hex(HX(10),1,50,'K',true,'PTES Regen B')
        plot_hex(HX(11),1,51,'K',true,'PTES Regen C')
        keyboard
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

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
HT = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
% Cold tanks
MT = run_tanks(MT,iL,fluidM,iM_out,iM_in,Load,T0);
% Atmospheric tanks
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);
% Liquid air tanks
for iT = 1:(Load.num+1)
    LAT.B(iT).p = pLA;
    LAT.B(iT).T = TLA;
    LAT.B(iT).h = hLA;
end
LAT = run_tanks(LAT,iL,gas1,iLA_out,iLA_in,Load,T0);

% PLOT CHARGE CYCLE
%PLOT_CYCLE
%keyboard