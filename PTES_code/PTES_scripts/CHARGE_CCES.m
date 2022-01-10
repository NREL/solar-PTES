%%% CCES PLANT LAYOUT %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAES SIDE
%  ---------
%  L1: Compressor inlet (atmosphere)     - CCMP(1)
%  L2: Hot HEX inlet                     - HX(1)
%  L3: Compressor inlet                  - CCMP(2)
%  L4: Hot HEX inlet                     - HX(2)
%  L5: Rejection unit inlet              - HX(3)
%  L6: High Coupler inlet (air side)     - HX(4)
%  L7: Cryo Coupler inlet (air side)     - HX(5)
%  L8: Expander inlet                    - CEXP(1)
%  L9: Liquid air tank inlet
%
%  PTES SIDE
%  ---------
%  P1: Compressor inlet                          - xa   - CCMP(3)
%  P2: Hot HEX inlet                             - xa   - HX(6)
%  P3: Auxiliary compressor inlet                - xa   - CCMP(4)
%  P4: Rejection unit inlet                      - xa   - HX(7)
%  P5: Separator inlet (after P4, sep. P6 & P12) - xa   - 
%  P6: Regenerator inlet (high p side)           - xb   - HX(8)
%  P7: Expander inlet (low temp exp)             - xb   - CEXP(2)
%  P8: Cryo Coupler inlet (neon side)            - xb
%  P9: Separator inlet (after P8, sep. P10 & P14)- xb
% P10: High Coupler inlet (neon side)            - xc
% P11: Mixer inlet (after P10, mix P16, res. P1) - xc
% P12: Expander inlet (ambient exp)              - xc   - CEXP(3)
% P13: Mixer inlet (after P12, mix P14, res. P15)- xc
% P14: Mixer inlet (after P9, mix P13, res. P15) - xd
% P15: Regenerator inlet (low p side)            - xb
% P16: Mixer inlet (after P15, mix P11, RES. P1) - xb
% 
% where:
% P5:   xa = xb + xc
% P9:   xb = xc + xd
% P15:  xb = xc + xd
% P1:   xa = xb + xc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set stage indices
iG1 = 1;  % keeps track of the gas stage number (LAES)
iP1 = 1;  % keeps track of the gas stage number (PTES)
iH  = 1;  % keeps track of the Hot fluid stream number
iM  = 1;  % keeps track of the Medium fluid stream number
iA  = 1;  % keeps track of the Air (heat rejection) stream number
iPMP = 1 ; % Keeps track of which pump is being used
iFAN = 1 ; % Keeps track of which fan is being used

% Initial guess of charge conditions (LAES side)
% Compressor inlet
gas1.state(iL,1).p = p0;
gas1.state(iL,1).T = T0;
mLA = Load.mdot(iL);
gas1.state(iL,1).mdot = mLA;
[gas1] = update(gas1,[iL,1],1);

% Obtain mass flow rate at PTES side. This requires computation of specific
% heat capacity of supercritical air between points L5 and L7.
Tm  = 150; nT = 100;
Tv1 = linspace(78,Tm,nT)';
Tv2 = linspace(Tm,T0,nT)';
pv  = pmax_LA.*ones(size(Tv1));
Cpv1 = RPN('PT_INPUTS',pv,Tv1,'C',gas1);
Cpv2 = RPN('PT_INPUTS',pv,Tv2,'C',gas1);
%figure(10); plot([Tv1;Tv2],[Cpv1;Cpv2]); xlabel('Temperature, K'); ylabel('Cp, J/kg/K');
%legend(sprintf('%.1f bar',pmax_LA/1e5));
Cp1 = mean(Cpv1);
Cp2 = mean(Cpv2);
CpN = RPN('PT_INPUTS',1e5,300,'C',gas2);
xb  = mLA*Cp1/CpN;
xc  = mLA*Cp2/CpN;
xa  = xb + xc;
xd  = xb - xc;
mN  = xa;

% Initial guess of charge conditions (PTES side).
% Compressor inlet (P1)
gas2.state(iL,iP1).p    = p0;
gas2.state(iL,iP1).T    = T0;
gas2.state(iL,iP1).mdot = mN;
[gas2] = update(gas2,[iL,iP1],1);
% Regenerator inlet (cold, low p side - P15)
iP15 = 15;
gas2.state(iL,iP15).p    = p0/(1-ploss);
gas2.state(iL,iP15).T    = Tm;
gas2.state(iL,iP15).mdot = xb;
[gas2] = update(gas2,[iL,iP15],1);
% High Coupler inlet (neon side - P10)
iP10 = 10;
gas2.state(iL,iP10).p    = p0/(1-ploss);
gas2.state(iL,iP10).T    = Tm;
gas2.state(iL,iP10).mdot = xc;
[gas2] = update(gas2,[iL,iP10],1);

% Set matrix of temperature and pressure points to test convergence
C_0 = [[gas1.state(iL,:).T],[gas1.state(iL,:).p],[gas2.state(iL,:).T],[gas2.state(iL,:).p]];
max_iter=100;
for counter=1:max_iter
    
    fprintf(1,['Charging CCES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    
    %%%%%%%%%%%%%%%%
    %%%%% LAES %%%%%
    %%%%%%%%%%%%%%%%    
    % FIRST COMPRESSION (L1-->L2)
    T_aim = Tmax;
    [CCMP(1),gas1,iG1] = compexp_func(CCMP(1),iL,gas1,iG1,'Taim',T_aim,1);
    
    % COOL (L2-->L3)
    % Set fluidH
    fluidH.state(iL,iH).T = HT.A(iL).T; fluidH.state(iL,iH).p = HT.A(iL).p;
    [fluidH] = update(fluidH,[iL,iH],1);
    % Run HEX
    [HX(1),gas1,iG1,fluidH,iH] = hex_func(HX(1),iL,gas1,iG1,fluidH,iH,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidH,iH] = compexp_func(CPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1); %#ok<*SAGROW>
    iH=iH+1;iPMP=iPMP+1;
    
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
    
    %print_states(gas1,iL,1:10,Load);
    %keyboard
    
    %%%%%%%%%%%%%%%%
    %%%%% PTES %%%%%
    %%%%%%%%%%%%%%%%
    % FIRST COMPRESSION (P1-->P2)
    T_aim = Tmax;
    [CCMP(3),gas2,iP2] = compexp_func(CCMP(3),iL,gas2,iP1,'Taim',T_aim,1);
    
    % COOL (P2-->P3)
    % Set fluidH
    fluidH.state(iL,iH).T = HT.A(iL).T; fluidH.state(iL,iH).p = HT.A(iL).p;
    [fluidH] = update(fluidH,[iL,iH],1);
    % Run HEX
    [HX(6),gas2,iP3,fluidH,iH] = hex_func(HX(6),iL,gas2,iP2,fluidH,iH,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidH,iH] = compexp_func(CPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
    iH=iH+1;iPMP=iPMP+1;
    
    % SECOND COMPRESSION (P3-->P4)
    p_aim = gas2.state(iL,iP3).p*1.2;
    [CCMP(4),gas2,iP4] = compexp_func (CCMP(4),iL,gas2,iP3,'Paim',p_aim, 1) ;
    
    % REJECT HEAT (external HEX) (P4-->P5)
    air.state(iL,iA).T = T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
    [HX(7), gas2, iP5, air, iA] = hex_func(HX(7),iL,gas2,iP4,air,iA,1,0.5);
    [CFAN(iFAN),air,iA] = compexp_func (CFAN(iFAN),iL,air,iA,'Paim',p0,1);
    iFAN=iFAN+1; iA=iA+1;
    
    % SEPARATE (P5-->P6 & P12)
    iP12 = 12;    % index of secondary stream born from split
    x12  = xc/xa; % mass flow rate fraction of secondary stream
    [gas2,iP6] = split_stream(gas2,iL,iP5,iP12,x12);
    
    % REGENERATE (neon-neon) (P6-->P7)
    [HX(8),gas2,iP7,~,~] = hex_func(HX(8),iL,gas2,iP6,gas2,iP15,0,0);
    
    % EXPANSION (P7-->P8)
    p_aim = p0/(1-ploss)^2;
    [CEXP(2),gas2,iP8] = compexp_func (CEXP(2),iL,gas2,iP7,'Paim',p_aim,1);
    
    %print_states(gas2,iL,1:17,Load);
    %keyboard
    
    %%%%%%%%%%%%%%%%%%%%
    %%%%% COUPLERS %%%%%
    %%%%%%%%%%%%%%%%%%%%
    
    % HIGH COUPLER (air-neon) (L5-->L6) (P10->P11)
    [HX(4),gas1,iG1,gas2,iP11] = hex_func(HX(4),iL,gas1,iG1,gas2,iP10,0,0);
    
    % CRYO COUPLER (air-neon) (L6-->L7) (P8->P9)
    [HX(5),gas1,iG1,gas2,iP9] = hex_func(HX(5),iL,gas1,iG1,gas2,iP8,0,0);

    
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%% CLOSE LAES %%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    % LIQUID AIR EXPANSION (L7-->L8)
    p_aim = max([RPN('QT_INPUTS',0.0,gas1.state(iL,iG1).T,'P',gas1)*1.01,p0]);
    [CEXP(1),gas1,iG1] = compexp_func (CEXP(1),iL,gas1,iG1,'Paim',p_aim, 1) ;
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%% CLOSE PTES %%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    % SEPARATE (P9-->P10 & P14)
    iP14 = 14;    % index of secondary stream born from split
    x14  = xd/xb; % mass flow rate fraction of secondary stream
    [gas2,iP10] = split_stream(gas2,iL,iP9,iP14,x14);
    
    % EXPANSION (P12-->P13)
    p_aim = p0/(1-ploss);
    [CEXP(3),gas2,iP13] = compexp_func (CEXP(3),iL,gas2,iP12,'Paim',p_aim,1);
    
    % MIX (P13 & P14-->P15)
    [MIX(1),gas2,iP15] = mix_streams2(MIX(1),gas2,[iL,iP14],[iL,iP13]);
    
    % REGENERATE (neon-neon) (P15-->P16)
    [HX(8),~,~,gas2,iP16] = hex_func(HX(8),iL,gas2,iP6,gas2,iP15,0,0);
    
    % MIX (P11 & P16-->P1)
    [MIX(2),gas2,iP17] = mix_streams2(MIX(2),gas2,[iL,iP16],[iL,iP11]);
    
    % Determine convergence and proceed
    C = [[gas1.state(iL,:).T],[gas1.state(iL,:).p],[gas2.state(iL,:).T],[gas2.state(iL,:).p]];
    err = abs((C(C~=0) - C_0(C~=0))./C(C~=0))*100; % percentage error
    fprintf('Convergence error = %f %%\n',max(err))
    convergence = all(err < 1e-1);
    %keyboard
    
    if convergence || counter==max_iter % is charge cycle converged?

        % Close working fluid streams
        gas1.stage(iL,iG1).type = 'end';
        gas1 = count_Nstg(gas1);
        gas2.stage(iL,iP17).type = 'end';
        gas2 = count_Nstg(gas2);
        
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
        %%{
        print_states(gas1,iL,1:gas1.Nstg(iL)+1,Load);
        print_states(gas2,iL,1:gas2.Nstg(iL)+1,Load);
        print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
        print_states(fluidM,iL,1:fluidM.Nstg(iL)+1,Load);
        print_states(air,iL,1:air.Nstg(iL)+1,Load);
        %keyboard
        %}
        
        % Uncomment these lines to prints HEXs
        %%{
        plot_hex(HX(3),1,40,'K',true) % Heat rejection
        plot_hex(HX(7),1,41,'K',true) % Heat rejection
        plot_hex(HX(4),1,42,'K',true) % High Coupler
        plot_hex(HX(5),1,43,'K',true) % Cryo Coupler
        plot_hex(HX(8),1,44,'K',true) % Regenerator
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

keyboard
% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
HT = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
% Cold tanks
MT = run_tanks(MT,iL,fluidM,iM_out,iM_in,Load,T0);
% Atmospheric tanks
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);
%}
keyboard