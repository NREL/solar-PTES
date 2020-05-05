%%% RANKINE PLANT LAYOUT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1: Boiler outlet
%  2: HPT outlet
%  3: Separator outlet (first steam extraction)
%  4: Reheater outlet
%  5: IPT outlet
%  6: Separator outlet (second stream extraction)
%  7: LPT outlet
%  8: Condenser outlet
%  9: Pump outlet (mixing with second steam extraction)
% 10: Mixer outlet
% 11: Pump outlet (mixing with first steam extraction)
% 12: Mixer outlet
% 13: Pump outlet
% 14: Boiler outlet (end main stream, same as point 1 at convergence)
% 15: Separator outlet (first steam extraction)  - secondary stream A
% 16: Mixer outlet (end secondary stream A)
% 17: Separator outlet (second steam extraction) - secondary stream B
% 18: Mixer outlet (end secondary stream B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set indices for first and second pump outlets
iP2 = 11; % same pressure as HPT outlet
iP1 = 9;  % same pressure as IPT outlet

% Set indices for secondary streams A and B
iSA = 15;
iSB = 17;

% Set main pressure and temperature levels along the cycle (apart from
% Ran_ptop and Ran_Tbot, defined in INPUTS file)
if Load.options.useCold(iL)
    Ran_Tbot = Ran_TbotC;    
else
    Ran_Tbot = Ran_Tbot0;
end
Ran_pbot  = RP1('QT_INPUTS',0.0,Ran_Tbot,'P',steam);

if Ran_pbot < Ran_pbotMIN % Check bottom pressure doesn't exceed lower bound
    Ran_pbot = Ran_pbotMIN ;
    Ran_Tbot  = RP1('PQ_INPUTS',Ran_pbot,0.0,'T',steam);
end

if design_mode == 1

    PR_dis    = Ran_ptop/Ran_pbot ; % Total pressure ratio
    PR_dis0   = PR_dis; % Total pressure ratio - deisgn point
    Ran_pmid1 = Ran_ptop/(PR_dis0)^(1/3);  % pressure at HPT outlet. First two stages pressure ratios are kept constant. Only final LP stage pressure ratio changes.
    Ran_pmid2 = Ran_pmid1/(PR_dis0)^(1/3); % pressure at IPT outlet
    Ran_Tmid1 = RP1('PQ_INPUTS',Ran_pmid1,0.0,'T',steam); %wet saturated temp1
    Ran_Tmid2 = RP1('PQ_INPUTS',Ran_pmid2,0.0,'T',steam); %wet saturated temp2
   
    % Initial guess of discharge conditions (Point 1)
    iG = 1;
    steam.state(iL,iG).T = HT.B(iL).T;
    steam.state(iL,iG).p = Ran_ptop;
    steam.state(iL,iG).mdot = Load.mdot(iL);
    [steam] = update(steam,[iL,iG],1);
    
    % Initial guess of conditions at point 11 (second pump outlet). This is
    % necessary to predict the x1 fraction to obtain wet saturated steam at
    % point 12.
    iG = 11;
    steam.state(iL,iG).T = Ran_Tmid2;
    steam.state(iL,iG).p = Ran_pmid1;
    [steam] = update(steam,[iL,iG],1);
    
    % Initial guess of conditions at point 9 (first pump outlet). This is
    % necessary to predict the x2 fraction to obtain wet saturated steam at
    % point 10.
    iG = 9;
    steam.state(iL,iG).T = Ran_Tbot;
    steam.state(iL,iG).p = Ran_pmid2;
    [steam] = update(steam,[iL,iG],1);
    
    
    % Set matrix of temperature and pressure points to test convergence
    A_0 = [[steam.state(iL,:).T];[steam.state(iL,:).p]];
    
else
    
    for ii = 1 : numel(A(A~=0))/2
        steam.state(iL,ii).T    = A(1,ii) ;
        steam.state(iL,ii).p    = A(2,ii) ;
        steam.state(iL,ii).mdot = Load.mdot(iL) ;
        
        % For inventory control, assume that the pressure scales with the off-design mass flow rate
        steam.state(iL,ii).p = steam.state(iL,ii).p * Load.mdot(iL) / DEXP(1).mdot0 ;
        [steam] = update(steam,[iL,ii],1);
    end 
    
    Ran_ptop  = Ran_ptop * Load.mdot(iL) / DEXP(1).mdot0 ; % This accounts for part-load operation
    PR_dis    = Ran_ptop/Ran_pbot; % Total pressure ratio
    Ran_pmid1 = Ran_ptop/(PR_dis0)^(1/3);  % pressure at HPT outlet. First two stages pressure ratios are kept constant. Only final LP stage pressure ratio changes.
    Ran_pmid2 = Ran_pmid1/(PR_dis0)^(1/3); % pressure at IPT outlet
    Ran_Tmid1 = RP1('PQ_INPUTS',Ran_pmid1,0.0,'T',steam); %wet saturated temp1
    Ran_Tmid2 = RP1('PQ_INPUTS',Ran_pmid2,0.0,'T',steam); %wet saturated temp2
    
    DEXP(3).mdot(iL) = DEXP(3).mdot0 ;
    DEXP(3).pr(iL) = 1.0 ;
    
    DEXP(3).eta0 = 0.9 * sqrt(Load.mdot(iL) / DEXP(1).mdot0) ;
    
end
max_iter=100;
for counter=1:max_iter
    fprintf(1,['Discharging RANKINE. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    
    % Set stage indices
    iG = 1;  % keeps track of the gas stage number
    iH = 1;  % keeps track of the Hot fluid stream number
    iC = 1;  % keeps track of the Cold fluid stream number
    iE = 1;  % keeps track of the Environment (heat rejection) stream number
    iA = 1;  % keeps track of the Air (heat rejection) stream number
    iPMP = 1 ; % Keeps track of which pump is being used
    
    % EXPAND (1-->2)
    [DEXP(1),steam,iG] = compexp_func (DEXP(1),iL,steam,iG,'Paim',Ran_pmid1,design_mode) ;
    
    % FIND x1 %
    %%%%%%%%%%%
    % x1 is the fraction of steam for the first steam extraction). It is
    % set to ensure that pump inlet temperature is just below wet saturated
    % conditions
    if ~design_mode
        steam.state(iL,iP2).p = steam.state(iL,iG).p ;
        [steam] = update(steam,[iL,iP2],1);
        steam.state(iL,iSA).p = steam.state(iL,iG).p ;
        [steam] = update(steam,[iL,iSA],1);
        
        x1 = DEXP(1).mdot0 - DEXP(2).mdot0 ;
    else
        f1 = @(x1) dT_from_mix_obj(x1,steam,iL,iG,iSA,iP2,Ran_Tmid1-1);
        %plot_function(f1,0.0,1.0,100,10)
        %TolX = 1e-6; %tolerance
        %options = optimset('TolX',TolX,'Display','iter');
        x1  = fzero(f1,[0.0,0.5]);%,options);
        %keyboard
        %%%%%%%%%%%    
    end
    
    % SEPARATE (2-->3)
    [steam,iG] = split_stream(steam,iL,iG,iSA,x1);
    
    % REHEAT (gas-liquid) (3-->4)
    fluidH.state(iL,iH).T = HT.B(iL).T; fluidH.state(iL,iH).p = HT.B(iL).p;
    [fluidH] = update(fluidH,[iL,iH],1);
    Taim = HT.A(iL).T;
    %[HX(ihx_JB+1),fluidH,iH,steam,iG] = set_hex(HX(ihx_JB+1),iL,fluidH,iH,steam,iG,4,Taim);
    [HX(ihx_JB+1),fluidH,iH,steam,iG] = hex_func(HX(ihx_JB+1),iL,fluidH,iH,steam,iG,4,Taim);
    [DPMP(iPMP),fluidH,iH] = compexp_func (DPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
    iH=iH+1; iPMP=iPMP+1;
    
    % EXPAND (4-->5)
    [DEXP(2),steam,iG] = compexp_func (DEXP(2),iL,steam,iG,'Paim',Ran_pmid2,design_mode) ;
    
    % FIND x2 %
    %%%%%%%%%%%
    % x2 is the fraction of steam for the second steam extraction). It is
    % set to ensure that pump inlet temperature is just below wet saturated
    % conditions
    if ~design_mode
        steam.state(iL,iP1).p = steam.state(iL,iG).p ;
        [steam] = update(steam,[iL,iP1],1);
        steam.state(iL,iSB).p = steam.state(iL,iG).p ;
        [steam] = update(steam,[iL,iSB],1);
        
        % The mass flow rate through the final stage is determined by
        % Stodola's ellipse if pcond < pcond0. If pcond >= pcond0, then
        % the mass flow is calculated in the original way, but the
        % efficiency of DEXP(3) decreases a bit as flow separates on the
        % blades.
        if DEXP(3).pr(iL) > 1.0
            x2 = DEXP(2).mdot(iL) - DEXP(3).mdot(iL) ;
        else
            f2 = @(x2) dT_from_mix_obj(x2,steam,iL,iG,iSB,iP1,Ran_Tmid2-1);
            x2  = fzero(f2,[0.0,0.5]);
        end
    else
        f2 = @(x2) dT_from_mix_obj(x2,steam,iL,iG,iSB,iP1,Ran_Tmid2-1);
        %plot_function(f1,0.0,1.0,100,10)
        %TolX = 1e-6; %tolerance
        %options = optimset('TolX',TolX,'Display','iter');
        x2  = fzero(f2,[0.0,0.5]);%,options);
        %keyboard
        %%%%%%%%%%%
       
    end
    
    % SEPARATE (5-->6)
    [steam,iG] = split_stream(steam,iL,iG,iSB,x2);
    % EXPAND (6-->7)
    p_aim     = Ran_pbot ; %steam.state(iL,iG).p/(PR_dis)^(1/3);
    [DEXP(3),steam,iG] = compexp_func (DEXP(3),iL,steam,iG,'Paim',p_aim,design_mode) ;
    
    if Load.options.useCold(iL)
        % COOL (condense using cold tanks)
        fluidC.state(iL,iC).T = CT.B(iL).T; fluidC.state(iL,iC).p = CT.B(iL).p; %#ok<*SAGROW>
        [fluidC] = update(fluidC,[iL,iC],1);
        T_aim = RP1('PQ_INPUTS',steam.state(iL,iG).p,0.0,'T',steam) - 1; %wet saturated
        %[HX(ihx_JB+2),steam,iG,fluidC,iC] = set_hex(HX(ihx_JB+2),iL,steam,iG,fluidC,iC,5,T_aim);
        [HX(ihx_JB+2),steam,iG,fluidC,iC] = hex_func(HX(ihx_JB+2),iL,steam,iG,fluidC,iC,5,T_aim);
        [DPMP(iPMP),fluidC,iC] = compexp_func (DPMP(iPMP),iL,fluidC,iC,'Paim',fluidC.state(iL,1).p,1);
        iC=iC+1; iPMP=iPMP+1;
    else
        % REJECT HEAT (external HEX) (7-->8)
        T_aim = Ran_Tbot - 1;
        %[steam,environ,iG,iE] = hex_set(steam,[iL,iG],environ,[iL,iE],T_aim,eff,ploss);
        air.state(iL,1).T = T0; air.state(iL,1).p = p0; air = update(air,[iL,1],1);
        %[HX(ihx_JB+3), steam, iG, air, iA] = set_hex(HX(ihx_JB+3),iL,steam,iG,air,iA,5,T_aim);
        [HX(ihx_JB+3), steam, iG, air, iA] = hex_func(HX(ihx_JB+3),iL,steam,iG,air,iA,5,T_aim);
        [DFAN(1),air,iA] = compexp_func (DFAN(1),iL,air,iA,'Paim',p0,1);
    end
    
    % COMPRESS (8-->9)
    p_aim = steam.state(iL,iSB).p;
    [DCMP(1),steam,iG] = compexp_func (DCMP(1),iL,steam,iG,'Paim',p_aim,1);
    
    % MIX (9-->10)
    [steam,iG,~] = mix_streams(steam,[iL,iG],[iL,iSB]);
    
    % COMPRESS (10-->11)
    p_aim = steam.state(iL,iSA).p;
    [DCMP(2),steam,iG] = compexp_func (DCMP(2),iL,steam,iG,'Paim',p_aim,1);
    
    % MIX (11-->12)
    [steam,iG,~] = mix_streams(steam,[iL,iG],[iL,iSA]);
    
    % COMPRESS (12-->13)
    p_aim = Ran_ptop;
    [DCMP(3),steam,iG] = compexp_func (DCMP(3),iL,steam,iG,'Paim',p_aim,1   );
    
    % HEAT (2-phase-liquid) (13-->1)
    fluidH.state(iL,iH).T = HT.B(iL).T; fluidH.state(iL,iH).p = HT.B(iL).p; %#ok<*SAGROW>
    [fluidH] = update(fluidH,[iL,iH],1);
    Taim = HT.A(iL).T;
    %[HX(ihx_JB+4),fluidH,iH,steam,iG] = set_hex(HX(ihx_JB+4),iL,fluidH,iH,steam,iG,4,Taim);
    [HX(ihx_JB+4),fluidH,iH,steam,iG] = hex_func(HX(ihx_JB+4),iL,fluidH,iH,steam,iG,4,Taim);
    [DPMP(iPMP),fluidH,iH] = compexp_func (DPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
    iH=iH+1; iPMP=iPMP+1;
    
    % Determine convergence and proceed
    A = [[steam.state(iL,:).T];[steam.state(iL,:).p]];

    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) || counter==max_iter % is discharge cycle converged?
        % Close working fluid streams
        steam.stage(iL,iG).type = 'end';
        steam.stage(iL,iSA+1).type = 'end';
        steam.stage(iL,iSB+1).type = 'end';
        steam = count_Nstg(steam);
        
        % Close air (heat rejection) streams
        iA_out = 1:5:(iA-1); iA_in  = iA_out + 2;
        for i=iA_in, air.stage(iL,i).type = 'end'; end
        air = count_Nstg(air);
        
        % Close storage fluid streams
        iH_out = 1:3:(iH-1); iH_in  = iH_out + 2;
        iC_out = 1:3:(iC-1); iC_in  = iC_out + 2;
        for i=iH_in, fluidH.stage(iL,i).type = 'end'; end
        for i=iC_in, fluidC.stage(iL,i).type = 'end'; end
        fluidH = count_Nstg(fluidH);
        fluidC = count_Nstg(fluidC);
        
        % Uncomment these lines to print states
        
        %print_states(steam,iL,1:steam.Nstg(iL)+1,Load);
        %print_states(air,iL,1:air.Nstg(iL)+1,Load);
        %print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
        %print_states(fluidC,iL,1:fluidC.Nstg(iL)+1,Load);
        %keyboard
        
        
        % Exit loop
        break
    else
        if ~design_mode
            
            %print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
%              fprintf('p1:   %13.8f \n',gas.state(iL,1).p/1e5)
%              fprintf('pEND: %13.8f \n\n',gas.state(iL,iG).p/1e5)
%              fprintf('T1:   %13.8f \n',gas.state(iL,1).T)
%              fprintf('TEND: %13.8f \n\n',gas.state(iL,iG).T)
            % Adjust inlet pressure to try to reach convergence. The
            % 'smoothing' factor has to be quite small (<0.1, say) for this to be stable
            gas.state(iL,1).p = gas.state(iL,1).p - 0.05 * (gas.state(iL,iG).p - gas.state(iL,1).p) ;
            gas.state(iL,1).T = gas.state(iL,1).T + 0.05 * (gas.state(iL,iG).T - gas.state(iL,1).T) ;
            gas.state(iL,1).mdot = Load.mdot(iL);
            [gas] = update(gas,[iL,1],1);
        else
            % Set new initial conditions
            steam.state(iL,1) = steam.state(iL,iG);
        end
        A_0 = A;
    end
end
if counter==max_iter
    warning('Exiting JB_CHARGE cycle without having reached convergence');
end

% Compute total mass flow rates of the hot and cold storage fluids and
% compute the end of discharge time (stop when one tank becomes empty)
[MdotH] = total_mdot(fluidH, iL, iH_out);
t_disH  = HT.B(iL).M/MdotH;
[MdotC] = total_mdot(fluidC, iL, iC_out);
t_disC  = CT.B(iL).M/MdotC;
Load.time(iL) = min([Load.time(iL),t_disH,t_disC])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
[HT] = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
% Cold tanks
if Load.options.useCold(iL)
    [CT] = run_tanks(CT,iL,fluidC,iC_out,iC_in,Load,T0);
else
    % If they are not used (heat rejected to the environment), set the cold
    % tanks to the same state as they were before the discharge process.
    CT.A(iL+1) = CT.A(iL);
    CT.B(iL+1) = CT.B(iL);
end
% Atmospheric tanks
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);


%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = dT_from_mix_obj(x1,steam,iL,i,iS,iP,T_obj)
% Obtain the difference in temperature between objective and computed
% values, after steam separation and feed water heating.

% Update mass flow rate at point iP
steam.state(iL,iP).mdot = steam.state(iL,i).mdot*(1-x1);

% Separate stream at point i, obtain stream at point iS
[steam,~] = split_stream(steam,iL,i,iS,x1);

% Mix stream at point iP with stream at point iS
[steam,i,~] = mix_streams(steam,[iL,iP],[iL,iS]);

% Compare temperature with objective temperature
T = steam.state(iL,i).T;
err = T - T_obj;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%