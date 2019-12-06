%%% PLANT LAYOUT %%%
%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%

% Set main pressure and temperature levels along the cycle (apart from
% Ran_ptop and Ran_Tbot, defined in INPUTS file)
if Load.options.useCold(iL)
    Ran_Tbot = Ran_TbotC;    
else
    Ran_Tbot = Ran_Tbot0;
end
Ran_pbot = CP1('QT_INPUTS',0.0,Ran_Tbot,'P',steam.handle);
PR_dis   = Ran_ptop/Ran_pbot;
Ran_pmid1 = Ran_ptop/(PR_dis)^(1/3);  % pressure at HPT outlet
Ran_pmid2 = Ran_pmid1/(PR_dis)^(1/3); % pressure at IPT outlet
Ran_Tmid1 = CP1('PQ_INPUTS',Ran_pmid1,0.0,'T',steam.handle); %wet saturated temp1
Ran_Tmid2 = CP1('PQ_INPUTS',Ran_pmid2,0.0,'T',steam.handle); %wet saturated temp2

% Set indices for first and second pump outlets
iP2 = 11; % same pressure as HPT outlet
iP1 = 9;  % same pressure as IPT outlet

% Set indices for secondary streams A and B
iSA = 15;
iSB = 17;

% Initial guess of discharge conditions (Point 1)
i = 1;
steam.state(iL,i).T = HT.B(iL).T;
steam.state(iL,i).p = Ran_ptop;
steam.state(iL,i).mdot = Load.mdot(iL);
[steam] = update(steam,[iL,i],1);

% Initial guess of conditions at point 11 (second pump outlet). This is
% necessary to predict the x1 fraction to obtain wet saturated steam at
% point 12.
i = 11;
steam.state(iL,i).T = Ran_Tmid2;
steam.state(iL,i).p = Ran_pmid1;
[steam] = update(steam,[iL,i],1);

% Initial guess of conditions at point 9 (first pump outlet). This is
% necessary to predict the x2 fraction to obtain wet saturated steam at
% point 10.
i = 9;
steam.state(iL,i).T = Ran_Tbot;
steam.state(iL,i).p = Ran_pmid2;
[steam] = update(steam,[iL,i],1);


% Set matrix of temperature and pressure points to test convergence
A_0 = [[steam.state(iL,:).T];[steam.state(iL,:).p]];
while 1
    fprintf(1,'Hello Rankine discharge PTES\n')
    
    % Set stage indices
    i  = 1;  % keeps track of the gas stage number
    iH = 1;  % keeps track of the Hot fluid stream number
    iC = 1;  % keeps track of the Cold fluid stream number
    iE = 1;  % keeps track of the Environment (heat rejection) stream number
    
    % EXPAND (1-->2)
    %[steam,i] = compexp(steam,[iL,i],eta,Ran_pmid1,4);
    [DEXP(1),steam,i] = compexp_func (DEXP(1),steam,[iL,i],'Paim',Ran_pmid1) ;
    
    % FIND x1 %
    %%%%%%%%%%%
    % x1 is the fraction of steam for the first steam extraction). It is
    % set to ensure that pump inlet temperature is just below wet saturated
    % conditions
    f1 = @(x1) dT_from_mix_obj(x1,steam,iL,i,iSA,iP2,Ran_Tmid1-1);
    %plot_function(f1,0.0,1.0,100,10)
    %TolX = 1e-6; %tolerance
    %options = optimset('TolX',TolX,'Display','iter');
    x1  = fzero(f1,[0.0,0.5]);%,options);
    %keyboard
    %%%%%%%%%%%
    
    % SEPARATE (2-->3)
    [steam,i] = split_stream(steam,iL,i,iSA,x1);
    
    % REHEAT (gas-liquid) (3-->4)
    fluidH(iH).state(iL,1).T = HT.B(iL).T; fluidH(iH).state(iL,1).p = HT.B(iL).p;
    [fluidH(iH)] = update(fluidH(iH),[iL,1],1);
    [steam,fluidH(iH),i,~,HX_REHEAT] = hex_TQ(steam,[iL,i],fluidH(iH),[iL,1],eff,ploss,'hex',2,1.0);
    iH = iH + 1;
    
    % EXPAND (4-->5)
    %[steam,i] = compexp(steam,[iL,i],eta,Ran_pmid2,4);
    [DEXP(2),steam,i] = compexp_func (DEXP(2),steam,[iL,i],'Paim',Ran_pmid2) ;
    
    % FIND x2 %
    %%%%%%%%%%%
    % x2 is the fraction of steam for the second steam extraction). It is
    % set to ensure that pump inlet temperature is just below wet saturated
    % conditions
    f2 = @(x2) dT_from_mix_obj(x2,steam,iL,i,iSB,iP1,Ran_Tmid2-1);
    %plot_function(f1,0.0,1.0,100,10)
    %TolX = 1e-6; %tolerance
    %options = optimset('TolX',TolX,'Display','iter');
    x2  = fzero(f2,[0.0,0.5]);%,options);
    %keyboard
    %%%%%%%%%%%
    
    % SEPARATE (5-->6)
    [steam,i] = split_stream(steam,iL,i,iSB,x2);
    
    % EXPAND (6-->7)
    p_aim     = steam.state(iL,i).p/(PR_dis)^(1/3);
    %[steam,i] = compexp(steam,[iL,i],eta,p_aim,4);
    [DEXP(3),steam,i] = compexp_func (DEXP(3),steam,[iL,i],'Paim',p_aim) ;
    
    if Load.options.useCold(iL)
        % COOL (condense using cold tanks)
        fluidC(iC).state(iL,1).T = CT.B(iL).T; fluidC(iC).state(iL,1).p = CT.B(iL).p; %#ok<*SAGROW>
        [fluidC(iC)] = update(fluidC(iC),[iL,1],1);
        T_aim = CP1('PQ_INPUTS',steam.state(iL,i).p,0.0,'T',steam.handle) - 1; %wet saturated
        HX_CONDEN.model = 'eff';
        HX_CONDEN.eff = eff;
        HX_CONDEN.ploss = ploss;
        HX_CONDEN.stage_type = 'hex';
        HX_CONDEN.NX = 100;
        [HX_CONDEN, steam, fluidC(iC), i, ~] = hex_func(HX_CONDEN,iL,steam,i,fluidC(iC),1,6,T_aim);
        iC=iC+1;
    else
        % REJECT HEAT (external HEX) (7-->8)
        T_aim = Ran_Tbot - 1;
        [steam,environ,i,iE] = hex_set(steam,[iL,i],environ,[iL,iE],T_aim,eff,ploss);
    end
    
    % COMPRESS (8-->9)
    p_aim = steam.state(iL,iSB).p;
    %[steam,i] = compexp(steam,[iL,i],eta,p_aim,5);
    [DCMP(1),steam,i] = compexp_func (DCMP(1),steam,[iL,i],'Paim',p_aim) ;
    
    % MIX (9-->10)
    [steam,i,~] = mix_streams(steam,[iL,i],[iL,iSB]);
    
    % COMPRESS (10-->11)
    p_aim = steam.state(iL,iSA).p;
    %[steam,i] = compexp(steam,[iL,i],eta,p_aim,5);
    [DCMP(2),steam,i] = compexp_func (DCMP(2),steam,[iL,i],'Paim',p_aim) ;
    
    % MIX (11-->12)
    [steam,i,~] = mix_streams(steam,[iL,i],[iL,iSA]);
    
    % COMPRESS (12-->13)
    p_aim = Ran_ptop;
    %[steam,i] = compexp(steam,[iL,i],eta,p_aim,5);
    [DCMP(3),steam,i] = compexp_func (DCMP(3),steam,[iL,i],'Paim',p_aim) ;
    
    % HEAT (2-phase-liquid) (13-->1)
    fluidH(iH).state(iL,1).T = HT.B(iL).T; fluidH(iH).state(iL,1).p = HT.B(iL).p; %#ok<*SAGROW>
    [fluidH(iH)] = update(fluidH(iH),[iL,1],1);
    Taim = HT.A(iL).T;
    [steam,fluidH(iH),i,~,HX_BOILER] = hex_TQ(steam,[iL,i],fluidH(iH),[iL,1],eff,ploss,'hex',4,Taim);
    iH = iH + 1;
    
    % Close cycle
    steam.stage(iL,i).type = 'end';
    steam.stage(iL,iSA+1).type = 'end';
    steam.stage(iL,iSB+1).type = 'end';
    steam = count_Nstg(steam);
    
    % Determine convergence and proceed
    A = [[steam.state(iL,:).T];[steam.state(iL,:).p]];

    %disp((A(A~=0) - A_0(A~=0))./A(A~=0)*100);
    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is discharge cycle converged?
        % Exit discharge cycle
        gas_min_rho_dis = min([steam.state(iL,1:steam.Nstg(iL)).rho]); %take data for power density calculation
        break
    else
        % Set new initial conditions
        steam.state(iL,1) = steam.state(iL,i);
        A_0 = A;
    end
end

% print_states(steam,iL,1:(steam.Nstg(iL)+1),Load)
% print_states(fluidH(1),iL,1:2,Load)
% print_states(fluidH(2),iL,1:2,Load)

% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of (at least) one of the source tanks to stablish
% discharge time.

% Find t_dis
[MdotH] = total_Mdot(fluidH,[iL,1]);
t_disH  = HT.B(iL).M/MdotH;
[MdotC] = total_Mdot(fluidC,[iL,1]);
t_disC  = CT.B(iL).M/MdotC;
Load.time(iL) = min([Load.time(iL),t_disH,t_disC])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
[HT] = run_tanks(HT,fluidH,iL,Load,T0);
% Cold tanks
[CT] = run_tanks(CT,fluidC,iL,Load,T0);
% % Cold tanks. If they are not used (heat rejected to the environment), set
% % the cold tanks to the same state as they were before this discharge
% % process.
% CT.A(iL+1) = CT.A(iL);
% CT.B(iL+1) = CT.B(iL);


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