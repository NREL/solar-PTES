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
%%%%%%%%%%%%%%%%%%%%

%%% INPUTS %%%
%%%%%%%%%%%%%%
Ran_ptop = 100e5;
Ran_Tbot = T0+15;
Ran_pbot = CP1('QT_INPUTS',0.0,Ran_Tbot,'P',steam(1).handle);
PR_dis   = Ran_ptop/Ran_pbot;
Ran_pmid1 = Ran_ptop/(PR_dis)^(1/3);  % pressure at HPT outlet
Ran_pmid2 = Ran_pmid1/(PR_dis)^(1/3); % pressure at IPT outlet
% Initial guesses
Ttop = HT.B(iL).T;
x1   = 0.10; % mass fraction of first steam extraction
x2   = 0.08; % mass fraction of second steam extraction
%%%%%%%%%%%%%%

% Initial guess of discharge conditions (Point 1)
i = 1;
steam(1).state(iL,i).T = Ttop;
steam(1).state(iL,i).p = Ran_ptop;
steam(1).state(iL,i).mdot = Load.mdot(iL);
[steam(1)] = update(steam(1),[iL,i],1);

% % Initial guess of conditions at points 8 (condenser outlet) and 9 (first
% % pump outlet) (not implemented yet)
% i = 8;
% steam(1).state(iL,i).T = Ran_Tbot - 2;
% steam(1).state(iL,i).p = Ran_pbot*(1 - ploss);
% steam(1).state(iL,i).mdot = Load.mdot(iL)*(1-x1)*(1-x2);
% [steam(1)] = update(steam(1),[iL,i],1);

% Set matrix of temperature and pressure points to test convergence
A_0 = [[steam(1).state(iL,:).T];[steam(1).state(iL,:).p]];
while 1
    fprintf(1,'Hello Rankine discharge PTES\n')
    
    % Set stage indices
    i  = 1;  % keeps track of the gas stage number
    iH = 1;  % keeps track of the Hot fluid stream number
    iE = 1;  % keeps track of the Environment (heat rejection) stream number
    
    % EXPAND (1-->2)
    [steam(1),i] = compexp(steam(1),[iL,i],eta,Ran_pmid1,4);
    
%     % Predict how much steam should be separated and mixed with water from
%     % point 11 (second pump outlet) in order to obtain wet saturated
%     % conditions (not implemented yet)
    
    
    % SEPARATE (2-->3)
    [steam(1),steam(2),i] = split_stream(steam(1),[iL,i],steam(2),[iL,1],x1);
    
    % REHEAT (gas-liquid) (3-->4)
    fluidH(iH).state(iL,1).T = HT.B(iL).T; fluidH(iH).state(iL,1).p = HT.B(iL).p;
    [fluidH(iH)] = update(fluidH(iH),[iL,1],1);
    [steam(1),fluidH(iH),i,~] = hex_TQ_cond(steam(1),[iL,i],fluidH(iH),[iL,1],eff,1.0,ploss,'hex',0,0);
    iH = iH + 1;
    
    % EXPAND (4-->5)
    [steam(1),i] = compexp(steam(1),[iL,i],eta,Ran_pmid2,4);
    
    % SEPARATE (5-->6)
    [steam(1),steam(3),i] = split_stream(steam(1),[iL,i],steam(2),[iL,1],x2);
    
    % EXPAND (6-->7)
    p_aim      = steam(1).state(iL,i).p/(PR_dis)^(1/3);
    [steam(1),i] = compexp(steam(1),[iL,i],eta,p_aim,4);
    
    % REJECT HEAT (external HEX) (7-->8)
    T_aim = Ran_Tbot - 2;
    [steam(1),environ,i,iE] = hex_set(steam(1),[iL,i],environ,[iL,iE],T_aim,eff,ploss);
    
    % COMPRESS (8-->9)
    p_aim = steam(3).state(iL,1).p;
    [steam(1),i] = compexp(steam(1),[iL,i],eta,p_aim,5);
    
    % MIX (9-->10)
    [steam(1),steam(3),i,~] = mix_streams(steam(1),[iL,i],steam(3),[iL,1]);
    
    % COMPRESS (10-->11)
    p_aim = steam(2).state(iL,1).p;
    [steam(1),i] = compexp(steam(1),[iL,i],eta,p_aim,5);
    
    % MIX (11-->12)
    [steam(1),steam(2),i,~] = mix_streams(steam(1),[iL,i],steam(2),[iL,1]);
    
    % COMPRESS (12-->13)
    p_aim = Ran_ptop;
    [steam(1),i] = compexp(steam(1),[iL,i],eta,p_aim,5);
    
    % HEAT (2-phase-liquid) (13-->1)
    fluidH(iH).state(iL,1).T = HT.B(iL).T; fluidH(iH).state(iL,1).p = HT.B(iL).p; %#ok<*SAGROW>
    [fluidH(iH)] = update(fluidH(iH),[iL,1],1);
    [steam(1),fluidH(iH),i,~] = hex_TQ_2p(steam(1),[iL,i],fluidH(iH),[iL,1],eff,ploss,'hex',2,1.20);
    iH = iH + 1;
    
    % Close cycle
    steam(1).stage(iL,i).type = steam(1).stage(iL,1).type;
    steam(2).stage(iL,2).type = 'end';
    steam(3).stage(iL,2).type = 'end';
    steam(1) = count_Nstg(steam(1));
    steam(2) = count_Nstg(steam(2));
    steam(3) = count_Nstg(steam(3));
        
    % Determine convergence and proceed
    A = [[steam(1).state(iL,:).T];[steam(1).state(iL,:).p]];

    %disp((A(A~=0) - A_0(A~=0))./A(A~=0)*100);
    if all(abs((A(A~=0) - A_0(A~=0))./A(A~=0))*100 < 1e-3) % is discharge cycle converged?
        % Exit discharge cycle
        gas_min_rho_dis = min([steam(1).state(iL,1:steam(1).Nstg(iL)).rho]); %take data for power density calculation
        break
    else
        % Set new initial conditions
        steam(1).state(iL,1) = steam(1).state(iL,i);
        A_0 = A;
    end
end

print_states(steam(1),iL,1:i,Load)
print_states(fluidH(1),iL,1:2,Load)
print_states(fluidH(2),iL,1:2,Load)

% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of (at least) one of the source tanks to stablish
% discharge time.

% Find t_dis
[MdotH] = total_Mdot(fluidH,[iL,1]);
t_dis  = HT.B(iL).M/MdotH;
Load.time(iL) = min([Load.time(iL),t_dis])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
[HT] = run_tanks(HT,fluidH,iL,Load,T0);
% Cold tanks. If they are not used (heat rejected to the environment), set
% the cold tanks to the same state as they were before this discharge
% process.
CT.A(iL+1) = CT.A(iL);
CT.B(iL+1) = CT.B(iL);