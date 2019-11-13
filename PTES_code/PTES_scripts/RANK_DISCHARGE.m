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

%%% INPUTS %%%
%%%%%%%%%%%%%%
Ran_ptop = 100e5;
Ran_Tbot = T0+15;
Ran_pbot = CP1('QT_INPUTS',0.0,Ran_Tbot,'P',steam.handle);
PR_dis   = Ran_ptop/Ran_pbot;
Ran_pmid1 = Ran_ptop/(PR_dis)^(1/3);  % pressure at HPT outlet
Ran_pmid2 = Ran_pmid1/(PR_dis)^(1/3); % pressure at IPT outlet
% Initial guesses
Ttop = HT.B(iL).T;
x1   = 0.10; % mass fraction of first steam extraction
x2   = 0.08; % mass fraction of second steam extraction
%%%%%%%%%%%%%%

% Set indices for secondary streams A and B
iSA = 15;
iSB = 17;

% Initial guess of discharge conditions (Point 1)
i = 1;
steam.state(iL,i).T = Ttop;
steam.state(iL,i).p = Ran_ptop;
steam.state(iL,i).mdot = Load.mdot(iL);
[steam] = update(steam,[iL,i],1);

% % Initial guess of conditions at points 8 (condenser outlet) and 9 (first
% % pump outlet) (not implemented yet)
% i = 8;
% steam.state(iL,i).T = Ran_Tbot - 2;
% steam.state(iL,i).p = Ran_pbot*(1 - ploss);
% steam.state(iL,i).mdot = Load.mdot(iL)*(1-x1)*(1-x2);
% [steam] = update(steam,[iL,i],1);

% Set matrix of temperature and pressure points to test convergence
A_0 = [[steam.state(iL,:).T];[steam.state(iL,:).p]];
while 1
    fprintf(1,'Hello Rankine discharge PTES\n')
    
    % Set stage indices
    i  = 1;  % keeps track of the gas stage number
    iH = 1;  % keeps track of the Hot fluid stream number
    iE = 1;  % keeps track of the Environment (heat rejection) stream number
    
    % EXPAND (1-->2)
    [steam,i] = compexp(steam,[iL,i],eta,Ran_pmid1,4);
    
%     % Predict how much steam should be separated and mixed with water from
%     % point 11 (second pump outlet) in order to obtain wet saturated
%     % conditions (not implemented yet)
    
    % SEPARATE (2-->3)
    [steam,i] = split_stream(steam,iL,i,iSA,x1);
    
    % REHEAT (gas-liquid) (3-->4)
    fluidH(iH).state(iL,1).T = HT.B(iL).T; fluidH(iH).state(iL,1).p = HT.B(iL).p;
    [fluidH(iH)] = update(fluidH(iH),[iL,1],1);
    [steam,fluidH(iH),i,~] = hex_TQ(steam,[iL,i],fluidH(iH),[iL,1],eff,ploss,'hex',2,1.0);
    iH = iH + 1;
    
    % EXPAND (4-->5)
    [steam,i] = compexp(steam,[iL,i],eta,Ran_pmid2,4);
    
    % SEPARATE (5-->6)
    [steam,i] = split_stream(steam,iL,i,iSB,x2);
    
    % EXPAND (6-->7)
    p_aim     = steam.state(iL,i).p/(PR_dis)^(1/3);
    [steam,i] = compexp(steam,[iL,i],eta,p_aim,4);
    
    % REJECT HEAT (external HEX) (7-->8)
    T_aim = Ran_Tbot - 2;
    [steam,environ,i,iE] = hex_set(steam,[iL,i],environ,[iL,iE],T_aim,eff,ploss);
    
    % COMPRESS (8-->9)
    p_aim = steam.state(iL,iSB).p;
    [steam,i] = compexp(steam,[iL,i],eta,p_aim,5);
    
    % MIX (9-->10)
    [steam,i,~] = mix_streams(steam,[iL,i],[iL,iSB]);
    
    % COMPRESS (10-->11)
    p_aim = steam.state(iL,iSA).p;
    [steam,i] = compexp(steam,[iL,i],eta,p_aim,5);
    
    % MIX (11-->12)
    [steam,i,~] = mix_streams(steam,[iL,i],[iL,iSA]);
    
    % COMPRESS (12-->13)
    p_aim = Ran_ptop;
    [steam,i] = compexp(steam,[iL,i],eta,p_aim,5);
    
    % HEAT (2-phase-liquid) (13-->1)
    fluidH(iH).state(iL,1).T = HT.B(iL).T; fluidH(iH).state(iL,1).p = HT.B(iL).p; %#ok<*SAGROW>
    [fluidH(iH)] = update(fluidH(iH),[iL,1],1);
    Taim = HT.A(iL).T;
    [steam,fluidH(iH),i,~] = hex_TQ(steam,[iL,i],fluidH(iH),[iL,1],eff,ploss,'hex',4,Taim);
    iH = iH + 1;
    
%     % Define heat exchanger geometry    
%     HX.shape     = 'circular';
%     HX.sigma     = 1e8;        % Maximum allowable stress, Pa
%     HX.L         = 3.0;        % Tube length, m
%     HX.D1        = 0.5e-2;     % Tube diameter, m
%     HX.t1        = 0.1*HX.D1;  % Tube thickness, m
%     HX.AfT       = 1.0;        % Total flow area, m2
%     HX.AfR       = 1.00;       % Ratio of flow areas, Af2/Af1, -
%     
%     % Code settings
%     HX.NX  = 100;               % Number of sections (grid)
%     HX.NI  = 1000;              % Maximum number of iterations
%     HX.TOL = 1e-2;              % Convergence tolerance, in %
%     
%     [steam,fluidH(iH),i,~] = hex_TQA(steam,[iL,i],fluidH(iH),[iL,1],HX,'hex',2,1.20);
%     iH = iH + 1;
    
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

print_states(steam,iL,1:(steam.Nstg(iL)+1),Load)
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