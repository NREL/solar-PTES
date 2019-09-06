% Set stage indices
i  = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the Environment (heat rejection) stream number

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
Ran_pbot = CP1('QT_INPUTS',0.0,Ran_Tbot,'P',steamA.handle);
PR_dis   = Ran_ptop/Ran_pbot;
% Initial guesses
Ttop = HT.B(iL).T;
x1   = 0.1; % mass fraction of first steam extraction
x2   = 0.1; % mass fraction of second steam extraction
%%%%%%%%%%%%%%

% Initial guess of discharge conditions (Point 1)
i = 1;
steamA.state(iL,i).T = Ttop;
steamA.state(iL,i).p = Ran_ptop;
steamA.state(iL,i).mdot = Load.mdot(iL);
[steamA] = update(steamA,[iL,i],1);

% EXPAND (1-->2)
p_aim      = steamA.state(iL,i).p/(PR_dis)^(1/3);
[steamA,i] = compexp(steamA,[iL,i],eta,p_aim,4);

% SEPARATE (2-->3)
[steamA,steamB,i] = split_stream(steamA,[iL,i],steamB,[iL,1],x1);

% REHEAT (gas-liquid) (3-->4)
fluidH(iH).state(iL,1).T = HT.B(iL).T; fluidH(iH).state(iL,1).p = HT.B(iL).p;
[fluidH(iH)] = update(fluidH(iH),[iL,1],1);
[steamA,fluidH(iH),i,~] = hex_TQ_cond(steamA,[iL,i],fluidH(iH),[iL,1],eff,1.0,ploss,'hex',0,0);

% EXPAND (4-->5)
p_aim      = steamA.state(iL,i).p/(PR_dis)^(1/3);
[steamA,i] = compexp(steamA,[iL,i],eta,p_aim,4);

% SEPARATE (5-->6)
[steamA,steamC,i] = split_stream(steamA,[iL,i],steamB,[iL,1],x2);

% EXPAND (6-->7)
p_aim      = steamA.state(iL,i).p/(PR_dis)^(1/3);
[steamA,i] = compexp(steamA,[iL,i],eta,p_aim,4);

% REJECT HEAT (external HEX) (7-->8)
T_aim = Ran_Tbot - 2;
[steamA,environ,i,iE] = hex_set(steamA,[iL,i],environ,[iL,iE],T_aim,eff,ploss);

% COMPRESS (8-->9)
p_aim = steamC.state(iL,1).p;
[steamA,i] = compexp(steamA,[iL,i],eta,p_aim,5);

% MIX (9-->10)
[steamA,steamC,i,~] = mix_streams(steamA,[iL,i],steamC,[iL,1]);

% COMPRESS (10-->11)
p_aim = steamB.state(iL,1).p;
[steamA,i] = compexp(steamA,[iL,i],eta,p_aim,5);

% MIX (11-->12)
[steamA,steamB,i,~] = mix_streams(steamA,[iL,i],steamB,[iL,1]);

% COMPRESS (12-->13)
p_aim = Ran_ptop;
[steamA,i] = compexp(steamA,[iL,i],eta,p_aim,5);

% HEAT (2-phase-liquid) (13-->1)
fluidH(iH).state(iL,1).T = HT.B(iL).T; fluidH(iH).state(iL,1).p = HT.B(iL).p;
[fluidH(iH)] = update(fluidH(iH),[iL,1],1);
[steamA,fluidH(iH),i,~] = hex_TQ_cond(steamA,[iL,i],fluidH(iH),[iL,1],eff,1.0,ploss,'hex',0,0);



print_states(steamA,iL,1:i,Load)
print_states(steamB,iL,1:2,Load)
print_states(steamC,iL,1:2,Load)

disp('hey')
%%% Not implemented yet %%%
keyboard
error('not implemented yet')
