fprintf(1,'Hello discharge LAES\n')

% Set stage indices
iG1 = 1;  % keeps track of the gas stage number (main stream)
iH  = 1;  % keeps track of the Hot fluid stream number
iC1 = 1;  % keeps track of the Cold fluid stream number
iC2 = 1;  % keeps track of the Cold fluid stream number
iE  = 1;  % keeps track of the heat rejection stream number

% Compute PR_dis based on charge pressure ratio and PR_fact
PR_dis = PR_fact*PR;

% Set the initial discharge conditions.
gas1.state(2,1).mdot = mdot*(1-Qual);
gas1.state(2,1).h    = LA(3).H/LA(3).M;
gas1.state(2,1).p    = LA(3).p;
[gas1] = update(gas1,[2,1],2);

% COMPRESS
p_aim = gas1.state(2,iG1).p*PR_dis*p0/pLA;
[gas1,iG1] = compexp(gas1,[2,iG1],eta,p_aim,5);

% HEAT (gas-fluid)
fluidC2.state(2,1).T = CT2.B(3).T; fluidC2.state(2,1).p = CT2.B(3).p; fluidC2.state(2,1).mdot = fluidC2.state(1,1).mdot;
[fluidC2] = update(fluidC2,[2,1],1);
[fluidC2,gas1,~,iG1] = hex_TQ_cond(fluidC2,[2,1],gas1,[2,iG1],eff,1.0,ploss,'hex',4,0);

% HEAT (gas-fluid)
fluidC1.state(2,1).T = CT1.B(3).T; fluidC1.state(2,1).p = CT1.B(3).p; fluidC1.state(2,1).mdot = fluidC1.state(1,1).mdot;
[fluidC1] = update(fluidC1,[2,1],1);
[fluidC1,gas1,~,iG1] = hex_TQ_cond(fluidC1,[2,1],gas1,[2,iG1],eff,1.0,ploss,'hex',4,0);

% REJECT HEAT (external HEX)
T_aim = environ.T0;
[gas1,environ,iG1,iE] = hex_set(gas1,[2,iG1],environ,[2,iE],T_aim,eff,ploss);

for iN = 1:Ne_dis
    % HEAT (gas-fluid)
    fluidH(iH).state(2,1).T = HT.B(3).T; fluidH(iH).state(2,1).p = HT.B(3).p; fluidH(iH).state(2,1).mdot = fluidH(iH).state(1,1).mdot; %#ok<*SAGROW>
    [fluidH(iH)] = update(fluidH(iH),[2,1],1);
    [fluidH(iH),gas1,~,iG1] = hex_TQ_cond(fluidH(iH),[2,1],gas1,[2,iG1],eff,1.0,ploss,'hex',4,0);
    iH=iH+1;
    
    % EXPAND
    PRe_dis = (gas1.state(2,iG1).p/p0)^(1/(Ne_dis+1-iN));  % expansion pressure ratio
    p_aim = gas1.state(2,iG1).p/PRe_dis;
    [gas1,iG1] = compexp(gas1,[2,iG1],eta,p_aim,1);
end


fprintf(1,'\n %10s %10s %10s %10s %10s %10s %10s %10s','T(K)','p(bar)','h(kJ/kg)','s(kJ/kg/K)','m(kg/s)','Q','Stage','Ind');
for i0=1:iG1, fprintf(1,'\n %10.1f %10.1f %10.1f %10.2f %10.1f %10.1f %10s %10d',gas1.state(2,i0).T,gas1.state(2,i0).p/1e5,gas1.state(2,i0).h/1e3,gas1.state(2,i0).s/1e3,gas1.state(2,i0).mdot,gas1.state(2,i0).Q,gas1.stage(2,i0).type,i0); end; fprintf(1,'\n');

% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of (at least) one of the source tanks to stablish
% discharge time.

% Find t_dis (minimum for both cycles to avoid depletion)
[MdotH]  = total_Mdot(fluidH,[2,1]);
t_dis    = HT.B(3).M/MdotH;
[MdotC1] = total_Mdot(fluidC1,[2,1]);
tC1_dis  = CT1.B(3).M/MdotC1;
[MdotC2] = total_Mdot(fluidC2,[2,1]);
tC2_dis  = CT2.B(3).M/MdotC2;
MdotLA   = gas1.state(2,1).mdot;
tLA_dis  = LA(3).M/MdotLA;
t_dis    = min([t_dis,tC1_dis,tC2_dis,tLA_dis]);

% During discharge, B is source and A is sink
% Hot tanks
[HT.B(4), HT.A(3), HT.A(4), WL_mixH_dis]  = liquid_tanks_compute(fluidH, 2,2,HT.B(3),t_dis,T0);
% Cold tanks
[CT1.B(4),CT1.A(3),CT1.A(4),WL_mixC1_dis] = liquid_tanks_compute(fluidC1,2,2,CT1.B(3),t_dis,T0);
[CT2.B(4),CT2.A(3),CT2.A(4),WL_mixC2_dis] = liquid_tanks_compute(fluidC2,2,2,CT2.B(3),t_dis,T0);
WL_mix_C_dis = WL_mixC1_dis + WL_mixC2_dis;
% Air source/sink
LA(4).T = LA(3).T;
LA(4).p = LA(3).p;
LA(4).M = LA(3).M - gas1.state(2,1).mdot*t_dis;
LA(4).V = 1/gas1.state(2,1).rho*LA(4).M;
LA(4).H = gas1.state(2,1).h*LA(4).M;
LA(4).S = gas1.state(2,1).s*LA(4).M;
LA(4).B = LA(4).H - T0*LA(4).S;
AA(4).T = gas1.state(2,iG1).T;
AA(4).p = gas1.state(2,iG1).p;
AA(4).M = gas1.state(2,iG1).mdot*t_dis;
AA(4).V = 1/gas1.state(2,iG1).rho*AA(4).M;
AA(4).H = gas1.state(2,iG1).h*AA(4).M;
AA(4).S = gas1.state(2,iG1).s*AA(4).M;
AA(4).B = AA(4).H - T0*AA(4).S;

