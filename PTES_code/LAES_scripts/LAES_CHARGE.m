fprintf(1,'Hello charge LAES\n')

% Set stage indices
iG1 = 1;  % keeps track of the gas stage number (main stream)
iG2 = 1;  % keeps track of the gas stage number (secondary)
iG3 = 1;  % keeps track of the gas stage number (secondary)
iH  = 1;  % keeps track of the Hot fluid stream number
iC1 = 1;  % keeps track of the Cold fluid stream number
iC2 = 1;  % keeps track of the Cold fluid stream number
iE  = 1;  % keeps track of the heat rejection stream number

% Define initial conditions
% Compressor inlet
gas1.state(1,1).p    = p0; gas1.state(1,1).T = T0;
gas1.state(1,1).mdot = mdot;
[gas1] = update(gas1,[1,1],1);
% Returning vapour stream
gas3.state(1,1).p = pLA;
gas3.state(1,1).h = CP1('PQ_INPUTS',pLA,1.0,'H',gas3.handle);
[gas3] = update(gas3,[1,1],2);

for iN = 1:Nc_ch
    % COMPRESS
    PRc = (pmax/gas1.state(1,iG1).p)^(1/(Nc_ch+1-iN));
    p_aim = gas1.state(1,iG1).p*PRc;
    [gas1,iG1] = compexp(gas1,[1,iG1],eta,p_aim,3);
    ptop  = gas1.state(1,iG1).p;
    
    % COOL (gas-liquid)
    fluidH(iH).state(1,1).T = HT.A(1).T; fluidH(iH).state(1,1).p = HT.A(1).p; %#ok<*SAGROW>
    [fluidH(iH)] = update(fluidH(iH),[1,1],1);
    [gas1,fluidH(iH),iG1,~] = hex_TQ_cond(gas1,[1,iG1],fluidH(iH),[1,1],eff,1.0,ploss,'hex',0,0);
    iH=iH+1;
    
    % REJECT HEAT (external HEX)
    T_aim = environ.T0;
    [gas1,environ,iG1,iE] = hex_set(gas1,[1,iG1],environ,[1,iE],T_aim,eff,ploss);
end

% Repeat until convergence
Qual   = xAir; %initial guess
Qual_0 = Qual;
iG1_0  = iG1; iG2_0 = iG2; iG3_0 = iG3;
iH_0   = iH;  iC1_0 = iC1; iC2_0 = iC2; iE_0 = iE;
while 1
    % SPLIT main and second streams
    [gas1,gas2,iG1] = split_stream(gas1,[1,iG1],gas2,[1,iG2],xAir);
    
    % COOL main stream (gas-liquid)
    fluidC1.state(1,1).T = CT1.A(1).T; fluidC1.state(1,1).p = CT1.A(1).p;
    [fluidC1] = update(fluidC1,[1,1],1);
    [fluidC1,gas1,~,iG1] = hex_TQ_cond(fluidC1,[1,1],gas1,[1,iG1],eff,1.0,ploss/2,'hex', 0, 0);
    
    % COOL main stream (gas-liquid)
    fluidC2.state(1,1).T = CT2.A(1).T; fluidC2.state(1,1).p = CT2.A(1).p;
    [fluidC2] = update(fluidC2,[1,1],1);
    [fluidC2,gas1,~,iG1] = hex_TQ_cond(fluidC2,[1,1],gas1,[1,iG1],eff,1.0,ploss/2,'hex', 0, 0);
    iSup = iG1;
    
    % REGENERATE (third stream cooling second stream)
    gas3.state(1,1).mdot = Qual*gas1.state(1,1).mdot;
    [gas2,gas3,iG2,iG3] = hex_TQ_cond(gas2,[1,iG2],gas3,[1,iG3],eff,0,(1-(1-ploss/2)^2),'regen',0,0);
    
    % MIX (main and second stream)
    [gas1,gas2,iG1,iG2] = mix_streams(gas1,[1,iG1],gas2,[1,iG2]);
    
    % EXPAND main stream
    p_aim = pLA;
    [gas1,iG1] = compexp(gas1,[1,iG1],eta,p_aim,4); %use isentropic efficiency
    
    % SEPARATE gas and liquid phases of main stream
    [gas1,gas3,iG1] = separate_phases(gas1,[1,iG1],gas3,[1,1]);
    
    %         mdot3 = gas3.state(1,iG3).mdot;
    %         Dh3   = gas1.state(1,1).h - gas3.state(1,1).h;
    %         Dh2   = gas2.state(1,1).h - gas1.state(1,iSup).h;
    %         mdot2 = mdot3*Dh3/Dh2;
    
    Qual  = gas3.state(1,1).mdot/mdot;
    if abs((Qual-Qual_0)/Qual) < 1e-6
        break
    else
        fprintf(1,'Qual = %.5f\n',Qual)
        Qual   = 0.8*Qual+0.2*Qual_0;
        Qual_0 = Qual;
        iG1 = iG1_0; iG2 = iG2_0; iG3 = iG3_0;
        iH  = iH_0;  iC1 = iC1_0; iC2 = iC2_0; iE = iE_0;
    end
end

fprintf(1,'\n %10s %10s %10s %10s %10s %10s %10s %10s','T(K)','p(bar)','h(kJ/kg)','s(kJ/kg/K)','m(kg/s)','Q','Stage','Ind');
for i0=1:iG1, fprintf(1,'\n %10.1f %10.1f %10.1f %10.2f %10.1f %10.1f %10s %10d',gas1.state(1,i0).T,gas1.state(1,i0).p/1e5,gas1.state(1,i0).h/1e3,gas1.state(1,i0).s/1e3,gas1.state(1,i0).mdot,gas1.state(1,i0).Q,gas1.stage(1,i0).type,i0); end; fprintf(1,'\n');
for i0=1:iG2, fprintf(1,'\n %10.1f %10.1f %10.1f %10.2f %10.1f %10.1f %10s %10d',gas2.state(1,i0).T,gas2.state(1,i0).p/1e5,gas2.state(1,i0).h/1e3,gas2.state(1,i0).s/1e3,gas2.state(1,i0).mdot,gas2.state(1,i0).Q,gas2.stage(1,i0).type,i0); end; fprintf(1,'\n');
for i0=1:iG3, fprintf(1,'\n %10.1f %10.1f %10.1f %10.2f %10.1f %10.1f %10s %10d',gas3.state(1,i0).T,gas3.state(1,i0).p/1e5,gas3.state(1,i0).h/1e3,gas3.state(1,i0).s/1e3,gas3.state(1,i0).mdot,gas3.state(1,i0).Q,gas3.stage(1,i0).type,i0); end; fprintf(1,'\n');

% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of source tanks to stablish initial mass
% Hot tanks
[HT.A(1)] = liquid_tank_start(fluidH,[1,1],HT.A(1),t_ch,T0);
[HT.A(2),HT.B(1),HT.B(2),WL_mixH_ch] = liquid_tanks_compute(fluidH,1,2,HT.A(1),t_ch,T0);

% Cold tanks
[CT1.A(1)] = liquid_tank_start(fluidC1,[1,1],CT1.A(1),t_ch,T0);
[CT1.A(2),CT1.B(1),CT1.B(2),WL_mixC1_ch] = liquid_tanks_compute(fluidC1,1,2,CT1.A(1),t_ch,T0);
[CT2.A(1)] = liquid_tank_start(fluidC2,[1,1],CT2.A(1),t_ch,T0);
[CT2.A(2),CT2.B(1),CT2.B(2),WL_mixC2_ch] = liquid_tanks_compute(fluidC2,1,2,CT2.A(1),t_ch,T0);
WL_mixC_ch = WL_mixC1_ch + WL_mixC2_ch;

% Liquid air tank
LA(2).T = gas1.state(1,iG1).T;
LA(2).p = gas1.state(1,iG1).p;
LA(2).M = gas1.state(1,iG1).mdot*t_ch;
LA(2).V = 1/gas1.state(1,iG1).rho*LA(2).M;
LA(2).H = gas1.state(1,iG1).h*LA(2).M;
LA(2).S = gas1.state(1,iG1).s*LA(2).M;
LA(2).B = LA(2).H - T0*LA(2).S;
LA(1).T = LA(2).T;
LA(1).p = LA(2).p;
LA(1).M = 0;
LA(1).V = 0;
LA(1).H = 0;
LA(1).S = 0;
LA(1).B = 0;

% Amospheric air source/sink
AA(1).T = gas1.state(1,1).T;
AA(1).p = gas1.state(1,1).p;
AA(1).M = (gas1.state(1,1).mdot - gas3.state(1,iG3).mdot)*t_ch;
AA(1).V = 1/gas1.state(1,1).rho*AA(1).M;
AA(1).H = gas1.state(1,1).h*AA(1).M;
AA(1).S = gas1.state(1,1).s*AA(1).M;
AA(1).B = AA(1).H - T0*AA(1).S;
AA(2).T = AA(1).T;
AA(2).p = AA(1).p;
AA(2).M = 0;
AA(2).V = 0;
AA(2).H = 0;
AA(2).S = 0;
AA(2).B = 0;