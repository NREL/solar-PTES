% First and second law balances

% Calculate work/heat/irreversibility terms for each component
% Charging compressors and discharging expanders
for ii = 1:length(CCMP)
    CCMP(ii) = compexp_energy(CCMP(ii),Load.time)  ;
end
% Charging expanders and discharging compressors
for ii = 1:length(CEXP)
    CEXP(ii) = compexp_energy(CEXP(ii),Load.time)  ;
end

switch Load.mode
    case {0,1,2,4,5,6}
        NC = Ne_ch ;
        NE = Nc_ch ;
    case 3
        NC = 3 ;
        NE = 3 ;
end

for ii = 1:length(DCMP)
    DCMP(ii) = compexp_energy(DCMP(ii),Load.time)  ;
end
for ii = 1:length(DEXP)
    DEXP(ii) = compexp_energy(DEXP(ii),Load.time)  ;
end

% Recompressor if specified for sCO2 cycle
if Load.mode == 4 || Load.mode == 5 || Load.mode == 6
    if Lrcmp
        RCMP = compexp_energy(RCMP,Load.time)  ;
    end
end

% FANS
for ii = 1:length(CFAN)
    CFAN(ii) = compexp_energy(CFAN(ii),Load.time)  ;
end
for ii = 1:length(DFAN)
    DFAN(ii) = compexp_energy(DFAN(ii),Load.time)  ;
end

% Adds up the contributions of the several stages to compute an energy
% balance
W_in_chg   = 0;
DH_chg     = 0;
W_lost_chg = 0;
QH_chg     = 0;  % heat to hot tanks
QH_dis     = 0;  % heat from hot tanks
QE_chg     = 0;  % heat rejected to environment

W_out_dis  = 0;
DH_dis     = 0;
W_lost_dis = 0;
QC_chg     = 0;  % heat from cold tanks
QC_dis     = 0;  % heat to cold tanks
QE_dis     = 0;  % heat rejected to environment

W_out_disRC  = 0;
QC_disRC     = 0;  
QH_disRC     = 0; 

nH = numel(fluidH);
nC = numel(fluidC);
for iL=1:Load.num
    
    if any(strcmp(Load.type(iL),{'chg','chgCO2','chgTSCO2'}))
        W_in_chg   = W_in_chg   -    sum([gas.stage(iL,:).w]   .*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        W_lost_chg = W_lost_chg + T0*sum([gas.stage(iL,:).sirr].*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        DH_chg     = DH_chg     +    sum([gas.stage(iL,:).Dh]  .*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        for i=1:nH
            QH_chg = QH_chg + sum([fluidH(i).stage(iL,:).q].*[fluidH(i).state(iL,1:(end-1)).mdot]*Load.time(iL));
        end
        for i=1:nC
            QC_chg = QC_chg - sum([fluidC(i).stage(iL,:).q].*[fluidC(i).state(iL,1:(end-1)).mdot]*Load.time(iL));
        end
        QE_chg = QE_chg + sum([environ.sink(iL,:).DHdot]*Load.time(iL));
        
    elseif any(strcmp(Load.type(iL),{'dis','disCO2','rcmpCO2','disTSCO2'}))
        W_out_dis  = W_out_dis  +    sum([gas.stage(iL,:).w]   .*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        W_lost_dis = W_lost_dis + T0*sum([gas.stage(iL,:).sirr].*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        DH_dis     = DH_dis     +    sum([gas.stage(iL,:).Dh]  .*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        for i=1:nH
            QH_dis = QH_dis - sum([fluidH(i).stage(iL,:).q].*[fluidH(i).state(iL,1:(end-1)).mdot]*Load.time(iL));
        end
        for i=1:nC
            QC_dis = QC_dis + sum([fluidC(i).stage(iL,:).q].*[fluidC(i).state(iL,1:(end-1)).mdot]*Load.time(iL));
        end
        QE_dis = QE_dis + sum([environ.sink(iL,:).DHdot]*Load.time(iL));
        
        % Also calculate the solar heat input to one tank
        if Load.mode == 6
            QH_sol = -sum([fluidH(1).stage(iL,:).q].*[fluidH(1).state(iL,1:(end-1)).mdot]*Load.time(iL));
            EX_sol = QH_sol - T0 * (fluidH(1).state(iL,1).s - fluidH(1).state(iL,2).s)*fluidH(1).state(iL,1).mdot*Load.time(iL);
        end
        
    elseif strcmp(Load.type(iL),'ran')
        W_out_dis  = W_out_dis  +    sum([steam.stage(iL,:).w]   .*[steam.state(iL,1:(end-1)).mdot]*Load.time(iL));
        W_lost_dis = W_lost_dis + T0*sum([steam.stage(iL,:).sirr].*[steam.state(iL,1:(end-1)).mdot]*Load.time(iL));
        DH_dis     = DH_dis     +    sum([steam.stage(iL,:).Dh]  .*[steam.state(iL,1:(end-1)).mdot]*Load.time(iL));
        for i=1:nH
            QH_dis = QH_dis - sum([fluidH(i).stage(iL,:).q].*[fluidH(i).state(iL,1:(end-1)).mdot]*Load.time(iL));
        end
        for i=1:nC
            QC_dis = QC_dis + sum([fluidC(i).stage(iL,:).q].*[fluidC(i).state(iL,1:(end-1)).mdot]*Load.time(iL));
        end
        QE_dis = QE_dis + sum([environ.sink(iL,:).DHdot]*Load.time(iL));
        
        % Also calculate heat and work terms for when there is only cooling from cold store
        if Load.options.useCold(iL) == 1
            W_out_disRC  = W_out_disRC  +    sum([steam.stage(iL,:).w]   .*[steam.state(iL,1:(end-1)).mdot]*Load.time(iL));
            for i=1:nH
                QH_disRC = QH_disRC - sum([fluidH(i).stage(iL,:).q].*[fluidH(i).state(iL,1:(end-1)).mdot]*Load.time(iL));
            end
            for i=1:nC
                QC_disRC = QC_disRC + sum([fluidC(i).stage(iL,:).q].*[fluidC(i).state(iL,1:(end-1)).mdot]*Load.time(iL));
            end
        end        
    end
    
    % Compute contributions from ambient air streams (heat rejection)
    if any(strcmp(Load.type(iL),{'chg','chgCO2','chgTSCO2'}))
        W_in_chg   = W_in_chg   -    sum([air.stage(iL,:).w]   .*[air.state(iL,1:(end-1)).mdot]*Load.time(iL));
        W_lost_chg = W_lost_chg + T0*sum([air.stage(iL,:).sirr].*[air.state(iL,1:(end-1)).mdot]*Load.time(iL));
        QE_chg     = QE_chg     +    sum([air.stage(iL,:).Dh]  .*[air.state(iL,1:(end-1)).mdot]*Load.time(iL));        
    elseif any(strcmp(Load.type(iL),{'dis','disCO2','ran','rcmpCO2','disTSCO2'}))
        W_out_dis  = W_out_dis  +    sum([air.stage(iL,:).w]   .*[air.state(iL,1:(end-1)).mdot]*Load.time(iL));
        W_lost_dis = W_lost_dis + T0*sum([air.stage(iL,:).sirr].*[air.state(iL,1:(end-1)).mdot]*Load.time(iL));
        QE_dis     = QE_dis     +    sum([air.stage(iL,:).Dh]  .*[air.state(iL,1:(end-1)).mdot]*Load.time(iL));      
    end
    
end

W_out_disNC = W_out_dis - W_out_disRC ;
QH_disNC = QH_dis - QH_disRC ;
QC_disNC = QC_dis - QC_disRC ;

fact    = 1e6*3600; % J to MWh
Net_chg = (+ W_in_chg  + QC_chg  - QH_chg  - DH_chg  - QE_chg);
Net_dis = (- W_out_dis - QC_dis  + QH_dis  - DH_dis  - QE_dis);

i_chg = Load.ind(any(Load.type == {'chg','chgCO2','chgTSCO2'},2));
i_dis = Load.ind(any(Load.type == {'dis','disCO2','ran','rcmpCO2','disTSCO2'},2));
i_ran = Load.ind(Load.type == 'ran');
i_gas = Load.ind(any(Load.type == {'chg','dis','chgCO2','disCO2','rcmpCO2','chgTSCO2','disTSCO2'},2));
i_act = Load.ind(any(Load.type == {'chg','dis','chgCO2','disCO2','ran','rcmpCO2','chgTSCO2','disTSCO2'},2));

t_chg = sum(Load.time(i_chg));
t_dis = sum(Load.time(i_dis));
ip1 = find(i_chg == 1,1,'first'); %index for printing (first charge period)
ip2 = find(i_dis == 1,1,'first'); %index for printing (first discharge period)

% Compute lost work on specific component types
% The WL_PTES_chg and WL_PTES_dis arrays are divided in 7 elements:
% 1: Losses in compressors
% 2: Losses in expanders
% 3: Losses in heat exchangers
% 4: Losses due to heat exchange with the environment
% 5: Mixing losses in storage tanks (liquids)
% 6: Mixing losses of the working fluid
% 7: Losses due to exergy leftover in tanks
WL_PTES_chg = zeros(1,7);
WL_PTES_dis = zeros(1,7);
for iL=1:Load.num
    if any(strcmp(Load.type(iL),{'chg','chgCO2','dis','disCO2','rcmpCO2','chgTSCO2','disTSCO2'}))
        for i0=1:gas.Nstg(iL)
            switch gas.stage(iL,i0).type
                case 'comp'
                    i1 = 1;
                case 'exp'
                    i1 = 2;
                case {'hex','regen'}
                    i1 = 3;
                case 'hex_reject'
                    i1 = 4;
                case 'mixing'
                    i1 = 6;
                case {'0','split','end'}
                otherwise
                    error('unknown stage type');
            end
            if any(strcmp(gas.stage(iL,i0).type,{'comp','exp','hex','regen','hex_reject','mixing'}))
                if any(strcmp(Load.type(iL),{'chg','chgCO2','chgTSCO2'}))
                    WL_PTES_chg(i1) = WL_PTES_chg(i1) + gas.stage(iL,i0).sirr*T0*gas.state(iL,i0).mdot*Load.time(iL);
                end
                if any(strcmp(Load.type(iL),{'dis','disCO2','rcmpCO2','disTSCO2'}))
                    WL_PTES_dis(i1) = WL_PTES_dis(i1) + gas.stage(iL,i0).sirr*T0*gas.state(iL,i0).mdot*Load.time(iL);
                end
            end
        end
    end
    if strcmp(Load.type(iL),'ran')
        for i0=1:steam.Nstg(iL)
            switch steam.stage(iL,i0).type
                case 'comp'
                    i1 = 1;
                case 'exp'
                    i1 = 2;
                case {'hex','regen'}
                    i1 = 3;
                case 'hex_reject'
                    i1 = 4;
                case 'mixing'
                    i1 = 6;
                case {'0','split','end'}
                otherwise
                    error('unknown stage type');
            end
            if any(strcmp(steam.stage(iL,i0).type,{'comp','exp','hex','regen','hex_reject','mixing'}))
                WL_PTES_dis(i1) = WL_PTES_dis(i1) + steam.stage(iL,i0).sirr*T0*steam.state(iL,i0).mdot*Load.time(iL);
            end
        end
    end
    for i0=1:air.Nstg(iL)
        switch air.stage(iL,i0).type
            case 'comp'
                i1 = 1;
            case 'exp'
                i1 = 2;
            case {'hex','regen'}
                i1 = 3;
            case 'hex_reject'
                i1 = 4;
            case 'mixing'
                i1 = 6;
            case {'0','split','end'}
            otherwise
                error('unknown stage type');
        end
        if any(strcmp(air.stage(iL,i0).type,{'comp','exp','hex','regen','hex_reject','mixing'}))
            if any(strcmp(Load.type(iL),{'chg','chgCO2','chgTSCO2'}))
                WL_PTES_chg(i1) = WL_PTES_chg(i1) + air.stage(iL,i0).sirr*T0*air.state(iL,i0).mdot*Load.time(iL);
            end
            if any(strcmp(Load.type(iL),{'dis','disCO2','ran','rcmpCO2','disTSCO2'}))
                WL_PTES_dis(i1) = WL_PTES_dis(i1) + air.stage(iL,i0).sirr*T0*air.state(iL,i0).mdot*Load.time(iL);
            end
        end
    end
    
end
switch Load.mode

    case {0,3,4,6} % PTES, Rankine, or sCO2-PTES
        
        for ii = 1 : Nhot
            WL_PTES_chg(5) = WL_PTES_chg(5) + HT(ii).WL_chg ;
            WL_PTES_dis(5) = WL_PTES_dis(5) + HT(ii).WL_dis ;
            WL_PTES_dis(7) = WL_PTES_dis(7) + HT(ii).A(end).B - HT(ii).A(1).B + HT(ii).B(end).B - HT(ii).B(1).B;
        end
        
        for ii = 1 : Ncld
            WL_PTES_chg(5) = WL_PTES_chg(5) + CT(ii).WL_chg ;
            WL_PTES_dis(5) = WL_PTES_dis(5) + CT(ii).WL_dis ;
            WL_PTES_dis(7) = WL_PTES_dis(7) + CT(ii).A(end).B - CT(ii).A(1).B + CT(ii).B(end).B - CT(ii).B(1).B;
        end
        
        WL_PTES_chg(5) = WL_PTES_chg(5) + AT.WL_chg ;
        WL_PTES_dis(5) = WL_PTES_dis(5) + AT.WL_dis ;

    case 1 % Heat pump only
        
        for ii = 1 : Nhot
            WL_PTES_chg(5) = WL_PTES_chg(5) + HT(ii).WL_chg ;
        end
           
        for ii = 1 : Ncld
            WL_PTES_chg(5) = WL_PTES_chg(5) + CT(ii).WL_chg ;
        end
        
        WL_PTES_chg(5) = WL_PTES_chg(5) + AT.WL_chg ;

    case {2,5} % Heat engine only
        WL_PTES_dis(5) = HT.WL_dis + AT.WL_dis ;
        WL_PTES_dis(7) = HT.A(end).B - HT.A(1).B + HT.B(end).B - HT.B(1).B;
end


% PRINT MAIN RESULTS ON SCREEN
if WM==1
    fprintf(1,'\n\n');
    fprintf(1,'PTES CYCLE\n');
    fprintf(1,'----------\n');
    
    % Print working fluid states
    fprintf(1,'Gas states:\n');
    for iL = i_gas
        print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
    end
    for iL = i_ran
        print_states(steam,iL,1:steam.Nstg(iL)+1,Load);
    end
    
    % Print hot streams
    fprintf(1,'\nHot fluid streams:\n');
    fprintf(1,'-->%s\n',fluidH(1).name);
    for iL = i_act
        for iH=1:numel(fluidH)
            print_states(fluidH(iH),iL,1:fluidH(iH).Nstg(iL)+1,Load)
        end
    end
    
    % Print cold streams
    fprintf(1,'\nCold fluid streams:\n');
    fprintf(1,'-->%s\n',fluidC(1).name);
    for iL = i_act
        for iC=1:numel(fluidC)
            print_states(fluidC(iC),iL,1:fluidC(iC).Nstg(iL)+1,Load)
        end
    end
    % Print hot tanks
    for ii = 1 : Nhot
        fprintf(1,'\nHot tank #%2i\n',ii);
        fprintf(1,'%10s %10s %13s %13s %13s %13s %8s ','A.T [K]','A.M [kg*1e6]','A.H [MWh]','B.T [K]','B.M [kg*1e6]','B.H [MWh]','state'); fprintf(1,'\n');
        for i0=1:(Load.num+1)
            fprintf(1,'%10.4g %13.3f %13.3f %10.4g %13.3f %13.3f %8d\n', HT(ii).A(i0).T,HT(ii).A(i0).M/1e6,HT(ii).A(i0).H/fact,HT(ii).B(i0).T,HT(ii).B(i0).M/1e6,HT(ii).B(i0).H/fact,i0)
        end
    end
    % Print cold tanks
    for ii = 1 : Ncld
        fprintf(1,'\nCold tank #%2i\n',ii);
        fprintf(1,'%10s %10s %13s %13s %13s %13s %8s ','A.T [K]','A.M [kg*1e6]','A.H [MWh]','B.T [K]','B.M [kg*1e6]','B.H [MWh]','state'); fprintf(1,'\n');
        for i0=1:(Load.num+1)
            fprintf(1,'%10.4g %13.3f %13.3f %10.4g %13.3f %13.3f %8d\n', CT(ii).A(i0).T,CT(ii).A(i0).M/1e6,CT(ii).A(i0).B/fact,CT(ii).B(i0).T,CT(ii).B(i0).M/1e6,CT(ii).B(i0).B/fact,i0)
        end
    end
    fprintf(1,'\n');
end

% Compute total loss from a second law perspective
WL_matrix = [ WL_PTES_chg ; WL_PTES_dis ; ];
Total_loss = sum(WL_matrix(:));

% Calculate tank stats
for ii = 1 : Nhot
    HT(ii) = tank_stats(HT(ii)) ;
end
for ii = 1 : Ncld
    CT(ii) = tank_stats(CT(ii)) ;
end

% Compute efficiencies, energy and power densities and errors
switch Load.mode
    case {0,3,4,6} % PTES, Rankine, or sCO2-PTES
        Heat_in_tanks = 0;
        for ii = 1 : Nhot
            Heat_in_tanks = Heat_in_tanks + (HT(ii).A(end).H - HT(ii).A(1).H) + (HT(ii).B(end).H - HT(ii).B(1).H) ;
        end        
        for ii = 1 : Ncld
            Heat_in_tanks = Heat_in_tanks + (CT(ii).A(end).H - CT(ii).A(1).H) + (CT(ii).B(end).H - CT(ii).B(1).H);
        end
        Heat_rejected = QE_chg + QE_dis;
        Total_Work_lost = W_in_chg - W_out_dis;
        First_law_error = (Heat_rejected + Heat_in_tanks - Total_Work_lost)/Total_Work_lost;
        Second_law_error = (Total_loss - Total_Work_lost)/Total_Work_lost;
        chi_PTES = W_out_dis/W_in_chg;
        vol = 0.0; % Calculate storage volume
        for ii = 1:Nhot; vol = vol + HT(ii).tank_volA + HT(ii).tank_volB ; end
        for ii = 1:Ncld; vol = vol + CT(ii).tank_volA + CT(ii).tank_volB ; end
        rhoE = W_out_dis/fact/vol*1e3; %kWh/m3
        
        % Calculate some special metrics for certain cycles
        if Load.mode == 3
            HEeff = 100. * W_out_dis / QH_dis ; % Heat engine average efficiency
            HEeffRC = 100. * W_out_disRC / QH_disRC ; % Efficiency for cycles where only cold tanks are used for condensing
            HEeffNC = 100. * W_out_disNC / QH_disNC ; % Efficiency for cycles where cold tanks are NOT used for condensing
        elseif Load.mode == 6
            SOLeff = 100. * W_out_dis / QH_sol ; % Solar conversion efficiency
            HEeff  = 100. * W_out_dis / QH_dis ; % Heat engine efficiency
            NETeff = 100. * (W_out_dis - W_in_chg) / QH_sol ; % Net efficiency
            EXeff  = 100. * W_out_dis / (W_in_chg + EX_sol) ; % Exergetic efficiency
        end
        
    case 1 % Heat pump only
        Heat_into_tanks   = (HT.A(end).H - HT.A(1).H) + (HT.B(end).H - HT.B(1).H) + (CT.A(end).H - CT.A(1).H) + (CT.B(end).H - CT.B(1).H);
        Exergy_into_tanks = (HT.A(end).B - HT.A(1).B) + (HT.B(end).B - HT.B(1).B) + (CT.A(end).B - CT.A(1).B) + (CT.B(end).B - CT.B(1).B);
        Heat_rejected = QE_chg;
        Total_Work_lost = W_in_chg - Exergy_into_tanks;
        First_law_error = (Heat_rejected + Heat_into_tanks - W_in_chg)/(W_in_chg);
        Second_law_error = (Total_loss - Total_Work_lost)/Total_Work_lost;
        chi_tot = Exergy_into_tanks/(W_in_chg);
        Exergy_into_hot_tanks = ((HT.A(end).B - HT.A(1).B) + (HT.B(end).B - HT.B(1).B));
        chi_hot = Exergy_into_hot_tanks/(W_in_chg);
        COP = QH_chg/W_in_chg;
        vol = 0.0; % Calculate storage volume
        for ii = 1:Nhot; vol = vol + HT(ii).tank_volA + HT(ii).tank_volB ; end
        rhoE = Exergy_into_hot_tanks/fact/vol*1e3; %kWh/m3
        
    case {2,5} % Heat engine only
        Heat_from_tanks   = (HT.A(2).H - HT.A(end).H) + (HT.B(2).H - HT.B(end).H);
        Exergy_from_tanks = (HT.A(2).B - HT.A(end).B) + (HT.B(2).B - HT.B(end).B);
        Heat_rejected = QE_dis;
        Total_Work_lost = Exergy_from_tanks - W_out_dis;
        First_law_error = (Heat_from_tanks - W_out_dis - Heat_rejected)/(W_out_dis);
        Second_law_error = (Total_loss - Total_Work_lost)/Total_Work_lost;
        chi_tot = (W_out_dis)/Exergy_from_tanks;
        EFF = W_out_dis/QH_dis;
        rhoE = W_out_dis/fact/(HT.A(end).V)*1e3; %kWh/m3
end


% Assume COLD rejection in heat pump mode!
if Load.mode == 1
    % Treat mixing_liquid loss and exergy of cold tanks as heat_reject
    WL_matrix(1,4) = WL_matrix(1,5) + (Exergy_into_tanks - Exergy_into_hot_tanks);
    WL_matrix(1,5) = 0;
elseif Load.mode == 5 || Load.mode == 6
    % Ignore tank losses from the solar tank
    WL_matrix(1,5) = WL_matrix(1,5) - HT(1).WL_chg ;
    WL_matrix(1,5) = WL_matrix(1,5) - HT(1).WL_dis ;
    WL_matrix(2,7) = WL_matrix(2,7) - (HT(1).A(end).B - HT(1).A(1).B + HT(1).B(end).B - HT(1).B(1).B);
end

switch Load.mode
    case {0,1,3,4,6}
        Exergy_in = W_in_chg;
    case {2,5}
        Exergy_in = Exergy_from_tanks;
end
WL_comp    = sum(WL_matrix(:,1))/Exergy_in*100;
WL_exp     = sum(WL_matrix(:,2))/Exergy_in*100;
WL_hexs    = sum(WL_matrix(:,3))/Exergy_in*100;
WL_reject  = sum(WL_matrix(:,4))/Exergy_in*100;
WL_mix_liq = sum(WL_matrix(:,5))/Exergy_in*100;
WL_mix_gas = sum(WL_matrix(:,6))/Exergy_in*100;
WL_tanks   = sum(WL_matrix(:,7))/Exergy_in*100;

if WM == 1
    fprintf(1,'\n\n');
    fprintf(1,'FIRST LAW BALANCE\n');
    fprintf(1,'-----------------\n');
    
    switch Load.mode
        case {0,1,3,4,6}
            fprintf(1,'CHARGE\n');
            fprintf(1,'Average power input:     %8.1f MW\n',W_in_chg/t_chg/1e6);
            fprintf(1,'Total charge time:       %8.1f h\n',t_chg/3600);
            fprintf(1,'Energy(el) input:        %8.1f MWh\n',W_in_chg/fact);
            fprintf(1,'Heat to hot tanks:       %8.1f MWh\n',QH_chg/fact);
            fprintf(1,'Heat from cold tanks:    %8.1f MWh\n',QC_chg/fact);
            fprintf(1,'DH working fluid:        %8.1f MWh\n',DH_chg/fact);
            fprintf(1,'Heat rejected:           %8.1f MWh\n',QE_chg/fact);
            fprintf(1,'NET:                     %8.1f MWh\n\n',Net_chg/fact);
    end
    
    switch Load.mode
        case {0,2,3,4,5,6}
            fprintf(1,'DISCHARGE\n');
            fprintf(1,'Average power output:    %8.1f MW\n',W_out_dis/t_dis/1e6);
            fprintf(1,'Total discharge time:    %8.1f h\n',t_dis/3600);
            fprintf(1,'Energy(el) output:       %8.1f MWh\n',W_out_dis/fact);
            fprintf(1,'Heat from hot tanks:     %8.1f MWh\n',QH_dis/fact);
            fprintf(1,'Heat to cold tanks:      %8.1f MWh\n',QC_dis/fact);
            fprintf(1,'DH working fluid:        %8.1f MWh\n',DH_dis/fact);
            fprintf(1,'Heat rejected:           %8.1f MWh\n',QE_dis/fact);
            fprintf(1,'NET:                     %8.1f MWh\n\n',Net_dis/fact);
    end
    
    switch Load.mode
        case 3
            fprintf(1,'DISCHARGE Efficiencies\n');
            fprintf(1,'Rankine cycle average efficiency:           %8.1f %%\n',HEeff);
            fprintf(1,'Rankine cycle efficiency NO cold stores:    %8.1f %%\n',HEeffNC);
            fprintf(1,'Rankine cycle efficiency using cold stores: %8.1f %%\n\n',HEeffRC);
        case 6
            fprintf(1,'EFFICIENCIES\n');
            fprintf(1,'Solar conversion efficiency:     %8.1f %%\n',SOLeff);
            fprintf(1,'Heat engine efficiency:          %8.1f %%\n',HEeff);
            fprintf(1,'Net efficiency:                  %8.1f %%\n',NETeff);
            fprintf(1,'Exergetic efficiency:            %8.1f %%\n\n',EXeff);
    end
    
    switch Load.mode
        case {0,3,4,6}
            fprintf(1,'COP:                     %8.2f \n',QH_chg/W_in_chg);
            fprintf(1,'Heat engine efficiency:  %8.2f %%\n',W_out_dis/QH_dis*100);
            fprintf(1,'Round trip efficiency:   %8.2f %%\n\n',chi_PTES*100);
            fprintf(1,'Exergy density:          %9.2f kWh/m3\n',rhoE);
            
            fprintf(1,'STORAGE MEDIA\n');
            for ii = 1 : Nhot
                fprintf(1,'%18s volume:%8.2f m3/MWh\n',fluidH(ii).name,HT(ii).fluid_volB/(W_out_dis/fact));
                fprintf(1,'%18s volume:%8.2f m3\n',fluidH(ii).name,HT(ii).fluid_volB);
                fprintf(1,'%18s mass:  %8.2f tons/MWh\n\n',fluidH(ii).name,HT(ii).fluid_mass/(W_out_dis/fact)/1e3);
            end
            for ii = 1 : Ncld
                fprintf(1,'%18s volume:%8.2f m3/MWh\n',fluidC(ii).name,CT(ii).fluid_volB/(W_out_dis/fact));
                fprintf(1,'%18s volume:%8.2f m3\n',fluidC(ii).name,CT(ii).fluid_volB);
                fprintf(1,'%18s mass:  %8.2f tons/MWh\n\n',fluidC(ii).name,CT(ii).fluid_mass/(W_out_dis/fact)/1e3);
            end
            
        case 1
            fprintf(1,'Exergetic efficiency:                %8.1f %%\n',chi_tot*100);
            fprintf(1,'Exergetic efficiency (cold reject):  %8.1f %%\n',chi_hot*100);
            fprintf(1,'Coefficient of Performance:           %8.2f\n\n',COP);
            fprintf(1,'Exergy density (hot tanks):          %9.2f kWh/m3\n',rhoE);
        case {2,5}
            fprintf(1,'Exergetic efficiency:      %7.1f %%\n',chi_tot*100);
            fprintf(1,'First law efficiency:      %7.1f %%\n\n',EFF*100);
            fprintf(1,'Exergy density:            %8.2f kWh/m3\n',rhoE);
    end
    
    fprintf(1,'First law error:   %8.5f %%\n',First_law_error*100);
    fprintf(1,'Second law error:  %8.5f %%\n\n',Second_law_error*100);
end

% Check for the case where storages are not returned to their original
% temperature. If stores are above T0 when discharged, then they must be
% returned to their original temp or greater. If stores are below T0 when
% discharged, then they must be returned to original temp or lower.
if Load.mode == any([0,4,6])
    problem = 0 ;
    
    % Hot tanks
    for ii = 1 : Nhot
       if HT(ii).A(1).T >= T0
           if HT(ii).A(4).T < HT(ii).A(1).T - 1.0; problem = 1 ; end
       else
           if HT(ii).A(4).T > HT(ii).A(1).T + 1.0; problem = 1 ; end
       end
    end
    if problem
        warning('Unsustainable discharge of a hot reservoir!')
    end
    problem = 0;
    % Cold tanks
    for ii = 1 : Ncld
       if CT(ii).A(1).T >= T0
           if CT(ii).A(4).T < CT(ii).A(1).T - 1.0; problem = 1 ; end
       else
           if CT(ii).A(4).T > CT(ii).A(1).T + 1.0; problem = 1 ; end
       end
    end
    if problem
        warning('Unsustainable discharge of a cold reservoir!')
    end
end

