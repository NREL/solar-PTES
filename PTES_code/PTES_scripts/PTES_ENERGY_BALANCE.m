% First and second law balances

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
nH = numel(fluidH);
nC = numel(fluidC);
for iL=1:Load.num
    
    % Calculate total work/heat/irreversibility terms
    % Charging compressors and discharging expanders
    for ii = 1 : Nc_ch
       CCMP(ii) = compexp_energy(CCMP(ii),Load.time(iL))  ;
    end
    % Charging expanders and discharging compressors
    for ii = 1 : Ne_ch
       CEXP(ii) = compexp_energy(CEXP(ii),Load.time(iL))  ;
    end
    
    switch Load.mode
        case {0,1,2,4}
            NC = Ne_ch ;
            NE = Nc_ch ;
        case 3
            NC = 3 ;
            NE = 3 ;
    end
    
    for ii = 1 : NC
        DCMP(ii) = compexp_energy(DCMP(ii),Load.time(iL))  ;
    end
    for ii = 1 : NE
        DEXP(ii) = compexp_energy(DEXP(ii),Load.time(iL))  ;
    end
    
    % Recompressor if specified for sCO2 cycle
    if Load.mode == 4
        if Lrcmp
            RCMP = compexp_energy(RCMP,Load.time(iL))  ;
        end
    end
    
    if any(strcmp(Load.type(iL),{'chg','chgCO2'}))
        W_in_chg   = W_in_chg   -    sum([gas.stage(iL,:).w]   .*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        W_lost_chg = W_lost_chg + T0*sum([gas.stage(iL,:).sirr].*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        DH_chg     = DH_chg     +    sum([gas.stage(iL,:).Dh]  .*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        for i=1:nH
            QH_chg = QH_chg + sum(fluidH(i).state(iL,1).mdot*(fluidH(i).state(iL,2).h-fluidH(i).state(iL,1).h)*Load.time(iL));
        end
        for i=1:nC
            QC_chg = QC_chg - sum(fluidC(i).state(iL,1).mdot*(fluidC(i).state(iL,2).h-fluidC(i).state(iL,1).h)*Load.time(iL));
        end
        QE_chg = QE_chg + sum([environ.sink(iL,:).DHdot]*Load.time(iL));
        
    elseif any(strcmp(Load.type(iL),{'dis','disCO2'}))
        W_out_dis  = W_out_dis  +    sum([gas.stage(iL,:).w]   .*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        W_lost_dis = W_lost_dis + T0*sum([gas.stage(iL,:).sirr].*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        DH_dis     = DH_dis     +    sum([gas.stage(iL,:).Dh]  .*[gas.state(iL,1:(end-1)).mdot]*Load.time(iL));
        for i=1:nH
            QH_dis = QH_dis - sum(fluidH(i).state(iL,1).mdot*(fluidH(i).state(iL,2).h-fluidH(i).state(iL,1).h)*Load.time(iL));
        end
        for i=1:nC
            QC_dis = QC_dis + sum(fluidC(i).state(iL,1).mdot*(fluidC(i).state(iL,2).h-fluidC(i).state(iL,1).h)*Load.time(iL));
        end
        QE_dis = QE_dis + sum([environ.sink(iL,:).DHdot]*Load.time(iL));
        
    elseif strcmp(Load.type(iL),'ran')
        W_out_dis  = W_out_dis  +    sum([steam.stage(iL,:).w]   .*[steam.state(iL,1:(end-1)).mdot]*Load.time(iL));
        W_lost_dis = W_lost_dis + T0*sum([steam.stage(iL,:).sirr].*[steam.state(iL,1:(end-1)).mdot]*Load.time(iL));
        DH_dis     = DH_dis     +    sum([steam.stage(iL,:).Dh]  .*[steam.state(iL,1:(end-1)).mdot]*Load.time(iL));
        for i=1:nH
            QH_dis = QH_dis - sum(fluidH(i).state(iL,1).mdot*(fluidH(i).state(iL,2).h-fluidH(i).state(iL,1).h)*Load.time(iL));
        end
        for i=1:nC
            QC_dis = QC_dis + sum(fluidC(i).state(iL,1).mdot*(fluidC(i).state(iL,2).h-fluidC(i).state(iL,1).h)*Load.time(iL));
        end
        QE_dis = QE_dis + sum([environ.sink(iL,:).DHdot]*Load.time(iL));
    end
end

fact    = 1e6*3600; % J to MWh
Net_chg = (+ W_in_chg  + QC_chg  - QH_chg  - DH_chg  - QE_chg);
Net_dis = (- W_out_dis - QC_dis  + QH_dis  - DH_dis  - QE_dis);

i_chg = Load.ind(any(Load.type == {'chg','chgCO2'},2));
i_dis = Load.ind(any(Load.type == {'dis','disCO2'},2));
i_ran = Load.ind(Load.type == 'ran');
i_gas = Load.ind(any(Load.type == {'chg','dis'},2));
i_act = Load.ind(any(Load.type == {'chg','dis','ran'},2));

t_chg = sum(Load.time(i_chg));
t_dis = sum(Load.time([i_dis,i_ran]));
ip1 = find(i_chg == 1,1,'first'); %index for printing (first charge period)
ip2 = find(i_dis == 1,1,'first'); %index for printing (first discharge period)

%Compute lost work on specific component types
WL_PTES_chg = zeros(1,7);
WL_PTES_dis = zeros(1,7);
for iL=1:Load.num
    if any(strcmp(Load.type(iL),{'chg','chgCO2'}))
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
            if any(strcmp(gas.stage(iL,i0).type,{'comp','exp','hex','regen','hex_reject'}))
                WL_PTES_chg(i1) = WL_PTES_chg(i1) + gas.stage(iL,i0).sirr*T0*gas.state(iL,i0).mdot*Load.time(iL);
            end
        end
    end
    if any(strcmp(Load.type(iL),{'dis','disCO2'}))
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
                WL_PTES_dis(i1) = WL_PTES_dis(i1) + gas.stage(iL,i0).sirr*T0*gas.state(iL,i0).mdot*Load.time(iL);
            end
        end
    end
    if strcmp(Load.type(iL),'ran')
        for is=1:numel(steam)
            for i0=1:steam(is).Nstg(iL)
                switch steam(is).stage(iL,i0).type
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
                if any(strcmp(steam(is).stage(iL,i0).type,{'comp','exp','hex','regen','hex_reject','mixing'}))
                    WL_PTES_dis(i1) = WL_PTES_dis(i1) + steam(is).stage(iL,i0).sirr*T0*steam(is).state(iL,i0).mdot*Load.time(iL);
                end
            end
        end
    end
end
switch Load.mode

    case {0,3,4} % PTES, Rankine, or sCO2-PTES

        WL_PTES_chg(5) = 0 ; WL_PTES_dis(5) = 0 ;
        WL_PTES_chg(7) = 0;  WL_PTES_dis(7) = 0;
        
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

    case 1 % Heat pump only
        WL_PTES_chg(5) = 0;
        WL_PTES_dis(5) = 0;
        WL_PTES_chg(7) = 0;
        WL_PTES_dis(7) = 0;
        
        for ii = 1 : Nhot
            WL_PTES_chg(5) = WL_PTES_chg(5) + HT(ii).WL_chg ;
        end
           
        for ii = 1 : Ncld
            WL_PTES_chg(5) = WL_PTES_chg(5) + CT(ii).WL_chg ;
        end

    case 2 % Heat engine only
        WL_PTES_chg(5)  = 0;
        WL_PTES_dis(5) = HT.WL_dis;
        WL_PTES_chg(7)  = 0;
        %WL_PTES_dis(7) = 0;
        WL_PTES_dis(7) = HT.A(end).B - HT.A(1).B + HT.B(end).B - HT.B(1).B;
end

% PRINT MAIN RESULTS ON SCREEN
if WM==1
    fprintf(1,'\n\n');
    fprintf(1,'PTES CYCLE\n');
    fprintf(1,'----------\n');
    fprintf(1,'Gas states:\n');
    fprintf(1,'%11s ','T [K]','p [bar]','h [MJ/kg]','s [kJ/kg/K]','mdot [kg/s]','Inlet of','Cycle'); fprintf(1,'\n');
    for iL = 1:Load.num
        for i0=1:gas.Nstg(iL)
            fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.1f %11s %11s\n',...
                gas.state(iL,i0).T, gas.state(iL,i0).p/1e5, gas.state(iL,i0).h/1e6,...
                gas.state(iL,i0).s/1e3, gas.state(iL,i0).mdot, gas.stage(iL,i0).type, Load.type(iL));
        end
        if any(strcmp(Load.type(iL),["chg","dis","chgCO2","disCO2"])), fprintf(1,'\n'); end
    end
    
    % Print hot streams
    for i0=1:nH
        fprintf(1,'\nHot fluid stream #%2i:\n',i0);
        fprintf(1,'-->%s\n',fluidH(i0).name);
        fprintf(1,'%11s ','Tin[K]','Tout[K]','Δh [MJ/kg]','Δs [kJ/kg/K]','mdot[kg/s]','Stream','Cycle'); fprintf(1,'\n');
        for iL = 1:Load.num
            if any(strcmp(Load.type(iL),["chg","dis","chgCO2","disCO2"]))
                fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.2f %11d %11s\n',...
                    fluidH(i0).state(iL,1).T, fluidH(i0).state(iL,2).T,...
                    (fluidH(i0).state(iL,2).h - fluidH(i0).state(iL,1).h)/1e6,...
                    (fluidH(i0).state(iL,2).s - fluidH(i0).state(iL,1).s)/1e3,...
                    fluidH(i0).state(iL,1).mdot, i0, Load.type(iL));
            end
        end
    end
    
    % Print cold streams
    for i0=1:nC
        fprintf(1,'\nCold fluid stream #%2i:\n',i0);
        fprintf(1,'-->%s\n',fluidC(1).name);
        fprintf(1,'%11s ','Tin[K]','Tout[K]','Δh [MJ/kg]','Δs [kJ/kg/K]','mdot[kg/s]','Stream','Cycle'); fprintf(1,'\n');
        for iL = 1:Load.num
            if any(strcmp(Load.type(iL),["chg","dis","chgCO2","disCO2"]))
                fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.2f %11d %11s\n',...
                    fluidC(i0).state(iL,1).T, fluidC(i0).state(iL,2).T,...
                    (fluidC(i0).state(iL,2).h - fluidC(i0).state(iL,1).h)/1e6,...
                    (fluidC(i0).state(iL,2).s - fluidC(i0).state(iL,1).s)/1e3,...
                    fluidC(i0).state(iL,1).mdot, i0, Load.type(iL));
            end
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

WL_matrix = [ WL_PTES_chg ; WL_PTES_dis ; ];
Total_loss = sum(WL_matrix(:));


% Compute efficiencies, energy and power densities and errors
switch Load.mode
    case {0,3,4} % PTES, Rankine, or sCO2-PTES
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
        for ii = 1:Nhot; vol = vol + HT(ii).A(1).V ; end
        for ii = 1:Ncld; vol = vol + CT(ii).A(1).V ; end
        rhoE = W_out_dis/fact/vol*1e3; %kWh/m3
        %rhoP_ch  = (W_in_chg/t_ch/(gas_min_rho_ch.mdot/gas_min_rho_ch.rho)/1e6); %MW/(m3/s)
        %rhoP_dis = (W_out_dis/t_dis/(gas_min_rho_dis.mdot/gas_min_rho_dis.rho)/1e6); %MW/(m3/s)
        
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
        for ii = 1:Nhot; vol = vol + HT(ii).A(1).V ; end
        rhoE = Exergy_into_hot_tanks/fact/vol*1e3; %kWh/m3
        %rhoP_ch  = (W_in_chg/t_ch/(gas_min_rho_ch.mdot/gas_min_rho_ch.rho)/1e6); %MW/(m3/s)
        
    case 2 % Heat engine only
        Heat_from_tanks   = (HT.A(2).H - HT.A(end).H) + (HT.B(2).H - HT.B(end).H);
        Exergy_from_tanks = (HT.A(2).B - HT.A(end).B) + (HT.B(2).B - HT.B(end).B);
        Heat_rejected = QE_dis;
        Total_Work_lost = Exergy_from_tanks - W_out_dis;
        First_law_error = (Heat_from_tanks - W_out_dis - Heat_rejected)/(W_out_dis);
        Second_law_error = (Total_loss - Total_Work_lost)/Total_Work_lost;
        chi_tot = (W_out_dis)/Exergy_from_tanks;
        EFF = W_out_dis/QH_dis;
        rhoE = W_out_dis/fact/(HT.A(end).V)*1e3; %kWh/m3
        %rhoP_dis = (W_out_dis/t_dis/(gas_min_rho_dis.mdot/gas_min_rho_dis.rho)/1e6); %MW/(m3/s)
end


% Assume COLD rejection in heat pump mode!
if Load.mode == 1
    % Treat mixing_liquid loss and exergy of cold tanks as heat_reject
    WL_matrix(1,4) = WL_matrix(1,5) + (Exergy_into_tanks - Exergy_into_hot_tanks);
    WL_matrix(1,5) = 0;
end

switch Load.mode
    case {0,1,3,4}
        Exergy_in = W_in_chg;
    case 2
        Exergy_in = Exergy_from_tanks;
end
WL_comp    = sum(WL_matrix(:,1))/Exergy_in*100;
WL_exp     = sum(WL_matrix(:,2))/Exergy_in*100;
WL_hexs    = sum(WL_matrix(:,3))/Exergy_in*100;
WL_reject  = sum(WL_matrix(:,4))/Exergy_in*100;
WL_mix_liq = sum(WL_matrix(:,5))/Exergy_in*100;
WL_mix_gas = sum(WL_matrix(:,6))/Exergy_in*100;
WL_tanks   = sum(WL_matrix(:,7))/Exergy_in*100;

%[WR,WR_dis]   = work_ratio(gas.stage,stages_ch,stages_dis);

%% THIS SEEMS TO BE A DUPLICATE?
%{
% PRINT MAIN RESULTS ON SCREEN
if WM==1
    fprintf(1,'\n\n');
    fprintf(1,'PTES CYCLE\n');
    fprintf(1,'----------\n');
    fprintf(1,'Gas states:\n');
    for iL = i_gas
        print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
    end
    for iL = i_ran
        print_states(steam,iL,1:steam.Nstg(iL)+1,Load);
    end
    
    fprintf(1,'\nHot fluid streams:\n');
    fprintf(1,'-->%s\n',fluidH(1).name);
    for iL = i_act
        for iH=1:numel(fluidH)
            print_states(fluidH(iH),iL,1:2,Load)
        end
    end
    fprintf(1,'\nCold fluid streams:\n');
    fprintf(1,'-->%s\n',fluidC(1).name);
    for iL = i_act
        for iC=1:numel(fluidC)
            print_states(fluidC(iC),iL,1:2,Load)
        end
    end
    fprintf(1,'\nHot tanks\n');
    fprintf(1,'%10s %10s %13s %13s %13s %13s %8s ','A.T [K]','A.M [kg*1e6]','A.H [MWh]','B.T [K]','B.M [kg*1e6]','B.H [MWh]','state'); fprintf(1,'\n');
    for i0=1:(Load.num+1)
        fprintf(1,'%10.4g %13.3f %13.3f %10.4g %13.3f %13.3f %8d\n', HT.A(i0).T,HT.A(i0).M/1e6,HT.A(i0).H/fact,HT.B(i0).T,HT.B(i0).M/1e6,HT.B(i0).H/fact,i0)
    end    
    fprintf(1,'\nCold tanks\n');
    fprintf(1,'%10s %10s %13s %13s %13s %13s %8s ','A.T [K]','A.M [kg*1e6]','A.H [MWh]','B.T [K]','B.M [kg*1e6]','B.H [MWh]','state'); fprintf(1,'\n');
    for i0=1:(Load.num+1)
        fprintf(1,'%10.4g %13.3f %13.3f %10.4g %13.3f %13.3f %8d\n', CT.A(i0).T,CT.A(i0).M/1e6,CT.A(i0).B/fact,CT.B(i0).T,CT.B(i0).M/1e6,CT.B(i0).B/fact,i0)
    end    
    fprintf(1,'\n');
end
%}
%%
if WM == 1
    fprintf(1,'\n\n');
    fprintf(1,'FIRST LAW BALANCE\n');
    fprintf(1,'-----------------\n');
    
    switch Load.mode
        case {0,1,3,4}
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
        case {0,2,3,4}
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
        case {0,3,4}
            fprintf(1,'Round trip efficiency:   %8.1f %%\n\n',chi_PTES*100);
            fprintf(1,'Exergy density:          %9.2f kWh/m3\n',rhoE);
            %fprintf(1,'Power density (charge):  %9.2f MW/(m3/s)\n',rhoP_ch);
            %fprintf(1,'Power density (disch):   %9.2f MW/(m3/s)\n\n',rhoP_dis);
            
            fprintf(1,'STORAGE MEDIA\n');
            for ii = 1 : Nhot
                fprintf(1,'%18s volume:%8.2f m3/MWh\n',fluidH(ii).name,HT(ii).A(1).V/(W_out_dis/fact));
                fprintf(1,'%18s mass:  %8.2f tons/MWh\n\n',fluidH(ii).name,HT(ii).A(1).M/(W_out_dis/fact)/1e3);
            end
            for ii = 1 : Ncld
                fprintf(1,'%18s volume:%8.2f m3/MWh\n',fluidC(ii).name,CT(ii).A(1).V/(W_out_dis/fact));
                fprintf(1,'%18s mass:  %8.2f tons/MWh\n\n',fluidC(ii).name,CT(ii).A(1).M/(W_out_dis/fact)/1e3);
            end
            
            %fprintf(1,'Work Ratios:  %6.3g (charge) %6.3g (discharge)\n\n',WR,WR_dis);
        case 1
            fprintf(1,'Exergetic efficiency:                %8.1f %%\n',chi_tot*100);
            fprintf(1,'Exergetic efficiency (cold reject):  %8.1f %%\n',chi_hot*100);
            fprintf(1,'Coefficient of Performance:           %8.2f\n\n',COP);
            fprintf(1,'Exergy density (hot tanks):          %9.2f kWh/m3\n',rhoE);
            %fprintf(1,'Power density (charge):              %9.2f MW/(m3/s)\n\n',rhoP_ch);
        case 2
            fprintf(1,'Exergetic efficiency:      %7.1f %%\n',chi_tot*100);
            fprintf(1,'First law efficiency:      %7.1f %%\n\n',EFF*100);
            fprintf(1,'Exergy density:            %8.2f kWh/m3\n',rhoE);
            %fprintf(1,'Power density (discharge): %8.2f MW/(m3/s)\n\n',rhoP_dis);
    end
    
    fprintf(1,'First law error:   %8.5f %%\n',First_law_error*100);
    fprintf(1,'Second law error:  %8.5f %%\n\n',Second_law_error*100);
end

% Check for the case where storages are not returned to their original
% temperature. If stores are above T0 when discharged, then they must be
% returned to their original temp or greater. If stores are below T0 when
% discharged, then they must be returned to original temp or lower.
if (Load.mode == 0 || Load.mode == 4)
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

