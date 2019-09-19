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
    if strcmp(Load.type(iL),'chg')
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
    elseif strcmp(Load.type(iL),'dis')
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
    end
end

fact    = 1e6*3600; % J to MWh
Net_chg = (+ W_in_chg  + QC_chg  - QH_chg  - DH_chg  - QE_chg);
Net_dis = (- W_out_dis - QC_dis  + QH_dis  - DH_dis  - QE_dis);
i_chg = Load.type == 'chg';
i_dis = Load.type == 'dis';
t_chg = sum(Load.time(i_chg));
t_dis = sum(Load.time(i_dis));
ip1 = find(i_chg == 1,1,'first'); %index for printing (first charge period)
ip2 = find(i_dis == 1,1,'first'); %index for printing (first discharge period)


%Compute lost work on specific component types
WL_PTES_chg = zeros(1,6);
WL_PTES_dis = zeros(1,6);
for iL=1:Load.num
    if strcmp(Load.type(iL),'chg')
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
                case '0'
                otherwise
                    error('unknown stage type');
            end
            if any(strcmp(gas.stage(iL,i0).type,{'comp','exp','hex','regen','hex_reject'}))
                WL_PTES_chg(i1) = WL_PTES_chg(i1) + gas.stage(iL,i0).sirr*T0*gas.state(iL,i0).mdot*Load.time(iL);
            end
        end
    end
    if strcmp(Load.type(iL),'dis')
        for i0=1:gas.Nstg(iL)
            switch gas.stage(iL,i0).type
                case 'comp'
                    i1 = 1;
                case 'exp'
                    i1 = 2;
                case {'hex','regen','mixing'}
                    i1 = 3;
                case 'hex_reject'
                    i1 = 4;
                case '0'
                otherwise
                    error('unknown stage type');
            end
            if any(strcmp(gas.stage(iL,i0).type,{'comp','exp','hex','regen','hex_reject'}))
                WL_PTES_dis(i1) = WL_PTES_dis(i1) + gas.stage(iL,i0).sirr*T0*gas.state(iL,i0).mdot*Load.time(iL);
            end
        end
    end
end
switch Load.mode
    case 0 % PTES
        WL_PTES_chg(5)  = HT.WL_chg + CT.WL_chg;
        WL_PTES_dis(5) = HT.WL_dis + CT.WL_dis;
        WL_PTES_chg(6)  = 0;
        WL_PTES_dis(6) = HT.A(end).B - HT.A(1).B + HT.B(end).B - HT.B(1).B + CT.A(end).B - CT.A(1).B + CT.B(end).B - CT.B(1).B;
    case 1 % Heat pump only
        WL_PTES_chg(5)  = HT.WL_chg + CT.WL_chg;
        WL_PTES_dis(5) = 0;
        WL_PTES_chg(6)  = 0;
        WL_PTES_dis(6) = 0;
    case 2 % Heat engine only
        WL_PTES_chg(5)  = 0;
        WL_PTES_dis(5) = HT.WL_dis;
        WL_PTES_chg(6)  = 0;
        %WL_PTES_dis(6) = 0;
        WL_PTES_dis(6) = HT.A(end).B - HT.A(1).B + HT.B(end).B - HT.B(1).B;
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
        if any(strcmp(Load.type(iL),["chg","dis"])), fprintf(1,'\n'); end
    end
    fprintf(1,'\nHot fluid streams:\n');
    fprintf(1,'-->%s\n',fluidH(1).name);
    fprintf(1,'%11s ','Tin[K]','Tout[K]','Δh [MJ/kg]','Δs [kJ/kg/K]','mdot[kg/s]','Stream','Cycle'); fprintf(1,'\n');
    for iL = 1:Load.num
        if any(strcmp(Load.type(iL),["chg","dis"]))
            for i0=1:nH
                fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.2f %11d %11s\n',...
                    fluidH(i0).state(iL,1).T, fluidH(i0).state(iL,2).T,...
                    (fluidH(i0).state(iL,2).h - fluidH(i0).state(iL,1).h)/1e6,...
                    (fluidH(i0).state(iL,2).s - fluidH(i0).state(iL,1).s)/1e3,...
                    fluidH(i0).state(iL,1).mdot, i0, Load.type(iL));
            end
        end
    end
    fprintf(1,'\nCold fluid streams:\n');
    fprintf(1,'-->%s\n',fluidC(1).name);
    fprintf(1,'%11s ','Tin[K]','Tout[K]','Δh [MJ/kg]','Δs [kJ/kg/K]','mdot[kg/s]','Stream','Cycle'); fprintf(1,'\n');
    for iL = 1:Load.num
        if any(strcmp(Load.type(iL),["chg","dis"]))
            for i0=1:nC
                fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.2f %11d %11s\n',...
                    fluidC(i0).state(iL,1).T, fluidC(i0).state(iL,2).T,...
                    (fluidC(i0).state(iL,2).h - fluidC(i0).state(iL,1).h)/1e6,...
                    (fluidC(i0).state(iL,2).s - fluidC(i0).state(iL,1).s)/1e3,...
                    fluidC(i0).state(iL,1).mdot, i0, Load.type(iL));
            end
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

WL_matrix = [ WL_PTES_chg ; WL_PTES_dis ; ];
Total_loss = sum(WL_matrix(:));


% Compute efficiencies, energy and power densities and errors
switch Load.mode
    case 0 % PTES
        Heat_in_tanks   = (HT.A(end).H - HT.A(1).H) + (HT.B(end).H - HT.B(1).H) + (CT.A(end).H - CT.A(1).H) + (CT.B(end).H - CT.B(1).H);
        Heat_rejected = QE_chg + QE_dis;
        Total_Work_lost = W_in_chg - W_out_dis;
        First_law_error = (Heat_rejected + Heat_in_tanks - Total_Work_lost)/Total_Work_lost;
        Second_law_error = (Total_loss - Total_Work_lost)/Total_Work_lost;
        chi_PTES = W_out_dis/W_in_chg;
        rhoE = W_out_dis/fact/(HT.A(1).V + CT.A(1).V)*1e3; %kWh/m3
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
        rhoE = Exergy_into_hot_tanks/fact/(HT.A(1).V)*1e3; %kWh/m3
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
    case {0,1}
        Exergy_in = W_in_chg;
    case 2
        Exergy_in = Exergy_from_tanks;
end
WL_comp    = sum(WL_matrix(:,1))/Exergy_in*100;
WL_exp     = sum(WL_matrix(:,2))/Exergy_in*100;
WL_hexs    = sum(WL_matrix(:,3))/Exergy_in*100;
WL_reject  = sum(WL_matrix(:,4))/Exergy_in*100;
WL_mix_liq = sum(WL_matrix(:,5))/Exergy_in*100;
WL_tanks   = sum(WL_matrix(:,6))/Exergy_in*100;

%[WR,WR_dis]   = work_ratio(gas.stage,stages_ch,stages_dis);

if WM == 1
    fprintf(1,'\n\n');
    fprintf(1,'FIRST LAW BALANCE\n');
    fprintf(1,'-----------------\n');
    
    switch Load.mode
        case {0,1}
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
        case {0,2}
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
        case 0 
            fprintf(1,'Round trip efficiency:   %8.1f %%\n\n',chi_PTES*100);
            fprintf(1,'Exergy density:          %9.2f kWh/m3\n',rhoE);
            %fprintf(1,'Power density (charge):  %9.2f MW/(m3/s)\n',rhoP_ch);
            %fprintf(1,'Power density (disch):   %9.2f MW/(m3/s)\n\n',rhoP_dis);
            
            fprintf(1,'STORAGE MEDIA\n');
            fprintf(1,'%18s volume:%8.2f m3/MWh\n',fluidH(1).name,HT.A(1).V/(W_out_dis/fact));
            fprintf(1,'%18s volume:%8.2f m3/MWh\n\n',fluidC(1).name,CT.A(1).V/(W_out_dis/fact));
            
            fprintf(1,'%18s mass:  %8.2f tons/MWh\n',fluidH(1).name,HT.A(1).M/(W_out_dis/fact)/1e3);
            fprintf(1,'%18s mass:  %8.2f tons/MWh\n',fluidC(1).name,CT.A(1).M/(W_out_dis/fact)/1e3);
            
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


% Manage warnings
if  Load.mode == 0 && ((CT.A(1).T > T0) && (CT.A(4).T < (CT.A(1).T - 1)))
    warning('Unsustainable discharge of cold reservoir!')
end
