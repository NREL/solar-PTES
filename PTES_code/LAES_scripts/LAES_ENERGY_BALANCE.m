% First and second law balances

% Select mode
switch mode
    case 0 % PTES
    case 1 % Heat pump only
        t_dis = 0;
    case 2 % Heat engine only
end


% Adds up the contributions of the several stages to compute an energy
% balance
W_in_ch   = 0;
DH_ch     = 0;
W_lost_ch = 0;
W_abs_ch  = 0;

W_out_dis  = 0;
DH_dis     = 0;
W_lost_dis = 0;
W_abs_dis  = 0;

stages_ch = find(~strcmp({gas.stage(1,:).type},'0'),1,'last')-1;
stages_dis = find(~strcmp({gas.stage(2,:).type},'0'),1,'last')-1;
n1 = max([stages_ch,stages_dis]);
for i=1:n1
    %Work
    W_in_ch   = W_in_ch   - gas.stage(1,i).w*gas.state(1,i).mdot*t_ch;    
    W_out_dis = W_out_dis + gas.stage(2,i).w*gas.state(2,i).mdot*t_dis;
    W_abs_ch  = W_abs_ch  + abs(gas.stage(1,i).w*gas.state(1,i).mdot*t_ch);
    W_abs_dis = W_abs_dis + abs(gas.stage(2,i).w*gas.state(2,i).mdot*t_dis);
    
    %Lost Work
    W_lost_ch  = W_lost_ch  + gas.stage(1,i).sirr*T0*gas.state(1,i).mdot*t_ch;
    W_lost_dis = W_lost_dis + gas.stage(2,i).sirr*T0*gas.state(2,i).mdot*t_dis;
    
    %Enthalpy change
    DH_ch  = DH_ch + gas.stage(1,i).Dh*gas.state(1,i).mdot*t_ch;
    DH_dis = DH_dis + gas.stage(2,i).Dh*gas.state(2,i).mdot*t_dis;
end

QH_ch = 0;   % heat to hot tanks
QH_dis = 0;  % heat from hot tanks
n2 = numel(fluidH);
for i=1:n2
    QH_ch  = QH_ch  + fluidH(i).state(1,1).mdot*(fluidH(i).state(1,2).h-fluidH(i).state(1,1).h)*t_ch;
    QH_dis = QH_dis - fluidH(i).state(2,1).mdot*(fluidH(i).state(2,2).h-fluidH(i).state(2,1).h)*t_dis;
end

QC_ch = 0;   % heat from cold tanks
QC_dis = 0;  % heat to cold tanks
n3 = numel(fluidC);
for i=1:n3
    QC_ch  = QC_ch  - fluidC(i).state(1,1).mdot*(fluidC(i).state(1,2).h-fluidC(i).state(1,1).h)*t_ch;
    QC_dis = QC_dis + fluidC(i).state(2,1).mdot*(fluidC(i).state(2,2).h-fluidC(i).state(2,1).h)*t_dis;
end

QExt_ch  = 0;  % heat rejected to environment
QExt_dis = 0;  % heat rejected to environment
[~,n4] = size(environ.sink);
for i=1:n4
    QExt_ch   = QExt_ch  + environ.sink(1,i).DHdot*t_ch;
    QExt_dis  = QExt_dis + environ.sink(2,i).DHdot*t_dis;
end

%Compute lost work on specific component types
WL_PTES_ch  = zeros(1,6);
WL_PTES_dis = zeros(1,6);
for i0=1:n1
    switch gas.stage(1,i0).type
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
    if any(strcmp(gas.stage(1,i0).type,{'comp','exp','hex','regen','hex_reject'}))
        WL_PTES_ch(i1) = WL_PTES_ch(i1) + gas.stage(1,i0).sirr*T0*gas.state(1,i0).mdot*t_ch;
    end
    switch gas.stage(2,i0).type
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
    if any(strcmp(gas.stage(2,i0).type,{'comp','exp','hex','regen','hex_reject'}))
        WL_PTES_dis(i1) = WL_PTES_dis(i1) + gas.stage(2,i0).sirr*T0*gas.state(2,i0).mdot*t_dis;
    end
end
switch mode
    case 0 % PTES
        WL_PTES_ch(5)  = WL_mixH_ch + WL_mixC_ch;
        WL_PTES_dis(5) = WL_mixH_dis + WL_mixC_dis;
        WL_PTES_ch(6)  = 0;
        WL_PTES_dis(6) = HT.A(4).B - HT.A(1).B + HT.B(4).B - HT.B(1).B + CT.A(4).B - CT.A(1).B + CT.B(4).B - CT.B(1).B;
    case 1 % Heat pump only
        WL_PTES_ch(5)  = WL_mixH_ch + WL_mixC_ch;
        WL_PTES_dis(5) = 0;
        WL_PTES_ch(6)  = 0;
        WL_PTES_dis(6) = 0;
    case 2 % Heat engine only
        WL_PTES_ch(5)  = 0;
        WL_PTES_dis(5) = WL_mixH_dis;
        WL_PTES_ch(6)  = 0;
        WL_PTES_dis(6) = 0;
end



% PRINT MAIN RESULTS ON SCREEN
if WM==1
    fprintf(1,'\n\n');
    fprintf(1,'PTES CYCLE\n');
    fprintf(1,'----------\n');
    fprintf(1,'Gas states:\n');
    fprintf(1,'%11s ','T [K]','p [bar]','h [MJ/kg]','s [kJ/kg/K]','mdot [kg/s]','Inlet of','Cycle'); fprintf(1,'\n');
    for i0=1:stages_ch+1,  fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.1f %11s %11s\n', gas.state(1,i0).T, gas.state(1,i0).p/1e5, gas.state(1,i0).h/1e6, gas.state(1,i0).s/1e3, gas.state(1,i0).mdot, gas.stage(1,i0).type, 'ch');end
    fprintf(1,'\n');
    for i0=1:stages_dis+1, fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.1f %11s %11s\n', gas.state(2,i0).T, gas.state(2,i0).p/1e5, gas.state(2,i0).h/1e6, gas.state(2,i0).s/1e3, gas.state(2,i0).mdot, gas.stage(2,i0).type, 'dis');end
    fprintf(1,'\nHot fluid streams:\n');
    fprintf(1,'-->%s\n',fluidH(1).name);
    fprintf(1,'%11s ','Tin[K]','Tout[K]','Δh [MJ/kg]','Δs [kJ/kg/K]','mdot[kg/s]','Stream','Cycle'); fprintf(1,'\n');
    for i0=1:n2, fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.2f %11d %11s\n', fluidH(i0).state(1,1).T, fluidH(i0).state(1,2).T, (fluidH(i0).state(1,2).h - fluidH(i0).state(1,1).h)/1e6, (fluidH(i0).state(1,2).s - fluidH(i0).state(1,1).s)/1e3, fluidH(i0).state(1,1).mdot, i0, 'ch');end
    for i0=1:n2, fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.2f %11d %11s\n', fluidH(i0).state(2,1).T, fluidH(i0).state(2,2).T, (fluidH(i0).state(2,2).h - fluidH(i0).state(2,1).h)/1e6, (fluidH(i0).state(2,2).s - fluidH(i0).state(2,1).s)/1e3, fluidH(i0).state(2,1).mdot, i0, 'dis');end
    fprintf(1,'\nCold fluid streams:\n');
    fprintf(1,'-->%s\n',fluidC(1).name);
    fprintf(1,'%11s ','Tin[K]','Tout[K]','Δh [MJ/kg]','Δs [kJ/kg/K]','mdot[kg/s]','Stream','Cycle'); fprintf(1,'\n');
    for i0=1:n3, fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.2f %11d %11s\n', fluidC(i0).state(1,1).T, fluidC(i0).state(1,2).T, (fluidC(i0).state(1,2).h - fluidC(i0).state(1,1).h)/1e6, (fluidC(i0).state(1,2).s - fluidC(i0).state(1,1).s)/1e3, fluidC(i0).state(1,1).mdot, i0, 'ch');end
    for i0=1:n3, fprintf(1,'%11.1f %11.1f %11.2f %11.2f %11.2f %11d %11s\n', fluidC(i0).state(2,1).T, fluidC(i0).state(2,2).T, (fluidC(i0).state(2,2).h - fluidC(i0).state(2,1).h)/1e6, (fluidC(i0).state(2,2).s - fluidC(i0).state(2,1).s)/1e3, fluidC(i0).state(2,1).mdot, i0, 'dis');end
    fprintf(1,'\nHot tanks\n');
    fprintf(1,'%10s %10s %13s %13s %13s %13s %8s ','A.T [K]','A.M [kg*1e6]','A.B [MWh]','B.T [K]','B.M [kg*1e6]','B.B [MWh]','state'); fprintf(1,'\n');
    for i0=1:4
        fprintf(1,'%10.4g %13.3f %13.3f %10.4g %13.3f %13.3f %8d\n', HT.A(i0).T,HT.A(i0).M/1e6,HT.A(i0).B/3600/1e6,HT.B(i0).T,HT.B(i0).M/1e6,HT.B(i0).B/3600/1e6,i0)
    end    
    fprintf(1,'\nCold tanks\n');
    fprintf(1,'%10s %10s %13s %13s %13s %13s %8s ','A.T [K]','A.M [kg*1e6]','A.B [MWh]','B.T [K]','B.M [kg*1e6]','B.B [MWh]','state'); fprintf(1,'\n');
    for i0=1:4
        fprintf(1,'%10.4g %13.3f %13.3f %10.4g %13.3f %13.3f %8d\n', CT.A(i0).T,CT.A(i0).M/1e6,CT.A(i0).B/3600/1e6,CT.B(i0).T,CT.B(i0).M/1e6,CT.B(i0).B/3600/1e6,i0)
    end    
    fprintf(1,'\n');
end

fact = 1e6*3600; % J to MWh

Net_ch = (W_in_ch  + QC_ch  - QH_ch  - DH_ch  - QExt_ch);
Net_dis = (QH_dis  - W_out_dis  - QC_dis   - DH_dis  - QExt_dis);

WL_matrix = [ WL_PTES_ch ; WL_PTES_dis ; ];
Total_loss = sum(WL_matrix(:));

WL_turbo   = sum((WL_matrix(:,1) + WL_matrix(:,2)))/W_in_ch*100;
WL_hexs    = sum(WL_matrix(:,3))/W_in_ch*100;
WL_reject  = sum(WL_matrix(:,4))/W_in_ch*100;
WL_mix_liq = sum(WL_matrix(:,5))/W_in_ch*100;
WL_tanks   = sum(WL_matrix(:,6))/W_in_ch*100;

[WR,WR_dis]   = work_ratio(gas.stage,stages_ch,stages_dis);

% Compute efficiencies, energy and power densities and errors
switch mode
    case 0 % PTES
        Heat_in_tanks   = ((HT.A(4).H - HT.A(1).H) + (HT.B(4).H - HT.B(1).H) + (CT.A(4).H - CT.A(1).H) + (CT.B(4).H - CT.B(1).H));
        Heat_rejected = QExt_ch + QExt_dis;
        Total_Work_lost = W_in_ch - W_out_dis;
        First_law_error = (Heat_rejected + Heat_in_tanks - Total_Work_lost)/Total_Work_lost;
        Second_law_error = (Total_loss - Total_Work_lost)/Total_Work_lost;
        chi_PTES = W_out_dis/W_in_ch;
        rhoE = W_out_dis/fact/(HT.A(1).V + CT.A(1).V)*1e3; %kWh/m3
        rhoP_ch  = (W_in_ch/t_ch/(gas_min_rho_ch.mdot/gas_min_rho_ch.rho)/1e6); %MW/(m3/s)
        rhoP_dis = (W_out_dis/t_dis/(gas_min_rho_dis.mdot/gas_min_rho_dis.rho)/1e6); %MW/(m3/s)
        
    case 1 % Heat pump only
        Heat_into_tanks   = ((HT.A(2).H - HT.A(1).H) + (HT.B(2).H - HT.B(1).H) + (CT.A(2).H - CT.A(1).H) + (CT.B(2).H - CT.B(1).H));
        Exergy_into_tanks = ((HT.A(2).B - HT.A(1).B) + (HT.B(2).B - HT.B(1).B) + (CT.A(2).B - CT.A(1).B) + (CT.B(2).B - CT.B(1).B));
        Heat_rejected = QExt_ch;
        Total_Work_lost = W_in_ch - Exergy_into_tanks;
        First_law_error = (Heat_rejected + Heat_into_tanks - W_in_ch)/(W_in_ch);
        Second_law_error = (Total_loss - Total_Work_lost)/Total_Work_lost;
        chi_tot = Exergy_into_tanks/(W_in_ch);
        Exergy_into_hot_tanks = ((HT.A(2).B - HT.A(1).B) + (HT.B(2).B - HT.B(1).B));
        chi_hot = Exergy_into_hot_tanks/(W_in_ch);
        COP = QH_ch/W_in_ch;
        rhoE = Exergy_into_hot_tanks/fact/(HT.A(1).V)*1e3; %kWh/m3
        rhoP_ch  = (W_in_ch/t_ch/(gas_min_rho_ch.mdot/gas_min_rho_ch.rho)/1e6); %MW/(m3/s)
        
    case 2 % Heat engine only
        Heat_from_tanks   = ((HT.A(3).H - HT.A(4).H) + (HT.B(3).H - HT.B(4).H));
        Exergy_from_tanks = ((HT.A(3).B - HT.A(4).B) + (HT.B(3).B - HT.B(4).B));
        Heat_rejected = QExt_dis;
        Total_Work_lost = Exergy_from_tanks - W_out_dis;
        First_law_error = (Heat_from_tanks - W_out_dis - Heat_rejected)/(W_out_dis);
        Second_law_error = (Total_loss - Total_Work_lost)/Total_Work_lost;
        chi_tot = (W_out_dis)/Exergy_from_tanks;
        EFF = W_out_dis/QH_dis;
        rhoE = Exergy_from_tanks/fact/(HT.A(4).V)*1e3; %kWh/m3
        rhoP_dis = (W_out_dis/t_dis/(gas_min_rho_dis.mdot/gas_min_rho_dis.rho)/1e6); %MW/(m3/s)
end


if WM==1
    fprintf(1,'\n\n');
    fprintf(1,'FIRST LAW BALANCE\n');
    fprintf(1,'-----------------\n');
    
    switch mode
        case {0,1}
            fprintf(1,'CHARGE\n');
            fprintf(1,'Net power input:         %8.1f MW\n',W_in_ch/t_ch/1e6);
            fprintf(1,'Charge time:             %8.1f h\n',t_ch/3600);
            fprintf(1,'Energy(el) input:        %8.1f MWh\n',W_in_ch/fact);
            fprintf(1,'Heat to hot tanks:       %8.1f MWh\n',QH_ch/fact);
            fprintf(1,'Heat from cold tanks:    %8.1f MWh\n',QC_ch/fact);
            fprintf(1,'DH working fluid:        %8.1f MWh\n',DH_ch/fact);
            fprintf(1,'Heat rejected:           %8.1f MWh\n',QExt_ch/fact);
            fprintf(1,'NET:                     %8.1f MWh\n\n',Net_ch/fact);
    end
    
    switch mode
        case {0,2}
            fprintf(1,'DISCHARGE\n');
            fprintf(1,'Net power output:        %8.1f MW\n',W_out_dis/t_dis/1e6);
            fprintf(1,'Discharge time:          %8.1f h\n',t_dis/3600);
            fprintf(1,'Energy(el) output:       %8.1f MWh\n',W_out_dis/fact);
            fprintf(1,'Heat from hot tanks:     %8.1f MWh\n',QH_dis/fact);
            fprintf(1,'Heat to cold tanks:      %8.1f MWh\n',QC_dis/fact);
            fprintf(1,'DH working fluid:        %8.1f MWh\n',DH_dis/fact);
            fprintf(1,'Heat rejected:           %8.1f MWh\n',QExt_dis/fact);
            fprintf(1,'NET:                     %8.1f MWh\n\n',Net_dis/fact);
    end
    
    switch mode
        case 0 
            fprintf(1,'Round trip efficiency:   %8.1f %%\n\n',chi_PTES*100);
            fprintf(1,'Exergy density:          %9.2f kWh/m3\n',rhoE);
            fprintf(1,'Power density (charge):  %9.2f MW/(m3/s)\n',rhoP_ch);
            fprintf(1,'Power density (disch):   %9.2f MW/(m3/s)\n\n',rhoP_dis);
            
            fprintf(1,'STORAGE MEDIA\n');
            fprintf(1,'%18s volume:%8.2f m3/MWh\n',fluidH(1).name,HT.A(1).V/(W_out_dis/fact));
            fprintf(1,'%18s volume:%8.2f m3/MWh\n\n',fluidC(1).name,CT.A(1).V/(W_out_dis/fact));
            
            fprintf(1,'%18s mass:  %8.2f tons/MWh\n',fluidH(1).name,HT.A(1).M/(W_out_dis/fact)/1e3);
            fprintf(1,'%18s mass:  %8.2f tons/MWh\n',fluidC(1).name,CT.A(1).M/(W_out_dis/fact)/1e3);
            
            fprintf(1,'Work Ratios:  %6.3g (charge) %6.3g (discharge)\n\n',WR,WR_dis);
        case 1
            fprintf(1,'Exergetic efficiency:                %8.1f %%\n',chi_tot*100);
            fprintf(1,'Exergetic efficiency (cold reject):  %8.1f %%\n',chi_hot*100);
            fprintf(1,'Coefficient of Performance:           %8.2f\n\n',COP);
            fprintf(1,'Exergy density (hot tanks):          %9.2f kWh/m3\n',rhoE);
            fprintf(1,'Power density (charge):              %9.2f MW/(m3/s)\n\n',rhoP_ch);
        case 2
            fprintf(1,'Exergetic efficiency:      %7.1f %%\n',chi_tot*100);
            fprintf(1,'First law efficiency:      %7.1f %%\n\n',EFF*100);
            fprintf(1,'Exergy density:            %8.2f kWh/m3\n',rhoE);
            fprintf(1,'Power density (discharge): %8.2f MW/(m3/s)\n\n',rhoP_dis);
    end
    
    fprintf(1,'First law error:   %8.5f %%\n',First_law_error*100);
    fprintf(1,'Second law error:  %8.5f %%\n\n',Second_law_error*100);
end


% Manage warnings
if  ((CT.A(1).T > T0) && (CT.A(4).T < (CT.A(1).T - 1))) && mode == 0
    warning('Unsustainable discharge of cold reservoir!')
end
