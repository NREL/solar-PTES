% First and second law balances

% Calculate work/heat/irreversibility terms for each component
% Charging compressors and discharging expanders
for ii = 1:length(CCMP)
    CCMP(ii) = compexp_energy(CCMP(ii),Load.time)  ; %#ok<*SAGROW>
end
% Charging expanders and discharging compressors
for ii = 1:length(CEXP)
    CEXP(ii) = compexp_energy(CEXP(ii),Load.time)  ;
end

for ii = 1:length(DCMP)
    DCMP(ii) = compexp_energy(DCMP(ii),Load.time)  ;
end
for ii = 1:length(DEXP)
    DEXP(ii) = compexp_energy(DEXP(ii),Load.time)  ;
end

% Recompressor if specified for sCO2 cycle
if any(Load.mode == [4,5,6])
    if Lrcmp
        RCMP = compexp_energy(RCMP,Load.time)  ;
    end
end

% Total energy flows for heat exchangers
for ii = 1:length(HX)
    HX(ii) = hx_energy(HX(ii),Load.time)  ;
end


% PARASITIC LOSSES
% FANS
for ii = 1:length(CFAN)
    CFAN(ii) = compexp_energy(CFAN(ii),Load.time)  ;
end
for ii = 1:length(DFAN)
    DFAN(ii) = compexp_energy(DFAN(ii),Load.time)  ;
end

% FLUID PUMPS
for ii = 1:length(CPMP)
    CPMP(ii) = compexp_energy(CPMP(ii),Load.time)  ;
end
for ii = 1:length(DPMP)
    DPMP(ii) = compexp_energy(DPMP(ii),Load.time)  ;
end

% Require similar routines for mixers etc.
% MIXERS
for ii = 1:length(MIX)
    MIX(ii) = misc_energy(MIX(ii),Load.time) ;
end

% Motor-Generator
switch Load.mode
    case {2,7}
        error('not implemented yet')
    case {3}
        GEN(1) = gen_power(GEN(1),CCMP,CEXP,DCMP,DEXP,Load.time);
        GEN(2) = gen_power(GEN(2),CCMP,CEXP,DCMP,DEXP,Load.time);
    case {4,5}
        GEN = gen_power(GEN,CCMP,CEXP,[DCMP RCMP],DEXP,Load.time);
    case {6}
        GEN(1) = gen_power(GEN(1),CCMP,CEXP,[DCMP RCMP],DEXP,Load.time);
        GEN(2) = gen_power(GEN(2),CCMP,CEXP,[DCMP RCMP],DEXP,Load.time);
    otherwise
        GEN = gen_power(GEN,CCMP,CEXP,DCMP,DEXP,Load.time);
end

switch Load.mode
    case {0,1,2,4,5,6}
        NC = Ne_ch ;
        NE = Nc_ch ;
    case {3,7}
        NC = 3 ;
        NE = 3 ;
end

% Adds up the contributions of the several stages to compute an energy balance
% SIGN CONVENTION. Work out is positive. Heat in is positive.
E_net_chg  = 0; % Electricity input to plant
E_in_chg   = 0; % Electricity input to motor (after fan and pump parasitics)
W_in_chg   = 0; % Mechanical work input to machines (after parasitics and motor losses)
W_CCMP     = 0; % Mechanical work input to compressors (after parasitics and motor losses)
W_CEXP     = 0; % Mechanical work output from turinbes (after parasitics and motor losses)
WL_mot_chg = 0; % Lost work by electrical motor
W_pmp_chg  = 0; % Parasitic work input - nice to separate this out
W_fan_chg  = 0; % Parasitic work input - nice to separate this out
DH_chg     = 0;
W_lost_chg = 0;
QH_chg     = 0;  % heat to hot tanks
QC_chg     = 0;  % heat from cold tanks
QE_chg     = 0;  % heat rejected to environment
QT_chg     = 0;  % heat to/from any heat exchanger

E_net_dis  = 0;  % Electricity output from plant (after generator and parasitic losses)
E_out_dis  = 0;  % Electricity output from generator (after generator losses)
W_out_dis  = 0;  % Mechanical work output from machines
W_DCMP     = 0; % Mechanical work input to compressors (after parasitics and motor losses)
W_DEXP     = 0; % Mechanical work output from turinbes (after parasitics and motor losses)
W_EXP_dis  = 0;  % Mechanical work output from turbines
WL_gen_dis = 0;  % Work lost by electrical generator
W_pmp_dis  = 0;  % Parasitic work input - nice to separate this out
W_fan_dis  = 0;  % Parasitic work input - nice to separate this out
DH_dis     = 0;
W_lost_dis = 0;
QH_dis     = 0;  % heat from hot tanks
QC_dis     = 0;  % heat to cold tanks
QE_dis     = 0;  % heat rejected to environment
QT_dis     = 0;  % heat to/from any heat exchanger

E_out_disRC = 0; % Electricity output (after generator losses) for cycle using cold stores
W_out_disRC = 0;
W_fan_disRC = 0;  % Parasitic work input - nice to separate this out
W_pmp_disRC = 0;  % Parasitic work input - nice to separate this out
QC_disRC    = 0;  
QH_disRC    = 0; 
t_disRC     = 0;

nH = numel(fluidH);
nM = numel(fluidM);
nC = numel(fluidC);

% Compute lost work on specific component types
% The WL_PTES_chg and WL_PTES_dis arrays are divided in 7 elements:
% 1: Losses in compressors
% 2: Losses in expanders
% 3: Losses in heat exchangers
% 4: Losses due to heat exchange with the environment
% 5: Mixing losses in storage tanks (liquids)
% 6: Mixing losses of the working fluid
% 7: Losses due to exergy leftover in tanks
% 8: Parasitic losses due to pumping fluid and heat rejection air.
WL_PTES_chg = zeros(1,9);
WL_PTES_dis = zeros(1,9);

for iL=1:Load.num
    
    % For charging load cycles, account for work, heat and irreversibilities
    if any(strcmp(Load.type(iL),{'chg','chgCO2','chgTSCO2','chgCC','chgICC','chgICC_PC'}))
        
        % Work into cycle
        for ii = 1:length(CCMP)
            W_in_chg       = W_in_chg   + CCMP(ii).W(iL) ; 
            W_CCMP         = W_CCMP     + CCMP(ii).W(iL) ; 
            WL_PTES_chg(1) = WL_PTES_chg(1) + T0 * CCMP(ii).Sirr(iL) ; 
        end
        
        % Work out of cycle
        for ii = 1:length(CEXP)
            W_in_chg   = W_in_chg   + CEXP(ii).W(iL) ; 
            W_CEXP     = W_CEXP     + CEXP(ii).W(iL) ; 
            WL_PTES_chg(2) = WL_PTES_chg(2) + T0 * CEXP(ii).Sirr(iL) ;
        end
        
        % Parasitic work into cycle
        for ii = 1:length(CFAN)
            W_fan_chg      = W_fan_chg + CFAN(ii).W(iL) ; 
            WL_PTES_chg(8) = WL_PTES_chg(8) + T0 * CFAN(ii).Sirr(iL) ;
        end
        for ii = 1:length(CPMP)
            W_pmp_chg      = W_pmp_chg + CPMP(ii).W(iL) ; 
            WL_PTES_chg(8) = WL_PTES_chg(8) + T0 * CPMP(ii).Sirr(iL) ;
        end
        
        % Heat flows in and out of cycle
        for ii = 1:length(HX)
            % Heat out of cycle into hot storage
            if strcmp(HX(ii).name,'hot')
                QH_chg = QH_chg + HX(ii).Q(iL,1) ;
            end
            
            % Heat into cycle from cold storage
            if strcmp(HX(ii).name,'cold')
                QC_chg = QC_chg + HX(ii).Q(iL,2) ;
            end
            
            % Heat out of cycle to the environment
            if strcmp(HX(ii).name,'rej')
                QE_chg = QE_chg + HX(ii).Q(iL,1);
                WL_PTES_chg(4) = WL_PTES_chg(4) + T0 * HX(ii).Sirr(iL,2) ; % Loss due to heat exchange with environment
            elseif strcmp(HX(ii).name,'hin')
                QE_chg = QE_chg + HX(ii).Q(iL,2) ;
                WL_PTES_chg(4) = WL_PTES_chg(4) + T0 * HX(ii).Sirr(iL,2) ; % Loss due to heat exchange with environment
            else
                WL_PTES_chg(3) = WL_PTES_chg(3) + T0 * HX(ii).Sirr(iL,2) ; % Loss due to heat exchange
            end
            
            % Heat exchanged on any heat exchanger
            QT_chg = QT_chg + abs(HX(ii).Q(iL,1));
        end
        
        % Electricity from Motor
        for ii = 1 : numel(GEN)
            E_in_chg   = E_in_chg   + GEN(ii).E(iL);
            WL_mot_chg = WL_mot_chg - GEN(ii).WL(iL);
        end
        
        % Work lost in other components (mixers, seperators, etc.)
        for ii = 1:length(MIX)
            WL_PTES_chg(6) = WL_PTES_chg(6) + T0 * MIX(ii).Sirr(iL,2) ;
        end
        
        
    % For discharging load cycles, account for work, heat and irreversibilities
    elseif any(strcmp(Load.type(iL),{'dis','ran','disCO2','rcmpCO2','disTSCO2','disCC','disICC','disICC_PC'}))
        % Work into cycle
        for ii = 1:length(DCMP)
            W_out_dis      = W_out_dis  + DCMP(ii).W(iL) ; 
            W_DCMP         = W_DCMP     + DCMP(ii).W(iL) ; 
            WL_PTES_dis(1) = WL_PTES_dis(1) + T0 * DCMP(ii).Sirr(iL) ; 
        end
        
        % Work out of cycle
        for ii = 1:length(DEXP)
            W_out_dis      = W_out_dis  + DEXP(ii).W(iL) ; 
            W_DEXP         = W_DEXP     + DEXP(ii).W(iL) ; 
            WL_PTES_dis(2) = WL_PTES_dis(2) + T0 * DEXP(ii).Sirr(iL) ; 
        end
        if any(Load.mode == [4,5,6])
            if Lrcmp
                W_out_dis      = W_out_dis  + RCMP.W(iL) ;
                WL_PTES_dis(1) = WL_PTES_dis(1) + T0 * RCMP.Sirr(iL) ; 
            end
        end
        
        % Parasitic work into cycle
        for ii = 1:length(DFAN)
            W_fan_dis      = W_fan_dis + DFAN(ii).W(iL) ; 
            WL_PTES_dis(8) = WL_PTES_dis(8) + T0 * DFAN(ii).Sirr(iL) ; 
        end
        for ii = 1:length(DPMP)
            W_pmp_dis      = W_pmp_dis + DPMP(ii).W(iL) ; 
            WL_PTES_dis(8) = WL_PTES_dis(8) + T0 * DPMP(ii).Sirr(iL) ; 
        end
        
        % Heat flows in and out of cycle
        for ii = 1:length(HX)
            % Heat into cycle from hot storage
            if strcmp(HX(ii).name,'hot')
                QH_dis = QH_dis + HX(ii).Q(iL,2) ;
            end
            
            % Heat out of cycle to cold storage
            if strcmp(HX(ii).name,'cold')
                QC_dis = QC_dis + HX(ii).Q(iL,1) ;
            end
            
            % Heat out of cycle to the environment
            if strcmp(HX(ii).name,'rej')
                QE_dis = QE_dis + HX(ii).Q(iL,1) ;
                WL_PTES_dis(4) = WL_PTES_dis(4) + T0 * HX(ii).Sirr(iL,2) ; % Loss due to heat exchange with environment
            else
                WL_PTES_dis(3) = WL_PTES_dis(3) + T0 * HX(ii).Sirr(iL,2) ; % Loss due to heat exchange
            end
            
        end
        
        % Electricity from Generator
        for ii = 1 : numel(GEN)
            E_out_dis  = E_out_dis  + GEN(ii).E(iL);
            WL_gen_dis = WL_gen_dis - GEN(ii).WL(iL);
        end
        
        % Work lost in other components
        for ii = 1:length(MIX)
            WL_PTES_dis(6) = WL_PTES_dis(6) + T0 * MIX(ii).Sirr(iL,2) ;
        end
        
        % Calculate heat and work terms separately for a Rankine cycle that uses cooling
        if strcmp(Load.type(iL),'ran') && Load.options.useCold(iL) == 1
                       
            % Work into cycle
            for ii = 1:length(DCMP)
                W_out_disRC = W_out_disRC + DCMP(ii).W(iL) ;
            end
            
            % Work out of cycle
            for ii = 1:length(DEXP)
                W_out_disRC = W_out_disRC + DEXP(ii).W(iL) ;
            end
            
            % Electricity from Generator
            for ii = 1 : numel(GEN)
                E_out_disRC  = E_out_disRC  + GEN(ii).E(iL);
            end
            
            % Parasitic work into cycle
            for ii = 1:length(DFAN)
                W_fan_disRC = W_fan_disRC + DFAN(ii).W(iL) ;
            end
            for ii = 1:length(DPMP)
                W_pmp_disRC = W_pmp_disRC + DPMP(ii).W(iL) ;
            end
            
            % Heat flows in and out of cycle
            for ii = 1:length(HX)
                % Heat into cycle from hot storage
                if strcmp(HX(ii).name,'hot')
                    QH_disRC = QH_disRC + HX(ii).Q(iL,2) ;
                end
                
                % Heat out of cycle to cold storage
                if strcmp(HX(ii).name,'cold')
                    QC_disRC = QC_disRC + HX(ii).Q(iL,1) ;
                end
            end
            
            t_disRC = t_disRC + Load.time(iL);
        end
    
        
    end
    
end

if Load.mode == 3
    % Calculate heat and work terms separately for a Rankine cycle that does not use cooling
    W_out_disNC = W_out_dis - W_out_disRC;
    E_out_disNC = E_out_dis - E_out_disRC;
    W_fan_disNC = W_fan_dis - W_fan_disRC;
    W_pmp_disNC = W_pmp_dis - W_pmp_disRC;
    QH_disNC    = QH_dis    - QH_disRC;
    QC_disNC    = QC_dis    - QC_disRC;
end

switch Load.mode
    case {20,22,24}
        % Calculate enthalpy flow from working fluid (air liquefaction)
        QWF_chg = -(LAT.DH_chg + AT.DH_chg)-W_fan_chg;
        QWF_dis = -(LAT.DH_dis + AT.DH_dis)-W_fan_dis;
        QE_chg  = 0;
        QE_dis  = 0;
    otherwise
        QWF_chg = 0;
        QWF_dis = 0;
end

% Add in losses from tanks (mixing and exergy left inside)
switch Load.mode

    case {0,3,4,6,20,22,24} % PTES, Rankine, or sCO2-PTES
        
        for ii = 1 : Nhot
            WL_PTES_chg(5) = WL_PTES_chg(5) + HT(ii).WL_chg ;
            WL_PTES_dis(5) = WL_PTES_dis(5) + HT(ii).WL_dis ;
            WL_PTES_dis(7) = WL_PTES_dis(7) + HT(ii).A(end).B - HT(ii).A(1).B + HT(ii).B(end).B - HT(ii).B(1).B;
        end
        
        for ii = 1 : Nmid
            WL_PTES_chg(5) = WL_PTES_chg(5) + MT(ii).WL_chg ;
            WL_PTES_dis(5) = WL_PTES_dis(5) + MT(ii).WL_dis ;
            WL_PTES_dis(7) = WL_PTES_dis(7) + MT(ii).A(end).B - MT(ii).A(1).B + MT(ii).B(end).B - MT(ii).B(1).B;
        end
        
        for ii = 1 : Ncld
            WL_PTES_chg(5) = WL_PTES_chg(5) + CT(ii).WL_chg ;
            WL_PTES_dis(5) = WL_PTES_dis(5) + CT(ii).WL_dis ;
            WL_PTES_dis(7) = WL_PTES_dis(7) + CT(ii).A(end).B - CT(ii).A(1).B + CT(ii).B(end).B - CT(ii).B(1).B;
        end
        
        for ii = 1 : Nlat
            WL_PTES_chg(5) = WL_PTES_chg(5) + LAT(ii).WL_chg ;
            WL_PTES_dis(5) = WL_PTES_dis(5) + LAT(ii).WL_dis ;
            WL_PTES_dis(7) = WL_PTES_dis(7) + LAT(ii).A(end).B - LAT(ii).A(1).B;
        end
        
        WL_PTES_chg(4) = WL_PTES_chg(4) + AT.WL_chg ;
        WL_PTES_dis(4) = WL_PTES_dis(4) + AT.WL_dis ;

    case 1 % Heat pump only
        
        for ii = 1 : Nhot
            WL_PTES_chg(5) = WL_PTES_chg(5) + HT(ii).WL_chg ;
        end
           
        for ii = 1 : Ncld
            WL_PTES_chg(5) = WL_PTES_chg(5) + CT(ii).WL_chg ;
        end
        
        WL_PTES_chg(5) = WL_PTES_chg(5) + AT.WL_chg ;

    case {2,5,7} % Heat engine only
        WL_PTES_dis(5) = HT.WL_dis + AT.WL_dis ;
        WL_PTES_dis(7) = HT.A(end).B - HT.A(1).B + HT.B(end).B - HT.B(1).B;
end

% Allocate losses from motor/generator
WL_PTES_chg(9) = abs(WL_mot_chg);
WL_PTES_dis(9) = abs(WL_gen_dis);

fact    = 1e6*3600; % J to MWh

i_chg = Load.ind(any(Load.type == {'chg','chgCO2','chgTSCO2','chgCC','chgICC','chgICC_PC'},2));
i_dis = Load.ind(any(Load.type == {'dis','disCO2','ran','rcmpCO2','disTSCO2','disCC','disICC','disICC_PC'},2));
i_ran = Load.ind(Load.type == 'ran');
i_gas = Load.ind(any(Load.type == {'chg','dis','chgCO2','disCO2','rcmpCO2','chgTSCO2','disTSCO2','chgCC','disCC','chgICC','disICC','chgICC_PC','disICC_PC'},2));
i_act = Load.ind(any(Load.type == {'chg','dis','chgCO2','disCO2','ran','rcmpCO2','chgTSCO2','disTSCO2','chgCC','disCC','chgICC','disICC','chgICC_PC','disICC_PC'},2));

t_chg = sum(Load.time(i_chg));
t_dis = sum(Load.time(i_dis));
ip1 = find(i_chg == 1,1,'first'); %index for printing (first charge period)
ip2 = find(i_dis == 1,1,'first'); %index for printing (first discharge period)
t_disNC = t_dis - t_disRC;

% FIRST LAW ENERGY BALANCES AND ACCOUNTING
if Load.mode==24
    QC_chg = - QC_chg;
    QC_dis = - QC_dis;
end
Net_chg = (- W_in_chg  + QC_chg  + QH_chg  + QE_chg + QWF_chg);
Net_dis = (- W_out_dis + QC_dis  + QH_dis  + QE_dis + QWF_dis);

% SECOND LAW ENERGY BALANCES AND ACCOUNTING
WL_matrix  = [ WL_PTES_chg ; WL_PTES_dis ; ];
Total_loss = sum(WL_matrix(:));

% Calculate tank stats - volumes and masses
for ii = 1 : Nhot
    HT(ii) = tank_stats(HT(ii)) ;
end
for ii = 1 : Nmid
    MT(ii) = tank_stats(MT(ii)) ;
end
for ii = 1 : Ncld
    CT(ii) = tank_stats(CT(ii)) ;
end
for ii = 1 : Nlat
    LAT(ii) = single_tank_stats(LAT(ii));
    LAT(ii).tank_volA = LAT(ii).fluid_volA;
end

% Compute efficiencies, energy and power densities and errors
switch Load.mode
    case {0,3,4,6,20,22,24} % PTES, Rankine, sCO2-PTES, or PTES-LAES CC
        Heat_in_tanks = 0;
        for ii = 1 : Nhot
            Heat_in_tanks = Heat_in_tanks + (HT(ii).A(end).H - HT(ii).A(1).H) + (HT(ii).B(end).H - HT(ii).B(1).H) ;
        end
        for ii = 1 : Nmid
            Heat_in_tanks = Heat_in_tanks + (MT(ii).A(end).H - MT(ii).A(1).H) + (MT(ii).B(end).H - MT(ii).B(1).H) ;
        end
        for ii = 1 : Ncld
            Heat_in_tanks = Heat_in_tanks + (CT(ii).A(end).H - CT(ii).A(1).H) + (CT(ii).B(end).H - CT(ii).B(1).H);
        end
        for ii = 1 : Nlat
            Heat_in_tanks = Heat_in_tanks + (LAT(ii).A(end).H -LAT(ii).A(1).H);
        end
        
        Heat_rejected    = QWF_chg + QWF_dis + QE_chg + QE_dis + W_fan_chg + W_fan_dis + WL_mot_chg + WL_gen_dis;
        Total_Work_lost  = E_in_chg + E_out_dis + W_fan_chg + W_pmp_chg + W_fan_dis + W_pmp_dis;
        First_law_error  = (-Heat_rejected + Heat_in_tanks + Total_Work_lost)/Total_Work_lost;
        Second_law_error = (Total_loss + Total_Work_lost)/Total_Work_lost;
        
        E_net_chg = E_in_chg  + W_fan_chg + W_pmp_chg;
        E_net_dis = E_out_dis + W_fan_dis + W_pmp_dis;
        
        chi_PTES = -W_out_dis/W_in_chg;
        chi_PTES_para = -E_net_dis / E_net_chg ;
        
        % Calculate storage volume
        volH = 0.0; % Calculate storage volume
        volM = 0.0; % Calculate storage volume
        volC = 0.0; % Calculate storage volume
        volL = 0.0; % Calculate storage volume
        for ii = 1:Nhot; volH = volH + HT(ii).tank_volA + HT(ii).tank_volB ; end
        for ii = 1:Nmid; volM = volM + MT(ii).tank_volA + MT(ii).tank_volB ; end
        for ii = 1:Ncld; volC = volC + CT(ii).tank_volA + CT(ii).tank_volB ; end
        for ii = 1:Nlat; volL = volL + LAT(ii).tank_volA ; end
        vol  = volH + volM + volC + volL;
        
        % Energy density
        rhoE = E_net_dis /fact/vol*1e3; %kWh/m3
        
        % Calculate some special metrics for certain cycles
        if Load.mode == 3
            if i_chg == 1 && i_dis(1) == 3 && i_dis(2) == 4
                Exergy_into_hot_tanks    = (HT.A(2).B - HT.A(1).B) + (HT.B(2).B - HT.B(1).B);
                Exergy_into_cold_tanks   = (CT.A(2).B - CT.A(1).B) + (CT.B(2).B - CT.B(1).B);
                Exergy_from_hot_tanks    = (HT.A(3).B - HT.A(5).B) + (HT.B(3).B - HT.B(5).B);
                Exergy_from_hot_tanksRC  = (HT.A(3).B - HT.A(4).B) + (HT.B(3).B - HT.B(4).B);
                Exergy_from_hot_tanksNC  = (HT.A(4).B - HT.A(5).B) + (HT.B(4).B - HT.B(5).B);
                Exergy_from_cold_tanks   = (CT.A(3).B - CT.A(5).B) + (CT.B(3).B - CT.B(5).B);
                Exergy_from_cold_tanksRC = (CT.A(3).B - CT.A(4).B) + (CT.B(3).B - CT.B(4).B);
                Exergy_from_cold_tanksNC = (CT.A(4).B - CT.A(5).B) + (CT.B(4).B - CT.B(5).B);
                Exergy_into_tanks   = Exergy_into_hot_tanks   + Exergy_into_cold_tanks;
                Exergy_from_tanks   = Exergy_from_hot_tanks   + Exergy_from_cold_tanks;
                Exergy_from_tanksRC = Exergy_from_hot_tanksRC + Exergy_from_cold_tanksRC;
                Exergy_from_tanksNC = Exergy_from_hot_tanksNC + Exergy_from_cold_tanksNC;
            else
                warning('Exergy efficiencies not computed')
                Exergy_into_hot_tanks  = NaN;
                Exergy_into_cold_tanks = NaN;
                Exergy_from_hot_tanks  = NaN;
                Exergy_from_cold_tanks = NaN;
            end
            
            if ~any(Load.options.useCold)
                vol  = volH;
                rhoE = E_net_dis /fact/vol*1e3; %kWh/m3
            end
            
            E_net_disRC = E_out_disRC + W_fan_disRC + W_pmp_disRC;
            E_net_disNC = E_out_disNC + W_fan_disNC + W_pmp_disNC;
            
            % COP and heat engine efficiency (including parasitics)
            COP     = QH_chg / E_net_chg;
            HEeff   = E_net_dis   / QH_dis ;   % Heat engine average efficiency
            HEeffRC = E_net_disRC / QH_disRC ; % Efficiency for cycles where only cold tanks are used for condensing
            HEeffNC = E_net_disNC / QH_disNC ; % Efficiency for cycles where cold tanks are NOT used for condensing
            
            % Heat Pump and Heat Engine exergetic efficiencies (including
            % parasitics)
            HPexergy_eff   = Exergy_into_tanks / (-E_net_chg); % Heat pump average exergetic efficiency
            HEexergy_eff   = E_net_dis / Exergy_from_tanks; % Heat engine average exergetic efficiency
            HEexergy_effRC = E_net_disRC / Exergy_from_tanksRC;
            HEexergy_effNC = E_net_disNC / Exergy_from_tanksNC;

        elseif Load.mode == 6
            QH_sol = 0;
            EX_sol = 0;
            for ii=1:Load.num
               if strcmp(Load.type(ii),'disTSCO2')
                  QH_sol = QH_sol + fluidH(1).state(ii,1).mdot * (fluidH(1).state(ii,1).h-fluidH(1).state(ii,2).h) * Load.time(ii) ;
                  EX_sol = EX_sol + fluidH(1).state(ii,1).mdot * (fluidH(1).state(ii,1).h-fluidH(1).state(ii,2).h-T0*(fluidH(1).state(ii,1).s-fluidH(1).state(ii,2).s)) * Load.time(ii) ;
               end
            end
            QH_chg  = (gas.state(1,4).h-gas.state(1,2).h)*gas.state(1,3).mdot*Load.time(1) ;
            SOLeff  = E_net_dis / QH_sol ; % Solar conversion efficiency
            HEeff   = E_net_dis / QH_dis ; % Heat engine efficiency
            NETeff  = (E_net_dis + E_net_chg) / QH_sol ; % Net efficiency
            EXeff   = E_net_dis / (-E_net_chg + EX_sol) ; % Exergetic efficiency
            
        elseif any(Load.mode == [22,24])
            W_ratio    = W_DEXP/(-W_DCMP);
            Qstr_W_ratio = (abs(QH_chg)+abs(QC_chg))/E_net_dis;
            QT_W_ratio = QT_chg/E_net_dis;
        end
        
    case 1 % Heat pump only
        Heat_into_tanks   = (HT.A(end).H - HT.A(1).H) + (HT.B(end).H - HT.B(1).H) + (CT.A(end).H - CT.A(1).H) + (CT.B(end).H - CT.B(1).H);
        Exergy_into_tanks = (HT.A(end).B - HT.A(1).B) + (HT.B(end).B - HT.B(1).B) + (CT.A(end).B - CT.A(1).B) + (CT.B(end).B - CT.B(1).B);
        
        Heat_rejected = QE_chg + W_fan_chg + WL_mot_chg;
        E_net_chg = E_in_chg + W_fan_chg + W_pmp_chg + WL_mot_chg;
        Total_Work_lost = -E_net_chg - Exergy_into_tanks;
        
        First_law_error = (-Heat_rejected + Heat_into_tanks + E_net_chg)/E_net_chg;
        Second_law_error = (Total_loss - Total_Work_lost)/Total_Work_lost;
                
        chi_tot               = Exergy_into_tanks/(W_in_chg);
        chi_tot_para          = Exergy_into_tanks/(E_net_chg);
        Exergy_into_hot_tanks = ((HT.A(end).B - HT.A(1).B) + (HT.B(end).B - HT.B(1).B));
        chi_hot               = Exergy_into_hot_tanks/(W_in_chg);
        chi_hot_para          = Exergy_into_hot_tanks/(E_net_chg);
        COP                   = QH_chg/W_in_chg;
        COP_para              = QH_chg/E_net_chg;
        vol = 0.0; % Calculate storage volume
        for ii = 1:Nhot; vol = vol + HT(ii).tank_volA + HT(ii).tank_volB ; end
        rhoE = Exergy_into_hot_tanks/fact/vol*1e3; %kWh/m3
        
    case {2,5,7} % Heat engine only
        Heat_from_tanks   = (HT.A(2).H - HT.A(end).H) + (HT.B(2).H - HT.B(end).H);
        Exergy_from_tanks = (HT.A(2).B - HT.A(end).B) + (HT.B(2).B - HT.B(end).B);
        
        Heat_rejected = QE_dis + W_fan_dis + WL_gen_dis;
        E_net_dis = E_out_dis + W_fan_dis + W_pmp_dis  + WL_gen_dis;
        Total_Work_lost = Exergy_from_tanks - E_net_dis;
        
        First_law_error = (Heat_from_tanks - E_net_dis + Heat_rejected)/(E_net_dis);
        Second_law_error = (Total_loss + Total_Work_lost)/Total_Work_lost;
        
        chi_tot      = (W_out_dis)/Exergy_from_tanks;
        chi_tot_para = (E_net_dis)/Exergy_from_tanks;
        EFF          = W_out_dis/QH_dis;
        EFF_para     = E_net_dis/QH_dis;
        rhoE         = E_net_dis/fact/(HT.A(end).V)*1e3; %kWh/m3
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
    case {0,1,3,4,6,20,22,24}
        Exergy_in = W_in_chg;
    case {2,5,7}
        Exergy_in = Exergy_from_tanks;
end

if WM==1
    fprintf(1,'\n\n');
    fprintf(1,'PTES CYCLE\n');
    fprintf(1,'----------\n');
    
    % Print working fluid states
    fprintf(1,'Gas states:\n');
    for ip = 1:numel(i_gas)
        iL = i_gas(ip);
        if Load.mode == 20
            print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
            print_states(gas2,iL,1:gas2.Nstg(iL)+1,Load);
        else
            print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        end
    end
    for ip = 1:numel(i_ran)
        iL = i_ran(ip);
        print_states(steam,iL,1:steam.Nstg(iL)+1,Load);
    end
    
    % Print hot streams
    fprintf(1,'\nHot fluid streams:\n');
    fprintf(1,'-->%s\n',fluidH(1).name);
    for ip = 1:numel(i_act)
        iL = i_act(ip);
        for iH=1:numel(fluidH)
            print_states(fluidH(iH),iL,1:fluidH(iH).Nstg(iL)+1,Load)
        end
    end
    
    if any(Load.mode == [20,22,24])
        % Print medium streams
        fprintf(1,'\nMedium fluid streams:\n');
        fprintf(1,'-->%s\n',fluidM(1).name);
        for ip = 1:numel(i_act)
            iL = i_act(ip);
            for iM=1:numel(fluidM)
                print_states(fluidM(iM),iL,1:fluidM(iM).Nstg(iL)+1,Load)
            end
        end
    end
        
    if ~any(Load.mode == [20,22])
        % Print cold streams
        fprintf(1,'\nCold fluid streams:\n');
        fprintf(1,'-->%s\n',fluidC(1).name);
        for ip = 1:numel(i_act)
            iL = i_act(ip);
            for iC=1:numel(fluidC)
                print_states(fluidC(iC),iL,1:fluidC(iC).Nstg(iL)+1,Load)
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
    
    if any(Load.mode == [20,22,24])
        
        % Print medium tanks
        for ii = 1 : Nmid
            fprintf(1,'\nMedium tank #%2i\n',ii);
            fprintf(1,'%10s %10s %13s %13s %13s %13s %8s ','A.T [K]','A.M [kg*1e6]','A.H [MWh]','B.T [K]','B.M [kg*1e6]','B.H [MWh]','state'); fprintf(1,'\n');
            for i0=1:(Load.num+1)
                fprintf(1,'%10.4g %13.3f %13.3f %10.4g %13.3f %13.3f %8d\n', MT(ii).A(i0).T,MT(ii).A(i0).M/1e6,MT(ii).A(i0).H/fact,MT(ii).B(i0).T,MT(ii).B(i0).M/1e6,MT(ii).B(i0).H/fact,i0)
            end
        end
        
        % Print liquid air tank
        fprintf(1,'\nLiquid air tank #%2i\n',ii);
        fprintf(1,'%10s %13s %13s %8s','A.T [K]','A.M [kg*1e6]','A.H [MWh]','state'); fprintf(1,'\n');
        for i0=1:(Load.num+1)
            fprintf(1,'%10.4g %13.3f %13.3f %8d\n', LAT.A(i0).T,LAT.A(i0).M/1e6,LAT.A(i0).H/fact,i0)
        end
        
        % Print atmospheric tank
        fprintf(1,'\nAtmospheric tank #%2i\n',ii);
        fprintf(1,'%10s %13s %13s %8s ','A.T [K]','A.M [kg*1e6]','A.H [MWh]','state'); fprintf(1,'\n');
        for i0=1:(Load.num+1)
            fprintf(1,'%10.4g %13.3f %13.3f %8d\n', AT.A(i0).T,AT.A(i0).M/1e6,AT.A(i0).H/fact,i0)
        end
        
    end
    
    if ~any(Load.mode == [20,22])
        
        % Print cold tanks
        for ii = 1 : Ncld
            fprintf(1,'\nCold tank #%2i\n',ii);
            fprintf(1,'%10s %10s %13s %13s %13s %13s %8s ','A.T [K]','A.M [kg*1e6]','A.H [MWh]','B.T [K]','B.M [kg*1e6]','B.H [MWh]','state'); fprintf(1,'\n');
            for i0=1:(Load.num+1)
                fprintf(1,'%10.4g %13.3f %13.3f %10.4g %13.3f %13.3f %8d\n', CT(ii).A(i0).T,CT(ii).A(i0).M/1e6,CT(ii).A(i0).H/fact,CT(ii).B(i0).T,CT(ii).B(i0).M/1e6,CT(ii).B(i0).H/fact,i0)
            end
        end
        
    end
    fprintf(1,'\n');
end

if WM == 1
    fprintf(1,'\n\n');
    fprintf(1,'FIRST LAW BALANCE\n');
    fprintf(1,'-----------------\n');
    
    switch Load.mode
        case {0,1,3,4,6,20,22,24}
            fprintf(1,'CHARGE\n');
            fprintf(1,'Average power input:     %8.2f MW\n',-E_net_chg/t_chg/1e6);
            fprintf(1,'Total charge time:       %8.2f h\n',t_chg/3600);
            fprintf(1,'Electricity input:       %8.2f MWh\n',-E_net_chg/fact);
            fprintf(1,'Work input:              %8.2f MWh\n',-W_in_chg/fact);
            fprintf(1,'Heat to hot tanks:       %8.2f MWh\n',-QH_chg/fact);
            if any(Load.mode == [20,22,24])
                fprintf(1,'Heat from working fluid: %8.2f MWh\n',QWF_chg/fact);
            end
            if ~any(Load.mode == [20,22])
                fprintf(1,'Heat from cold tanks:    %8.2f MWh\n',QC_chg/fact);
            end
            fprintf(1,'Heat rejected:           %8.2f MWh\n',-(QE_chg+W_fan_chg+WL_mot_chg)/fact);
            fprintf(1,'NET:                     %8.2f MWh\n\n',Net_chg/fact);
    end
    
    switch Load.mode
        case {0,2,3,4,5,6,7,20,22,24}
            fprintf(1,'DISCHARGE\n');
            fprintf(1,'Average power output:    %8.2f MW\n',E_net_dis/t_dis/1e6);
            fprintf(1,'Total discharge time:    %8.2f h\n',t_dis/3600);
            fprintf(1,'Electricity output:      %8.2f MWh\n',E_net_dis/fact);
            fprintf(1,'Work output:             %8.2f MWh\n',W_out_dis/fact);
            fprintf(1,'Heat from hot tanks:     %8.2f MWh\n',QH_dis/fact);
            if any(Load.mode == [20,22,24])
                fprintf(1,'Heat to working fluid:   %8.2f MWh\n',-QWF_dis/fact);
            end
            if ~any(Load.mode == [20,22])
                fprintf(1,'Heat to cold tanks:      %8.2f MWh\n',-QC_dis/fact);
            end
            fprintf(1,'Heat rejected:           %8.2f MWh\n',-(QE_dis+W_fan_dis+WL_gen_dis)/fact);
            fprintf(1,'NET:                     %8.2f MWh\n\n',Net_dis/fact);
    end
    
    switch Load.mode
        case 3
            fprintf(1,'DISCHARGE DETAILS:\n');
            fprintf(1,'Rankine cycle average power output:          %8.2f MW\n',W_out_dis/t_dis/1e6);
            fprintf(1,'Rankine cycle power output NO cold stores:   %8.2f MW\n',W_out_disNC/(t_dis-t_disRC)/1e6);
            fprintf(1,'Rankine cycle power output WITH cold stores: %8.2f MW\n',W_out_disRC/t_disRC/1e6);
            
            fprintf(1,'Rankine cycle average efficiency:            %8.2f %%\n',HEeff*100);
            fprintf(1,'Rankine cycle efficiency NO cold stores:     %8.2f %%\n',HEeffNC*100);
            fprintf(1,'Rankine cycle efficiency using cold stores:  %8.2f %%\n\n',HEeffRC*100);
        case 6
            fprintf(1,'EFFICIENCIES\n');
            fprintf(1,'Solar conversion efficiency:     %8.1f %%\n',SOLeff*100);
            fprintf(1,'Heat engine efficiency:          %8.1f %%\n',HEeff*100);
            fprintf(1,'Net efficiency:                  %8.1f %%\n',NETeff*100);
            fprintf(1,'Exergetic efficiency:            %8.1f %%\n\n',EXeff*100);
    end
    
    switch Load.mode
        case {0,3,4,6,20,22,24}
            fprintf(1,'COP:                                 %8.2f \n',QH_chg/W_in_chg);
            fprintf(1,'COP (including parasitics):          %8.2f \n',QH_chg/E_net_chg);
            fprintf(1,'Heat engine efficiency:              %8.2f %%\n',W_out_dis/QH_dis*100);
            fprintf(1,'Heat engine eff. (inc. parasitics):  %8.2f %%\n\n',E_net_dis/QH_dis*100);
            fprintf(1,'Round trip efficiency:               %8.2f %%\n',chi_PTES*100);
            fprintf(1,'Round trip eff. (inc. parasitics):   %8.2f %%\n\n',chi_PTES_para*100);
            fprintf(1,'Exergy density:                      %8.2f kWh/m3\n\n',rhoE);
            
            fprintf(1,'STORAGE MEDIA\n');
            for ii = 1 : Nhot
                fprintf(1,'%18s volume:%8.2f m3/MWh\n',fluidH(ii).name,HT(ii).fluid_volB/(W_out_dis/fact));
                fprintf(1,'%18s volume:%8.2f m3\n',fluidH(ii).name,HT(ii).fluid_volB);
                fprintf(1,'%18s mass:  %8.2f tons/MWh\n\n',fluidH(ii).name,HT(ii).fluid_mass/(W_out_dis/fact)/1e3);
            end
            for ii = 1 : Nmid
                fprintf(1,'%18s volume:%8.2f m3/MWh\n',fluidM(ii).name,MT(ii).fluid_volB/(W_out_dis/fact));
                fprintf(1,'%18s volume:%8.2f m3\n',fluidM(ii).name,MT(ii).fluid_volB);
                fprintf(1,'%18s mass:  %8.2f tons/MWh\n\n',fluidM(ii).name,MT(ii).fluid_mass/(W_out_dis/fact)/1e3);
            end
            for ii = 1 : Ncld
                fprintf(1,'%18s volume:%8.2f m3/MWh\n',fluidC(ii).name,CT(ii).fluid_volB/(W_out_dis/fact));
                fprintf(1,'%18s volume:%8.2f m3\n',fluidC(ii).name,CT(ii).fluid_volB);
                fprintf(1,'%18s mass:  %8.2f tons/MWh\n\n',fluidC(ii).name,CT(ii).fluid_mass/(W_out_dis/fact)/1e3);
            end
            for ii = 1 : Nlat
                fprintf(1,'%18s volume:%8.2f m3/MWh\n','Liquid air',LAT(ii).fluid_volA/(W_out_dis/fact));
                fprintf(1,'%18s volume:%8.2f m3\n','Liquid air',LAT(ii).fluid_volA);
                fprintf(1,'%18s mass:  %8.2f tons/MWh\n\n','Liquid air',LAT(ii).fluid_mass/(W_out_dis/fact)/1e3);
            end
            
        case 1
            fprintf(1,'Exergetic efficiency:                         %8.1f %%\n',-chi_tot*100);
            fprintf(1,'Exergetic efficiency (inc. parasitics):       %8.1f %%\n',-chi_tot_para*100);
            fprintf(1,'Exergetic efficiency (cold reject):           %8.1f %%\n',-chi_hot*100);
            fprintf(1,'Coefficient of Performance:                   %8.3f\n',COP);
            fprintf(1,'Coefficient of Performance (inc. parasitics): %8.3f\n\n',COP_para);
            fprintf(1,'Exergy density (hot tanks):                   %8.3f kWh/m3\n',rhoE);
        case {2,5,7}
            fprintf(1,'Exergetic efficiency:             %8.3f %%\n',chi_tot*100);
            fprintf(1,'Exergetic eff. (inc. parasitics): %8.3f %%\n',chi_tot_para*100);
            fprintf(1,'First law efficiency:             %8.3f %%\n',EFF*100);
            fprintf(1,'First law eff. (inc. parasitics): %8.3f %%\n\n',EFF_para*100);
            fprintf(1,'Exergy density:                   %8.3f kWh/m3\n\n',rhoE);
    end
    
    fprintf(1,'First law error:   %8.5f %%\n',First_law_error*100);
    fprintf(1,'Second law error:  %8.5f %%\n\n',Second_law_error*100);
end

% Check for the case where storages are not returned to their original
% temperature. If stores are above T0 when discharged, then they must be
% returned to their original temp or greater. If stores are below T0 when
% discharged, then they must be returned to original temp or lower.
heater_in = 0;
chi_PTES_true = chi_PTES_para ;
if any(Load.mode ==[0,4,6,20,22,24])
    % Hot tanks
    problem = 0 ;
    for ii = 1 : Nhot
       if HT(ii).A(1).T >= T0
           if HT(ii).A(end).T < HT(ii).A(1).T - 1.0; problem = 1 ; end
       else
           if HT(ii).A(end).T > HT(ii).A(1).T + 1.0; problem = 1 ; end
       end
    end
    if problem
        warning('Unsustainable discharge of a hot reservoir!')
    end
    % Medium tanks
    problem = 0 ;
    for ii = 1 : Nmid
       if MT(ii).A(1).T >= T0
           if MT(ii).A(end).T < MT(ii).A(1).T - 1.0; problem = 1 ; end
       else
           if MT(ii).A(end).T > MT(ii).A(1).T + 1.0; problem = 1 ; end
       end
    end
    if problem
        warning('Unsustainable discharge of a medium reservoir!')
    end
    % Cold tanks
    problem = 0;
    for ii = 1 : Ncld
       if CT(ii).A(1).T >= T0
           if CT(ii).A(end).T < CT(ii).A(1).T - 1.0; problem = 1 ; end
       else
           if CT(ii).A(end).T > CT(ii).A(1).T + 1.0; problem = 1 ; end
       end
    end
    if problem
        warning('Unsustainable discharge of a cold reservoir!')
    end
    
    % If the hot fluid gets cooled down to a temperature below its original
    % value, calculate the heat required to boost it back. Then find the
    % 'true' round-trip efficiency assuming this heat is provided by an
    % electrical heater
    if fluidH(1).state(end,3).T < fluidH(1).state(1,1).T
        warning('Unsustainable discharge of a hot reservoir!')
        warning('Hot fluid cooled down too much in discharge. Calculating new roundtrip efficiency')
        heater_in =  (fluidH.state(1,1).h - fluidH.state(end,3).h) * fluidH.state(end,3).mdot * t_dis;
        chi_PTES_true = -(E_net_dis - heater_in) / E_net_chg ;
        fprintf(1,'TRUE round trip eff. (inc. heating):   %8.2f %%\n\n',chi_PTES_true*100);
    end
end

switch Load.mode
    case {20,22,24}
        
    otherwise
        
        % PRINT HEXs
        % If the heat exchanger was employed with the 'eff' or 'DT' modes, the
        % required geometry is computed now
        for ii = 1 : numel(HX)
            if any(strcmp(HX(ii).model,{'eff','DT'})) && (~HX(ii).Lgeom_set)
                HX(ii)   = hex_set_geom(HX(ii));
            end
        end
        %%{
        fprintf('Heat exchanger summary\n');
        if Load.mode==3
            print_hexs(HX,i_chg,'Charge:');
            print_hexs(HX,Load.ind(Load.type == 'ran' & logical(Load.options.useCold)),...
                'Discharge using cold stores:');
            print_hexs(HX,Load.ind(Load.type == 'ran' & ~logical(Load.options.useCold)),...
                'Discharge without cold stores:');
        else
            print_hexs(HX,i_chg,'Charge:');
            print_hexs(HX,i_dis,'Discharge:');
        end
        %%}
        
        if Lreadload
            ENERGY_LOAD
        end
        
end