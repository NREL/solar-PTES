% This script calculates the work and efficiencies for each load cycle
% This is meant for off-design calculations or consecutive load cycles

% First calculate some of the design parameters

% Net work input to charging cycle
Wchg0 = 0 ;
for ii = 1:length(CCMP)
    Wchg0 = Wchg0 - CCMP(ii).W0 ;
end
for ii = 1:length(CEXP)
    Wchg0 = Wchg0 + CEXP(ii).W0 ;
end
for ii = 1:length(CFAN)
    Wchg0 = Wchg0 - CFAN(ii).W0 ;
end
for ii = 1:length(CPMP)
    Wchg0 = Wchg0 - CPMP(ii).W0 ;
end

% Heat into hot storage during charge
QH_chg0 = 0;
for ii = 1:length(HX)
    if strcmp(HX(ii).name,'hot')
        QH_chg0 = QH_chg0 + HX(ii).Q(1,1) ;
    end
end

COP0 = QH_chg0 / (Wchg0 * 3600) ;

% Net work output from discharging cycle
Wdis0 = 0 ;
for ii = 1:length(DCMP)
    Wdis0 = Wdis0 - DCMP(ii).W0 ;
end
for ii = 1:length(DEXP)
    Wdis0 = Wdis0 + DEXP(ii).W0 ;
end
for ii = 1:length(DFAN)
    Wdis0 = Wdis0 - DFAN(ii).W0 ;
end
for ii = 1:length(DPMP)
    Wdis0 = Wdis0 - DPMP(ii).W0 ;
end

% Heat out of hot storage during discharge
QH_dis0 = 0;
for ii = 1:length(HX)
    if strcmp(HX(ii).name,'hot')
        QH_dis0 = QH_dis0 - HX(ii).Q(2,1) ;
    end
end

HEeff0 = Wdis0 * 3600. / QH_dis0 ;
RTeff0 = -Wdis0 / Wchg0 ;


% Now do the same for each load cycle
WchgL = zeros(Load.num,1) ;
WdisL = zeros(Load.num,1) ;
QH_chgL = zeros(Load.num,1) ;
QH_disL = zeros(Load.num,1) ;
COPL     = zeros(Load.num,1) ;
HEeffL   = zeros(Load.num,1) ;

for jj = 1 : Load.num
    % Net work input to charging cycle
    for ii = 1:length(CCMP)
        WchgL(jj) = WchgL(jj) + CCMP(ii).W(jj) ;
    end
    for ii = 1:length(CEXP)
        WchgL(jj) = WchgL(jj) + CEXP(ii).W(jj) ;
    end
    for ii = 1:length(CFAN)
        WchgL(jj) = WchgL(jj) + CFAN(ii).W(jj) ;
    end
    for ii = 1:length(CPMP)
        WchgL(jj) = WchgL(jj) + CPMP(ii).W(jj) ;
    end
    
    % Heat into hot storage during charge
    for ii = 1:length(HX)
        if strcmp(HX(ii).name,'hot')
            QH_chgL(jj) = QH_chgL(jj) + HX(ii).Q(jj,1) ;
        end
    end
    
    COPL(jj) = QH_chgL(jj) / WchgL(jj) ;
    
    % Net work output from discharging cycle
    for ii = 1:length(DCMP)
        WdisL(jj) = WdisL(jj) + DCMP(ii).W(jj) ;
    end
    for ii = 1:length(DEXP)
        WdisL(jj) = WdisL(jj) + DEXP(ii).W(jj) ;
    end
    for ii = 1:length(DFAN)
        WdisL(jj) = WdisL(jj) + DFAN(ii).W(jj) ;
    end
    for ii = 1:length(DPMP)
        WdisL(jj) = WdisL(jj) + DPMP(ii).W(jj) ;
    end
    
    % Heat out of hot storage during discharge
    for ii = 1:length(HX)
        if strcmp(HX(ii).name,'hot')
            QH_disL(jj) = QH_disL(jj) - HX(ii).Q(jj,1) ;
        end
    end
    
    HEeffL(jj) = WdisL(jj) / QH_disL(jj) ;

    
end


RTeffAVE = sum(WdisL) / sum(WchgL) ;