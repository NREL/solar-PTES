classdef CYCLE 
    % CYCLE is a construct for the overall cycle metrics
    
    properties
        name    % name
        note    % Any other details
        WinC    % Work input in charge
        WoutC   % Work output in charge
        WnetC   % Net work in charge
        WinD    % Work input in discharge
        WoutD   % Work output in discharge
        WnetD   % Net work in discharge
        Qsol    % Solar heat added
        QhotC   % Heat to hot storage during charging
        QcldC   % Heat to cold storage during charging
        QhotD   % Heat to hot storage during charging
        QcldD   % Heat to cold storage during charging
        etaHS   % Recovery efficiency of HS
        etaCS   % Recovery efficiency of CS
        Qin     % Heat input to the cycle
        Qrej    % Heat rejected from power cycle
        QrejH1  % Heat rejected from HS1
        QrejH2  % Heat rejected from HS2
        QrehC1  % Make-up heat for CS1
        QrehC2  % Make-up heat for CS1
        COP     % Coefficient of performance of charge
        eta     % Heat engine efficiency of discharge
        eta_NP  % Heat engine effieincy with no reheat penalty
        RTeff   % Round-trip efficiency,
        RTeff_NP
        Exsol   % Exergy added by solar heat
        RTeff2
        eta_net % Net efficiency
        powden  % The power density. Net work * minimum density (at expander outlet)
        
        dRCP_qh % Heat available on hot side of recuperation during discharge
        dRCP_qc % Heat required on cold side of recuperation during discharge
        
    end
    
    % Constructor function
    methods
        function obj = CYCLE(cyc_dat) 
           obj.name = cyc_dat.name ;
        end
    end
    
end