% Plot T-s diagrams

% Set properties of storage temperature lines;
NameArray  = {'LineWidth','LineStyle'};
ValueArray = {0.8,'-'};
plot2 = false ; % Plot charge and discharge on separate figures

switch Load.mode
    case {0,4,6} % PTES
        
        if ~plot2
            % Plot states in Ts diagram
            pl1 = plot_Ts_diag(gas,Load,{'chg','chgCO2','chgTSCO2'},1,100,'k-','k-o',true,false);
            pl2 = plot_Ts_diag(gas,Load,{'dis','disCO2','disTSCO2'},1,100,'k:','k:s',true,true);
        
            % Plot storage tank temperatures
            plot_T_storage(HT(1),Load,{'chg','chgCO2','chgTSCO2'},1,{'r','r'},true,NameArray,ValueArray);
            plot_T_storage(CT(1),Load,{'chg','chgCO2','chgTSCO2'},1,{'b','b'},true,NameArray,ValueArray);
            plot_T_storage(HT(1),Load,{'dis','disCO2','disTSCO2'},1,{'r','r'},true,NameArray,ValueArray);
            plot_T_storage(CT(1),Load,{'dis','disCO2','disTSCO2'},1,{'b','b'},true,NameArray,ValueArray);
        
            % Set legends
            figure(1); hold off;
            legend([pl1, pl2],{'charge','discharge'},'Location','best');
            
        else
            
            % Plot states in Ts diagram
            %{
            pl1 = plot_Ts_diag(gas,Load,{'chg','chgCO2','chgTSCO2'},1,100,'k-','k-o',true,true);
            % Plot storage tank temperatures
            %plot_T_storage(HT(1),Load,{'chg','chgCO2','chgTSCO2'},1,{'r','r'},true,NameArray,ValueArray);
            %plot_T_storage(CT(1),Load,{'chg','chgCO2','disTSCO2'},1,{'b','b'},true,NameArray,ValueArray);
            % Set legends
            figure(1); hold off;
            legend([pl1],{'charge'},'Location','best');
            %}
            pl2 = plot_Ts_diag(gas,Load,{'dis','disCO2','disTSCO2'},2,100,'k:','k:s',true,true);
            %plot_T_storage(HT(1),Load,{'chg','chgCO2','chgTSCO2'},2,{'r','r'},true,NameArray,ValueArray);
            %plot_T_storage(CT(1),Load,{'chg','chgCO2','disTSCO2'},2,{'b','b'},true,NameArray,ValueArray);
            % Set legends
            figure(2); hold off;
            legend(pl2,{'discharge'},'Location','best');
            %}
        end
        
    case 1 % Heat pump only
        
        % Plot states in Ts diagram
        pl1 = plot_Ts_diag(gas,Load,{'chg','chgCO2'},1,100,'k-','k-o',true,false);
        
        % Plot storage tank temperatures
        plot_T_storage(HT(1),Load,{'chg','chgCO2'},1,{'r','r'},true,NameArray,ValueArray);
        plot_T_storage(CT(1),Load,{'chg','chgCO2'},1,{'b','b'},true,NameArray,ValueArray);
        
        % Set legend
        figure(1); hold off;
        legend(pl1,'heat pump','Location','best');
        
    case {2,5} % Heat engine only
        
        % Plot states in Ts diagram
        pl2 = plot_Ts_diag(gas,Load,{'dis','disCO2','rcmpCO2'},2,100,'k:','k:s',true,true);
        
        % Plot storage tank temperatures
        plot_T_storage(HT(1),Load,{'dis','disCO2','rcmpCO2'},2,{'r','r'},true,NameArray,ValueArray);
        plot_T_storage(CT(1),Load,{'dis','disCO2','rcmpCO2'},2,{'b','b'},true,NameArray,ValueArray);
        
        % Set legend
        figure(2); hold off;
        legend(pl2,'heat engine','Location','best');
        
    case 3 %JB-charge, Rankine Discharge
        
        % Plot states in Ts diagram
        pl1 = plot_Ts_diag(gas,  Load,'chg',1,100,'k-','k-o',true,false);        
        pl2 = plot_Ts_diag(steam,Load,'ran',2,100,'k:','k:s',true,true);
        
        % Plot storage tank temperatures
        plot_T_storage(HT,Load,'chg',1,{'r','r'},true,NameArray,ValueArray);
        plot_T_storage(CT,Load,'chg',1,{'b','b'},true,NameArray,ValueArray);
        plot_T_storage(HT,Load,'ran',2,{'r','r'},true,NameArray,ValueArray);
        plot_T_storage(CT,Load,'ran',2,{'b','b'},true,NameArray,ValueArray);
        
        % Set legends
        figure(1); hold off;
        legend(pl1,'charge','Location','best');
        figure(2); hold off;
        legend(pl2,'discharge','Location','best');
        
    case 7 %EH-charge, Rankine Discharge
        
        % Plot states in Ts diagram      
        pl2 = plot_Ts_diag(steam,Load,'ran',2,100,'k:','k:s',true,true);
        
        % Plot storage tank temperatures
        plot_T_storage(HT,Load,'ran',2,{'r','r'},true,NameArray,ValueArray);
        plot_T_storage(CT,Load,'ran',2,{'b','b'},true,NameArray,ValueArray);
        
        % Set legends
        figure(2); hold off;
        legend(pl2,'discharge','Location','best');
        
    case 20 %PTES-LAES Combined Cycle
        if iL==1 || iL>1
            % CHARGE
            % Plot states in Ts diagram
            pl1 = plot_Ts_diag(gas1,Load,'chgCC',1,100,'k-','k-o',true,true);
            pl2 = plot_Ts_diag(gas2,Load,'chgCC',2,100,'k-','k-s',true,true);
            
            % Plot storage tank temperatures
            plot_T_storage(MT,Load,'chgCC',1,{'r','r'},true,NameArray,ValueArray);
            plot_T_storage(HT,Load,'chgCC',2,{'r','r'},true,NameArray,ValueArray);
            
            % Set legends
            figure(1); hold off;
            legend(pl1,'LAES charge','Location','best');
            figure(2); hold off;
            legend(pl2,'PTES charge','Location','best');
        end
        if iL>1
            % DISCHARGE
            % Plot states in Ts diagram
            pl3 = plot_Ts_diag(gas1,Load,'disCC',3,100,'k-','k-o',true,true);
            pl4 = plot_Ts_diag(gas2,Load,'disCC',4,100,'k-','k-s',true,true);
            
            % Plot storage tank temperatures
            plot_T_storage(MT,Load,'disCC',3,{'r','r'},true,NameArray,ValueArray);
            plot_T_storage(HT,Load,'disCC',4,{'r','r'},true,NameArray,ValueArray);
            
            % Set legends
            figure(3); hold off;
            legend(pl3,'LAES discharge','Location','best');
            figure(4); hold off;
            legend(pl4,'PTES discharge','Location','best');
        end
        
    case 22 %Integrated PTES-LAES Combined Cycle
        if iL==1 || iL>1
            % CHARGE
            % Plot states in Ts diagram
            pl1 = plot_Ts_diag(gas,Load,'chgICC',1,100,'k-','k-o',true,true);
            
            % Plot storage tank temperatures
            plot_T_storage(MT,Load,'chgICC',1,{'r','r'},true,NameArray,ValueArray);
            %plot_T_storage(HT,Load,'chgICC',1,{'r','r'},true,NameArray,ValueArray);
            
            % Set legends
            figure(1); hold off;
            %legend(pl1,'charge','Location','best');
        end
        if iL>1
            % DISCHARGE
            % Plot states in Ts diagram
            pl2 = plot_Ts_diag(gas,Load,'disICC',2,100,'k-','k-o',true,true);
            
            % Plot storage tank temperatures
            plot_T_storage(MT,Load,'disICC',2,{'r','r'},true,NameArray,ValueArray);
            %plot_T_storage(HT,Load,'disICC',2,{'r','r'},true,NameArray,ValueArray);
            
            % Set legends
            figure(2); hold off;
            %legend(pl2,'discharge','Location','best');
            %title('Integrated discharging cycle')
        end
        
    case 24 %Integrated PTES-LAES Combined Cycle with Pre-cooling system
        if iL==1 || iL>1
            % CHARGE
            % Plot states in Ts diagram
            pl1 = plot_Ts_diag(gas,Load,'chgICC_PC',1,100,'k-','k-o',true,true);
            
            % Plot storage tank temperatures
            plot_T_storage(HT,Load,'chgICC_PC',1,{'r','r'},true,NameArray,ValueArray);
            plot_T_storage(MT,Load,'chgICC_PC',1,{'r','r'},true,NameArray,ValueArray);
            plot_T_storage(CT,Load,'chgICC_PC',1,{'b','b'},true,NameArray,ValueArray);
            
            % Set legends
            figure(1); hold off;
            %legend(pl1,'charge','Location','best');
        end
        if iL>1
            % DISCHARGE
            % Plot states in Ts diagram
            pl2 = plot_Ts_diag(gas,Load,'disICC_PC',2,100,'k-','k-o',true,true);
            
            % Plot storage tank temperatures
            plot_T_storage(HT,Load,'disICC_PC',2,{'r','r'},true,NameArray,ValueArray);
            plot_T_storage(MT,Load,'disICC_PC',2,{'r','r'},true,NameArray,ValueArray);            
            plot_T_storage(CT,Load,'disICC_PC',2,{'b','b'},true,NameArray,ValueArray);
            
            % Set legends
            figure(2); hold off;
            %legend(pl2,'discharge','Location','best');
        end
                
end