classdef TS
    % TS is a construct for making T-s diagrams
    
    properties
        name    % name
        note    % Any other details
        TC      % Charging temps
        sC      % Charging entropies
        SC_pts  % Charging entropy points
        sC_l    % Charging entropies (lower pressure)
        sC_h    % Charging entropies (higher pressure)
        PC      % Charging pressures
        PC_l    % Charging pressures (lower pressure)
        PC_h    % Charging pressures (higher pressure)
        TD      % Discharging temps
        sD      % Discharging entropies curves
        SD_pts  % Discharging entropy points
        PD      % Discharging pressures
        PD_l    % Discharging pressures (lower pressure)
        PD_h    % Discharging pressures (higher pressure)
        CHG     % Charging point properties
        DIS     % Discharging point properties
        type    % 'Joule' for all electric heaters *Could use this to sort through the different cycles?*
    end
    
    % Constructor function
    methods
        function obj = TS(TS_dat1,TS_dat2,TS_dat3,TS_dat4,TS_dat5)
            
            if TS_dat5 == 'Sat'
                   %plot the dome
                   %Google the fluid's critical Temp!
                   Tcrit = 304.25 ;
                   T_sat = 274 : Tcrit ;
                   for x = 1 : length(T_sat)
                        s_satl(x) = CoolProp.PropsSI('S','Q',0,'T',T_sat(x),TS_dat4) ;
                        s_satg(x) = CoolProp.PropsSI('S','Q',1,'T',T_sat(x),TS_dat4) ;
                   end
                   plot(s_satl,T_sat - 273.15,'k',s_satg,T_sat - 273.15,'k','Linewidth',1.5), hold on
            else
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %DISCHARGING
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            obj.DIS = TS_dat2((2:length(TS_dat2(:,1))),:) ;
               
              for y = 1:length(obj.DIS(:,1))/2
                   obj.TD     = [obj.TD, obj.DIS(2*y-1,3)] ;
                   obj.PD     = [obj.PD, obj.DIS(2*y-1,4)] ;
                   obj.SD_pts = [obj.SD_pts, obj.DIS(2*y-1,6)] ;
              end
              
              obj.TD     = [obj.TD, obj.DIS(length(obj.DIS(:,1)),3)] + 273.15 ;
              obj.PD     = [obj.PD, obj.DIS(length(obj.DIS(:,1)),4)] * 1.e5  ;
              obj.SD_pts = [obj.SD_pts, obj.DIS(length(obj.DIS(:,1)),6)] * 1.e3 ; 
               
              ND = length(obj.TD) ;
              
              for x = 1 : ND - 1
                    if max([obj.PD(x), obj.PD(x+1)]) / min([obj.PD(x), obj.PD(x+1)]) > 1.1
                        p_D(x) = plot(obj.SD_pts(x:x+1), obj.TD(x:x+1) - 273.15, 'r','Linewidth',1.5) ; hold on
                    else
                        TD  = min(obj.TD(x:x+1)) : max(obj.TD(x:x+1)) ;
                        if obj.PD(x) - obj.PD(x+1) == 0
                            PD = obj.PD(x) ; 
                            s_D = [] ;
                            for   y = 1 : length(TD)
                                s_D(y) = CoolProp.PropsSI('S','P',PD,'T',TD(y),TS_dat4) ;
                            end
                        else
                        PD  = min(obj.PD(x:x+1)) : length(min(obj.PD(x:x+1)) : max(obj.PD(x:x+1))) / length(TD) : max(obj.PD(x:x+1)) ;
                        s_D = [] ;
                        for   y = 1 : length(TD)
                            s_D(y) = CoolProp.PropsSI('S','P',PD(y),'T',TD(y),TS_dat4) ;
                        end
                        end
                        p_D(x) = plot(s_D, TD - 273.15, 'r','Linewidth',1.5) ; hold on
                    end
               end
               
               if TS_dat3 == 0
               
               else
                    if TS_dat5 == 'Sat'
                    p_D(x) = plot([obj.SD_pts(ND),obj.SD_pts(ND-1)],[TD - 273.15,TD - 273.15],'r','Linewidth',1.5) ;
                    else
                    if max([obj.PD(1), obj.PD(ND)]) / min([obj.PD(1), obj.PD(ND)]) > 1.1
                        p_D(x) = plot([obj.SD_pts(1),obj.SD_pts(ND)], [obj.TD(1) - 273.15,obj.TD(ND) - 273.15], 'r','Linewidth',1.5) ; hold on
                    else
                        TD  = min([obj.TD(1),obj.TD(ND)]) : max([obj.TD(1),obj.TD(ND)]) ;
                        PD  = min([obj.PD(1),obj.PD(ND)]) : length(min([obj.PD(1),obj.PD(ND)]) : max([obj.PD(1),obj.PD(ND)])) / length(TD) : max([obj.PD(1),obj.PD(ND)]) ;
                        s_D = [];
                        for   y = 1 : length(TD)
                            s_D(y) = CoolProp.PropsSI('S','P',PD(y),'T',TD(y),TS_dat4) ;
                        end
                        p_D(x) = plot(s_D, TD - 273.15, 'r','Linewidth',1.5) ; hold on
                    end
                    end
               
               end  
            
            if TS_dat1 == 0
            else
                
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
               %CHARGING
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               obj.CHG = TS_dat1((2:length(TS_dat1(:,1))),:) ;

               for x = 1:length(obj.CHG(:,1))/2
                    obj.TC     = [obj.TC, obj.CHG(2*x-1,3)] ;
                    obj.PC     = [obj.PC, obj.CHG(2*x-1,4)] ;
                    obj.SC_pts = [obj.SC_pts, obj.CHG(2*x-1,6)] ;
               end
   
               obj.TC     = [obj.TC, obj.CHG(length(obj.CHG(:,1)),3)] + 273.15;
               obj.PC     = [obj.PC, obj.CHG(length(obj.CHG(:,1)),4)] * 1.e5 ;
               obj.SC_pts = [obj.SC_pts, obj.CHG(length(obj.CHG(:,1)),6)] * 1.e3;

               NC = length(obj.TC) ;
               
               for x = 1 : NC - 1
                    if max([obj.PC(x), obj.PC(x+1)]) / min([obj.PC(x), obj.PC(x+1)]) > 1.1
                        p_C(x) = plot(obj.SC_pts(x:x+1), obj.TC(x:x+1) - 273.15, '--b','Linewidth',1.5) ; hold on
                    else
                        TC  = min(obj.TC(x:x+1)) : max(obj.TC(x:x+1)) ;
                        if obj.PC(x) - obj.PC(x+1) == 0
                            PC = obj.PC ;
                            s_C = [] ;
                            for   y = 1 : length(TC)
                                s_C(y) = CoolProp.PropsSI('S','P',PC,'T',TC(y),TS_dat4) ;
                            end
                        else
                        PC  = min(obj.PC(x:x+1)) : length(min(obj.PC(x:x+1)) : max(obj.PC(x:x+1))) / length(TC) : max(obj.PC(x:x+1)) ;
                        s_C = [] ;
                        for   y = 1 : length(TC)
                            s_C(y) = CoolProp.PropsSI('S','P',PC(y),'T',TC(y),TS_dat4) ;
                        end
                        end
                        p_C(x) = plot(s_C, TC - 273.15, '--b','Linewidth',1.5) ; hold on
                    end
               end
               
               
               
               if TS_dat3 == 0
               else
              
               if TS_dat5 == 'Sat'
                    p_C(x) = plot([obj.SC_pts(NC),obj.SC_pts(NC-1)],[TC - 273.15,TC - 273.15],'Linewidth',1.5) ;
               else
                    if max([obj.PC(1), obj.PC(NC)]) / min([obj.PC(1), obj.PC(NC)]) > 1.1
                        p_C(x) = plot([obj.SC_pts(1),obj.SC_pts(NC)], [obj.TC(1) - 273.15,obj.TC(NC) - 273.15], '--b','Linewidth',1.5) ; hold on
                    else
                        TC  = min([obj.TC(1),obj.TC(NC)]) : max([obj.TC(1),obj.TC(NC)]) ;
                        PC  = min([obj.PC(1),obj.PC(NC)]) : length(min([obj.PC(1),obj.PC(NC)]) : max([obj.PC(1),obj.PC(NC)])) / length(TC) : max([obj.PC(1),obj.PC(NC)]) ;
                        s_C = [];
                        for   y = 1 : length(TC)
                            s_C(y) = CoolProp.PropsSI('S','P',PC(y),'T',TC(y),TS_dat4) ;
                        end
                            p_C(x) = plot(s_C, TC - 273.15, '--b','Linewidth',1.5) ;
                    end
               end              
               end
            end
            
               if TS_dat1 == 0
                   scatter(obj.SD_pts, obj.TD  - 273.15, 'k')
               else  
               legend([p_C(1) p_D(1)],'Charge','Discharge','Location','southeast','AutoUpdate', 'off')
%                scatter([obj.SC_pts, obj.SD_pts], [obj.TC - 273.15, obj.TD - 273.15], 'k','x')
               scatter(obj.SC_pts, obj.TC - 273.15, 'k','x')
               scatter(obj.SD_pts, obj.TD - 273.15, 'k')
               end
               grid()
        end
    end
end