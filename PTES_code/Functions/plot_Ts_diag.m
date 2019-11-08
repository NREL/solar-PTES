function [ref] = plot_Ts_diag( fluid, Load, load_types, fignum, num, line_spec, point_spec, Celcius, Sat_curve)
% PLOT A TS DIAGRAM OF THE CYCLE.
% Use the fluid structure to obtain the fluid states at each point along
% the cycle. The Load structure is used to obtained the index iL of the
% load period which is to be plot, as specified in load_type. Only the
% first load period that matches load_types is plotted (i.e. if there are
% several 'chg' cycles within Load.type and load_types = {chg}, only the
% first one is plot).
% Fignum is the figure number. num is the number of points used to plot
% each stage. line_spec and point_spec give the graphic specifications of
% the lines and points to be plot.
% Celcius is a logical variable that determines whether to use degrees C
% (rather than K) in the plots. Sat_curve is another logical variable that
% determines whether the saturation curve of the fluid should also be
% plotted.
% The "ref" value that is returned can be used when creating the legend.

% Usage example:
% pl1 = plot_stage(gas,Load,{'chg','chgCO2'},1,100,'k-','k-o',true,false);


figure(fignum)

% Set temperature in Celcius?
switch Celcius
    case false
        K_C = 0;
    case true
        K_C = -273.15;
end

% Screen the different Load periods
for iL=1:Load.num
    
    % Plot the first Load period that matches load_types, then break
    if any(strcmp(Load.type(iL),load_types))
        
        for i=1:fluid.Nstg(iL)
            
            % Import fluid.state and fluid.stage
            state_in  = fluid.state(iL,i);
            state_out = fluid.state(iL,i+1);
            stage     = fluid.stage(iL,i);
            
            % Create enthalpy and pressure arrays
            if any(strcmp(stage.type,{'hex','hex_reject','regen','split','mixing','separate'}))
                type = 1;
                
            elseif any(strcmp(stage.type,{'comp','exp'}))
                type = 3;
                
            elseif strcmp(stage.type,{'end'})
                type = 0;
                
            else
                warning('not implemented');
            end
            
            switch type
                case 1 % nearly-isobaric line
                    p_vect = logspace(log10(state_in.p),log10(state_out.p),num);
                    h_vect = linspace(state_in.h,state_out.h,num);
                    
                case 2 % polytropic line
                    phi   = log(state_out.p/state_in.p)/log(state_out.h/state_in.h);
                    h_vect = linspace(state_in.h,state_out.h,num);
                    p_vect = state_in.p*(h_vect./state_in.h).^phi;
                    
                case 3 % straight line (on T-s diagram)
                    T_vect = linspace(state_in.T,state_out.T,num);
                    s_vect = linspace(state_in.s,state_out.s,num);
            end
            
            switch type
                case {1,2}
                    Pcrit = CP1(0,0,0,'Pcrit',fluid.handle);
                    
                    if any(p_vect>0.99*Pcrit & p_vect<1.00*Pcrit)
                        [s1,T1,~,~,~] = CP5('HmassP_INPUTS',h_vect,p_vect*1.01,'S','T','P','P','P',fluid.handle);
                        [s2,T2,~,~,~] = CP5('HmassP_INPUTS',h_vect,p_vect*0.99,'S','T','P','P','P',fluid.handle);
                        s_vect = 0.5*(s1 + s2);
                        T_vect = 0.5*(T1 + T2);
                    else
                        [s_vect,T_vect,~,~,~] = CP5('HmassP_INPUTS',h_vect,p_vect,'S','T','P','P','P',fluid.handle);
                    end
                    
                    a0 = (s_vect == 0 & T_vect == 0); %logical array of zero elements
                    if any(a0)
                        s_vect(a0) = NaN;
                        T_vect(a0) = NaN;
                        s_vect = fillmissing(s_vect,'pchip');
                        T_vect = fillmissing(T_vect,'pchip');
                    end
            end
            
            if type ~= 0
                plot(s_vect/1e3,T_vect+K_C,line_spec,'LineWidth',2); hold on;
                ref = plot([fluid.state(iL,i).s]/1e3,[fluid.state(iL,i).T]+K_C,point_spec,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
            end            
            
        end
        
        % Obtain minimum and maximum values of entropy to be plotted
        xv   = [fluid.state(iL,1:fluid.Nstg(iL)).s];
        xmin = min(xv)/1e3;
        xmax = max(xv)/1e3;
        
        break
        
    end
end

if Sat_curve    
    % Plot saturation curve
    [T_vect_c,sL,sG,sQ] = plot_sat_curve(fignum,num,10,fluid.handle);
    for i0=1:size(sQ,1)
        plot(sQ(i0,:)/1e3,T_vect_c+K_C,'-k','LineWidth',1); hold on;
    end
    plot(sL/1000,T_vect_c+K_C,'-k');
    plot(sG/1000,T_vect_c+K_C,'-k');
end

% Set x axis
Dx = xmax-xmin;
xlim([(xmin-0.1*Dx) (xmax+0.1*Dx)]);
xlabel('Specific Entropy [kJ/kg.K]');

% Set y axis
switch Celcius
    case false
        ylabel('Temperature [K]');
        ylim([0 1000])
        yticks(0:100:1000)
    case true
        ylabel('Temperature [$$^{\circ}$$C]');
        ylim([-150 700])        
end
grid off;

end