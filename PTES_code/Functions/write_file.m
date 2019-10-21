function [ ] = write_file( fluid, ind, fileID, num )

% Import fluid.state and fluid.stage
fluid_in  = fluid.state(ind(1),ind(2));
fluid_out = fluid.state(ind(1),ind(2)+1);
stage     = fluid.stage(ind(1),ind(2));

% Create enthalpy and pressure arrays
if any(strcmp(stage.type,{'hex','hex_reject','regen','split','mixing','separate'}))
    type = 1;
    
elseif any(strcmp(stage.type,{'comp','exp'}))
    type = 3;
    
else
    error('not implemented');
end

switch type
    case 1 % nearly-isobaric line
        p_vect = logspace(log10(fluid_in.p),log10(fluid_out.p),num);
        h_vect = linspace(fluid_in.h,fluid_out.h,num);
        
    case 2 % polytropic line
        phi   = log(fluid_out.p/fluid_in.p)/log(fluid_out.h/fluid_in.h);
        h_vect = linspace(fluid_in.h,fluid_out.h,num);
        p_vect = fluid_in.p*(h_vect./fluid_in.h).^phi;
        
    case 3 % straight line (on T-s diagram)
        T_vect = linspace(fluid_in.T,fluid_out.T,num);
        s_vect = linspace(fluid_in.s,fluid_out.s,num);
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
        
    case 3
        Tcrit = CP1(0,0,0,'Tcrit',fluid.handle);

        if any(T_vect>0.99*Tcrit & T_vect<1.00*Tcrit)
            p1 = CP1('SmassT_INPUTS',s_vect,T_vect*1.01,'P',fluid.handle);
            p2 = CP1('SmassT_INPUTS',s_vect,T_vect*0.99,'P',fluid.handle);
            p_vect = 0.5*(p1 + p2);
        else
            p_vect = CP1('SmassT_INPUTS',s_vect,T_vect,'P',fluid.handle);
        end
        
        a0 = (p_vect == 0); %logical array of zero elements
        if any(a0)
            p_vect(a0) = NaN;
            p_vect = fillmissing(p_vect,'pchip');
        end

    otherwise
end


% for i0 = 1:num
%     fprintf(1,'i0 = %3d, h = %5.1f J/kg/K, p = %6.1f bar\n',i0,h_vect(i0),p_vect(i0)/1e5);
% end





% Write the plotting vectors
for i=1:num
    fprintf(fileID,' %12.1f %12.1f %12.1f %12.1f  %%%10s\n',T_vect(i), (T_vect(i) -273.15), p_vect(i)/1e5, s_vect(i), stage.type);
end

end