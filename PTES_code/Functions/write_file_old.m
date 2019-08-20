function [ ] = write_file( fluid, ind, fileID, num )

% Import fluid.state and fluid.stage
fluid_in  = fluid.state(ind(1),ind(2));
fluid_out = fluid.state(ind(1),ind(2)+1);
stage     = fluid.stage(ind(1),ind(2));

% Create temperature and pressure arrays
T_vect = linspace(fluid_in.T,fluid_out.T,num);
if any(strcmp(stage.type,{'hex','hex_reject','regen','split','mixing','separate'}))
    p_vect = linspace(fluid_in.p,fluid_out.p,num);
elseif any(strcmp(stage.type,{'comp','exp'}))
    phi   = log(fluid_out.p/fluid_in.p)/log(fluid_out.T/fluid_in.T);
    p_vect = fluid_in.p*(T_vect./fluid_in.T).^phi;
else
    error('not implemented');
end

% Obtain critical pressure and temperature
Tcrit = CP1(0,0,0,'Tcrit',fluid.handle);
Pcrit = CP1(0,0,0,'Pcrit',fluid.handle);

% Compute entropy array
s_vect = zeros(1,num);
for i0 = 1:num
    %fprintf(1,'i0 = %3d, T = %5.1f K, p = %6.1f bar\n',i0,T_vect(i0),p_vect(i0)/1e5);
    if (p_vect(i0)>0.99*Pcrit && p_vect(i0)<1.00*Pcrit)
        s1 = CP1('PT_INPUTS',p_vect(i0)*1.01,T_vect(i0),'S',fluid.handle);
        s2 = CP1('PT_INPUTS',p_vect(i0)*0.99,T_vect(i0),'S',fluid.handle);
        s_vect(i0) = 0.5*(s1 + s2);
    elseif (T_vect(i0)>0.99*Tcrit && T_vect(i0)<1.01*Tcrit)
        s1 = CP1('PT_INPUTS',p_vect(i0),T_vect(i0)*1.02,'S',fluid.handle);
        s2 = CP1('PT_INPUTS',p_vect(i0),T_vect(i0)*0.98,'S',fluid.handle);
        s_vect(i0) = 0.5*(s1 + s2);
    else
        s_vect(i0) = CP1('PT_INPUTS',p_vect(i0),T_vect(i0),'S',fluid.handle);
    end
end

% Write the plotting vectors
for i=1:num
    fprintf(fileID,'%f %f %f %f\n',T_vect(i), (T_vect(i) -273), p_vect(i), s_vect(i));
end

end

