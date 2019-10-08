function [T_vect_c,sL,sG,sQ] = plot_sat_curve(fignum,n_points,n_curves,fluid_handle)

% Plot the saturation curve of a fluid using the CoolProp library.
% fignum   = figure handle for plotting
% n_points = number of plotting points for saturation curve
% n_curves = number of curves with different qualities

figure(fignum);

Tmin  = CP1(0,0,0,'Tmin',fluid_handle) + 1;
Tcrit = CP1(0,0,0,'Tcrit',fluid_handle);
T_vect_c = linspace(Tmin,Tcrit, n_points); % from freezing up to critical temperature
Qv = linspace(0.0,1.0,n_curves);
sQ = zeros(n_curves,n_points);
sL = CP1('QT_INPUTS',0.0*ones(size(T_vect_c)),T_vect_c,'S',fluid_handle);
sG = CP1('QT_INPUTS',1.0*ones(size(T_vect_c)),T_vect_c,'S',fluid_handle);
for i01=1:n_points
    for i02=1:n_curves
        sQ(i02,i01) = Qv(i02)*sG(i01) + (1-Qv(i02))*sL(i01);
    end
end



end

