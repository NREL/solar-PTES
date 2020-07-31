function [T_vect_c,sL,sG,sQ] = plot_sat_curve(fignum,n_points,n_curves,fluid)

% Plot the saturation curve of a fluid using the CoolProp library.
% fignum   = figure handle for plotting
% n_points = number of plotting points for saturation curve
% n_curves = number of curves with different qualities

figure(fignum);

Tmin  = RPN(0,0,0,'Tmin',fluid) + 1;
Tcrit = RPN(0,0,0,'Tcrit',fluid);
T_vect_c = linspace(Tmin,Tcrit, n_points); % from freezing up to critical temperature
Qv = linspace(0.0,1.0,n_curves);
sQ = zeros(n_curves,n_points);
sL = RPN('QT_INPUTS',0.0*ones(size(T_vect_c)),T_vect_c,'S',fluid);
sG = RPN('QT_INPUTS',1.0*ones(size(T_vect_c)),T_vect_c,'S',fluid);
for i01=1:n_points
    for i02=1:n_curves
        sQ(i02,i01) = Qv(i02)*sG(i01) + (1-Qv(i02))*sL(i01);
    end
end



end

