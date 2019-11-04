function [] = plot_function(fun,xmin,xmax,nx,fignum)
% Plots a one-variable function between xmin and xmax, using nx points.

figure(fignum)
x = linspace(xmin,xmax,nx);
y = zeros(size(x));
for i=1:nx
    y(i) = fun(x(i));
end
plot(x,y)

end

