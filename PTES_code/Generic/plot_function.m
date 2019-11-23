function [] = plot_function(fun,xmin,xmax,nx,fignum,varargin)
% Plots a one-variable function between xmin and xmax, using nx points.

if ~any(nargin == [5,6])
    error('Error in number of inputs')
end

figure(fignum)

if nargin == 5
    
    x = linspace(xmin,xmax,nx);
    y = zeros(size(x));
    for i=1:nx
        y(i) = fun(x(i));
    end
    plot(x,y)
    
elseif nargin == 6
    
    config = varargin{1};
    
    if any(strcmp(config,{'semilogx','loglog'}))
        x = logspace(log10(xmin),log10(xmax),nx);
    else
        x = linspace(xmin,xmax,nx);
    end
    y = zeros(size(x));
    for i=1:nx
        y(i) = fun(x(i));
    end
    
    switch config
        case 'semilogx'
            semilogx(x,y)
            
        case 'semilogy'
            semilogy(x,y)
            
        case 'loglog'
            loglog(x,y)
            
        otherwise
            warning('error in plot mode, using default')
            plot(x,y)
    end
    
    keyboard
end


end

