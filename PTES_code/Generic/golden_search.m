function [ x, y, iter ] = golden_search( f, a, b, tol, MinOrMax )
%GOLDEN SEARCH SOLVER
%   Uses the Golden Ratio to find the maximum or minimum of function f
%   between limits a and b. tol is the tolerance, x the result and y=f(x0).
%   iter returns the number of iterations.

fprintf('\nInside GOLDEN SEARCH subroutine!!!\n')

% inputs
% clear;
% a   = -4;
% b   = 8;
% tol = 0.01;
% MinOrMax = 'Max';
% f = @(x) -x.^2 + 5*x + 50;
% plot(linspace(a,b,1e3),f(linspace(a,b,1e3)))

make_plot = 1;
if make_plot
    a0 = a;
    b0 = b;
end


% function
if strcmp(MinOrMax,'Max')
    sign = 1;
elseif strcmp (MinOrMax,'Min')
    sign = -1;
else
    error('Error in MinOrMax selection in Golden Search function')
end

GR  = (sqrt(5)-1)/2; %Golden Ratio
D   = GR*(b - a);
x2  = a + D;
x1  = b - D;

fx1 = f(x1);
fx2 = f(x2);

if make_plot
    xv = zeros(1,100);
    fv = zeros(1,100);
    xv(1) = x1;  xv(2) = x2;
    fv(1) = fx1; fv(2) = fx2;
end

iter = 0;
err  = D;
while err>tol
    iter = iter + 1;
    
    if fx1*sign > fx2*sign
        x   = x1;
        y   = fx1;
        b   = x2;
        x2  = x1;
        fx2 = fx1;
        D   = GR*(b - a);
        x1  = b - D;
        fx1 = f(x1);
        if make_plot
            xv(2 + iter) = x1;
            fv(2 + iter) = fx1;
        end
    else
        x   = x2;
        y   = fx2;
        a   = x1;
        x1  = x2;
        fx1 = fx2;
        D   = GR*(b - a);
        x2  = a + D;
        fx2 = f(x2);
        if make_plot
            xv(2 + iter) = x2;
            fv(2 + iter) = fx2;
        end
    end
    err = D;    
end

if make_plot
    figure(4)
    plot(xv(1:(2+iter)),1-fv(1:(2+iter)),'^k')
    title('Golden search')
    xlabel('x variable')
    ylabel('function')
    xlim([a0 b0])
end

end

