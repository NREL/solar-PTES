function [ x, y, xv, yv, iter ] = golden_search( f, a, b, tol, MinOrMax, ne)
%GOLDEN SEARCH SOLVER
%   Uses the Golden Ratio to find the *absolute* maximum or minimum of
%   function f between limits a and b. tol is the tolerance, x the result
%   and y=f(x). iter returns the number of iterations. MinOrMax describes
%   the objective of the search. ne is the maximum number of iterations and
%   the number of elements of the arrays xv and yv, which keep track of the
%   convergence procedure.
%   
%   Example call:
%   f = @(x1) my_function(x1,x2,x3);
%   [x1_opt,y_opt,x1v,yv,iter] = golden_search(f,-20,20,1e-3,'Min',100);


if strcmp(MinOrMax,'Max')
    sign = 1;
elseif strcmp (MinOrMax,'Min')
    sign = -1;
end

GR  = (sqrt(5)-1)/2; %Golden Ratio
D   = GR*(b - a);
x2  = a + D;
x1  = b - D;

fx1 = f(x1);
fx2 = f(x2);

xv = zeros(ne,1);
yv = zeros(ne,1);
xv(1) = x1;  xv(2) = x2;
yv(1) = fx1; yv(2) = fx2;

iter = 2;
err  = D;
while err>tol
    iter = iter + 1;
    
    if iter == ne
        warning('gold_search exiting at the maximum number of iterations')
        break
    end
    
    if fx1*sign > fx2*sign
        x   = x1;
        y   = fx1;
        b   = x2;
        x2  = x1;
        fx2 = fx1;
        D   = GR*(b - a);
        x1  = b - D;
        fx1 = f(x1);
        
        xv(iter) = x1;
        yv(iter) = fx1;
    else
        x   = x2;
        y   = fx2;
        a   = x1;
        x1  = x2;
        fx1 = fx2;
        D   = GR*(b - a);
        x2  = a + D;
        fx2 = f(x2);
        
        xv(iter) = x2;
        yv(iter) = fx2;
    end
    err = D;
    
end

end

