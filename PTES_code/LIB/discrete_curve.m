function [x,fx] = discrete_curve(fun0,xa,xb,xlim,n,upwards)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% function
ni    = n*10;
xr0   = (xb/xa)^(1/(n-1));
x     = zeros(ni,1);
fx    = zeros(ni,1);
x(1)  = xa; % initial point
fx(1) = fun0(x(1));
%fprintf(1,'\n');
for i=1:ni
    if upwards
        Dx = min([x(i)*(xr0-1)*(xlim-x(i))/xlim,xb-x(i)]);
    else
        Dx = max([x(i)*(xr0-1)*(xlim-x(i))/xlim,xb-x(i)]);
    end
    x(i+1)  = x(i) + Dx;
    fx(i+1) = fun0(x(i+1));
    %Df      = fx(i+1) - fx(i);
    %fprintf(1,'x=%8.3g,  fx=%8.3g,  Dx=%8.3g,  Df=%8.3g,  i=%2d\n',x(i),fx(i),Dx,Df,i);
    
    if abs(x(i+1)-xb)/xb < 1e-10
        %fprintf(1,'x=%8.3g,  fx=%8.3g, hurray!\n',x(i+1),fx(i+1));
        break
    end
end
if abs(x(i+1)-xb)/xb < 1e-10
    i=i+1;
end

if upwards
    x  = x(1:i);
    fx = fx(1:i);
else
    x  = x(i:-1:1);
    fx = fx(i:-1:1);
end

end

