function [ y ] = rtab_1D_1out( xv, yv, x, mode)
%  RTAB Read Table
%  Read y(x) function from table by interpolation. xv and yv are data
%  arrays. x is the query point(s). y is the interpolated answer(s).

%  xv data points NEEDS TO BE ORDERED from small to big values.
%  If mode == 0, xv data points must be given in REGULAR INTERVALS.
%  If mode == 1, xv data points do not need to be given in REGULAR
%  INTERVALS (but interpolation is slower)

% Set boundary indices and pre-allocate y array
i1 = 1;
i2 = numel(xv);
y  = zeros(size(x));

% Deal with boundary and outside-of-boundary points
y(x == xv(i1)) = yv(i1);
y(x == xv(i2)) = yv(i2);
y((x < xv(i1)) | (x > xv(i2))) = NaN;

% INSIDE query points
% Find indices of x array of inside query points
ib = (x > xv(i1)) & (x < xv(i2));
target = x(ib);
% Find i_low and i_high incides of xv array corresponding to x(ib) query
% points
if mode == 0
    fx1 = xv(i1) - target;
    fx2 = xv(i2) - target;
    i3 = (i1*fx2 - i2*fx1)./(fx2 - fx1);
    i_low  = floor(i3);
    i_high = i_low + 1;
elseif mode ==1
    i_low  = zeros(size(target));
    i_high = zeros(size(target));
    for i = 1:length(target)
        i_low(i)  = find(xv <= target(i),1,'last');
        i_high(i) = i_low(i) + 1;
    end    
else
    error('mode not implemented')
end
% Determine surrounding xv and yv values for interpolation
x_low   = xv(i_low);
x_high  = xv(i_high);
y_low  = yv(i_low);
y_high = yv(i_high);
Dx = x_high-x_low;
Dy = y_high-y_low;
% Interpolate
y(ib) = y_low + (x(ib)-x_low).*Dy./Dx;

end