function [ y1, y2, y3, y4, y5 ] = rtab_1D_5out( xv, y1v, y2v, y3v, y4v, y5v, x, mode)
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
y1  = zeros(size(x));
y2  = zeros(size(x));
y3  = zeros(size(x));
y4  = zeros(size(x));
y5  = zeros(size(x));

% Deal with boundary and outside-of-boundary points
cond1 = x == xv(i1);
cond2 = x == xv(i2);
cond3 = (x < xv(i1)) | (x > xv(i2));

y1(cond1) = y1v(i1);
y1(cond2) = y1v(i2);
y1(cond3) = NaN;

y2(cond1) = y2v(i1);
y2(cond2) = y2v(i2);
y2(cond3) = NaN;

y3(cond1) = y3v(i1);
y3(cond2) = y3v(i2);
y3(cond3) = NaN;

y4(cond1) = y4v(i1);
y4(cond2) = y4v(i2);
y4(cond3) = NaN;

y5(cond1) = y5v(i1);
y5(cond2) = y5v(i2);
y5(cond3) = NaN;

% Proceed for inside-of-boundary query points
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
Dx      = x_high-x_low;

y1_low  = y1v(i_low);
y1_high = y1v(i_high);
Dy1     = y1_high-y1_low;

y2_low  = y2v(i_low);
y2_high = y2v(i_high);
Dy2     = y2_high-y2_low;

y3_low  = y3v(i_low);
y3_high = y3v(i_high);
Dy3     = y3_high-y3_low;

y4_low  = y4v(i_low);
y4_high = y4v(i_high);
Dy4     = y4_high-y4_low;

y5_low  = y5v(i_low);
y5_high = y5v(i_high);
Dy5     = y5_high-y5_low;

% Interpolate
y1(ib) = y1_low + (x(ib)-x_low).*Dy1./Dx;
y2(ib) = y2_low + (x(ib)-x_low).*Dy2./Dx;
y3(ib) = y3_low + (x(ib)-x_low).*Dy3./Dx;
y4(ib) = y4_low + (x(ib)-x_low).*Dy4./Dx;
y5(ib) = y5_low + (x(ib)-x_low).*Dy5./Dx;

end