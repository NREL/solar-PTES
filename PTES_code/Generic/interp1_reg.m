function vq = interp1_reg(x,v,xq,scale,method)

switch method
    case 'linear'
    otherwise
        error('not implemented')
end

switch scale
    case 'lin'
        % Obtain xmin, dx and numx
        xmin  = min(x);
        dx    = x(2) - x(1);
        numx  = length(x);
        
        % Obtain the x-indices (and values) around the xq
        % query points
        ix1 = floor(1 + (xq-xmin)./dx);
        ix1(ix1<=1) = 1;
        ix1(ix1>=numx) = numx - 1;
        ix2 = ix1 + 1;
        x1  = x(ix1); % x1 = xmin + (ix1 - 1).*dx
        x2  = x(ix2); % x2 = xmin + (ix2 - 1).*dx
        
    case 'log'
        
        % Obtain xmin, xr and numx
        xmin  = min(x);
        xr    = x(2)/x(1);
        numx  = length(x);
        
        % Obtain the y-indices (and values) around the Yq query
        % points.
        ix1 = floor(1 + log(xq/xmin)./log(xr));
        ix1(ix1<=1) = 1;
        ix1(ix1>=numx) = numx - 1;
        ix2 = ix1 + 1;
        x1  = x(ix1); %x1 = xmin*xr.^(ix1 - 1);
        x2  = x(ix2); %x2 = xmin*xr.^(ix2 - 1);
        
    otherwise
        error('not implemented')
end

% Extract coefficients for interpolation formula (i.e.
% values of v on corners around query point).
v1 = v(ix1);
v2 = v(ix2);

% Apply interpolation formula
vq  = ((x2-xq).*v1 + (xq-x1).*v2)./(x2-x1);

end