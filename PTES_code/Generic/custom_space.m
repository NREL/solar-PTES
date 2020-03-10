function x = custom_space(x1,x2,n,par)
%CUSTOM_SPACE Createe an array of points with customised spacing.
%
%   USAGE:
%   custom_space(x1,x2,n,par)
%
%   CUSTOM_SPACE Employs Matlab's in-built LINSPACE function to generate a
%   base array which is then modified according to 'par'.
%
%   'x1' and 'x2' are the starting and ending points of the array. 'n' is
%   the number of points.
%
%   'par' determines the spacing, which can be 'lin' (linear), 'log'
%   (logarithmic), or a customised positive scalar. A positive scalar
%   concentrates the points towards x1. A negative scalar concentrates the
%   points towards x2.
%
%   'par'=0 produces same results as 'par'='lin'.
%   
%   E.g.:
%   custom_space(-15,25.5,10,'lin')
%   custom_space(-15,25.5,10,1.0)
%   custom_space(1,2000,20,'log')
%   custom_space(1,2,11,0.1)
%   custom_space(1,2,11,10)

switch par
    case {'lin','log'}
    otherwise
        if ~isscalar(par)
            error(['Input "par" must be either "lin", "log" ',...
                'or a numeric scalar.'])
        end
end

switch par
    case 'lin'
        x = linspace(x1,x2,n);
        
    case 'log'
        x1 = log(x1);
        x2 = log(x2);
        x  = linspace(x1,x2,n);
        x  = exp(x);
        
    otherwise
        b  = abs(par)+1;
        x0 = linspace(0,1,n);
        x0 = x0.^(b);
        if par < 0
            x0 = 1-x0;
            x0 = x0(end:-1:1);
        end
        x  = x1 + x0*(x2-x1);
end

end

