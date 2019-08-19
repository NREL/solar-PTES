function [ pp, np ] = find_pinch_points(TvH,CpvH,TvC,CpvC,mH,mC)
% Find pinch points
% This applies the bisection rule to find the root for which mH*CpH = mC*CpC
% Valid for the case with DTp=0

% TvH and CpvH are the x and y arrays containing discrete values of CpH(T)
% TvC and CpvC are the x and y arrays containing discrete values of CpC(T)
% mH and mC are the mass flow rates of the hot and cold fluids

n = length(TvC); %this is assumed to be the same for TvC and TvH
TC1 = max([TvC(1),TvH(1)]);
TH2 = min([TvC(n),TvH(n)]);

Tv  = linspace(TC1,TH2,n)';
CpH = rtab1(TvH,CpvH,Tv(:),0);
CpC = rtab1(TvC,CpvC,Tv(:),0);

fx  = mH.*CpH-mC.*CpC;

%Visualize pinch points?
% figure(19)
% plot(Tv,fx)

fx0 = fx(1);
np  = 0; %number of pinch points found
num = length(Tv);
ip  = zeros(num,1);
for i=2:num
    if fx(i)*fx0<0
        np = np + 1;
        ip(np) = i; %pinch point index
        fx0 = fx(i);
    end
end

if np==0
    % Intermediate pinch points not found. Using minimum instead
    np = 1;
    [~,i] = min(abs(fx));
    pp = Tv(i);
    return
end

pp = zeros(np,1);

for i=1:np
    x1 = Tv(ip(i)-1); %right before pinch point
    x2 = Tv(ip(i));   %right after pinch point
    fx1 = fx(ip(i)-1);
    fx2 = fx(ip(i));
    x3  = (x1 + x2)/2;
    fx3 = rtab1(Tv,mH.*CpH,x3,0) - mC*rtab1(Tv,mC.*CpC,x3,0);
    for i2=1:10
        if (fx1*fx3 < 0)
            x2  = x3;
            fx2 = fx3;
        elseif (fx2*fx3 < 0)
            x1  = x3;
            fx1 = fx3;
        else
            error('***Root not found***')
        end
        x3  = (x1 + x2)/2;
        fx3 = rtab1(Tv,mH.*CpH,x3,0) - rtab1(Tv,mC.*CpC,x3,0);
    end
    pp(i) = x3;
end

end

