function [fit]=test(x)
         g=sum(x(1:2)-0.5).^2;          %DTLZ2
         f1= cos(x(1)*pi/2)*(1+g);
         f2= sin(x(1)*pi/2)*(1+g);
         fit = [f1 f2];
  
end