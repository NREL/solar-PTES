function [input,total] = SUM_ENERGY(input,total,mode,N) 

% Modes 1 and 2 - add to the total
% Mode 1 for compressors, expanders, heat exchangers
% Mode 2 for the gas.stage class

% Modes 3 and 4 - assign to the final CCMP, EXP classes etc.
switch mode
    
    case 1
        total.w    = total.w + input.w(N) ;
        total.q    = total.q + input.q(N) ;
        total.Dh   = total.Dh + input.Dh(N) ;
        total.sirr = total.sirr + input.sirr(N) ;
    
    case 2
        for i = 1 : N
           total(i).w    = total(i).w + input(i).w ;
           total(i).q    = total(i).q + input(i).q ;
           total(i).Dh   = total(i).Dh + input(i).Dh ;
           total(i).sirr = total(i).sirr + input(i).sirr ;
        end
        
    case 3
        input.w(N) = total.w ;
        input.q(N) = total.q ;
        input.Dh(N) = total.Dh ;
        input.sirr(N) = total.sirr ;
        
    case 4
        for i = 1 : N
           input(i).w = total(i).w ;
           input(i).q = total(i).q ;
           input(i).Dh = total(i).Dh ;
           input(i).sirr = total(i).sirr ;
        end
        
end

end