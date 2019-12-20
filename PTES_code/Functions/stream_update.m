function [S] = stream_update(fluid,S,mode)
% STREAM_UPDATE takes two arrays of pressure and enthalpy points and
% obtains the corresponding arrays of T, k, mu and Pr

if mode == 1 % Enthalpy and pressure
    
    if strcmp(fluid.read,'CP') %read from CoolProp
        
        [S.T,S.k,S.mu,S.Pr,S.rho] = CP5('HmassP_INPUTS',S.h,S.p,'T','CONDUCTIVITY','VISCOSITY','PRANDTL','DMASS',fluid.handle);
        S.v  = 1./S.rho;
        %S.Pr = S.Cp.*S.mu./S.k;
        
    elseif strcmp(fluid.read,'TAB') %read from table
        
        hv  = fluid.TAB(:,2);
        Tv  = fluid.TAB(:,1);
        kv  = fluid.TAB(:,6);
        muv = fluid.TAB(:,7);
        Prv = fluid.TAB(:,8);
        vv  = fluid.TAB(:,3);
        [S.T, S.k, S.mu, S.Pr, S.v] = rtab_1D_5out( hv, Tv, kv, muv, Prv, vv, S.h, 1);
        S.rho = 1./S.v;
                
    end
    
    S.Cp  = S.Pr.*S.k./S.mu;
    
    Cp_mean = median(S.Cp);
    Pr_mean = median(S.Pr);
    cond1 = S.Cp < 0 | S.Cp > 100*Cp_mean | S.Cp < Cp_mean/100;
    cond2 = S.Pr < 0 | S.Pr > 100*Pr_mean | S.Pr < Pr_mean/100;
     
    if any([cond1;cond2])
%         figure(31)
%         plot(S.T)
%         figure(32)
%         semilogy(S.Pr)
%         figure(33)
%         semilogy(S.Cp)
        
        Cp_mean = median(S.Cp(~cond1));
        S.Cp(cond1) = Cp_mean;
        Pr_mean = median(S.Pr(~cond2));
        S.Pr(cond2) = Pr_mean;
        
%         figure(34)
%         semilogy(S.Pr)
%         figure(35)
%         semilogy(S.Cp)
%         keyboard
    end
    
else
    error('not implemented')
end


end

