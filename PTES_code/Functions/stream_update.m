function [S] = stream_update(fluid,S,mode)
% STREAM_UPDATE takes two arrays of pressure and enthalpy points and
% obtains the corresponding arrays of T, k, mu and Pr

if mode == 1 % Enthalpy and pressure
    
    if strcmp(fluid.read,'CP') %read from CoolProp
        
        [S.T,S.k,S.mu,S.Pr,S.rho] = CP5('HmassP_INPUTS',S.h,S.p,'T','CONDUCTIVITY','VISCOSITY','PRANDTL','DMASS',fluid.handle);
        
    elseif strcmp(fluid.read,'TAB') %read from table
        
        hv  = fluid.TAB(:,2);
        Tv  = fluid.TAB(:,1);
        kv  = fluid.TAB(:,6);
        muv = fluid.TAB(:,7);
        Prv = fluid.TAB(:,8);
        vv  = fluid.TAB(:,3);
        [S.T, S.k, S.mu, S.Pr, v] = rtab_1D_5out( hv, Tv, kv, muv, Prv, vv, S.h, 1);
        S.rho = 1./v;
                
    end
    
    S.Cp = S.Pr.*S.k./S.mu;
    
else
    error('not implemented')
end


end

