function [state] = update_state_hybrid(state,fluid,mode)

switch mode
    case 1 %Temperature and pressure are known
        [state.rho,state.h,state.s,state.Q,state.Cp,state.mu,state.k_therm] =...
            RPN('PT_INPUTS',state.p,state.T,{'D','H','S','Q','CPMASS','V','L'},fluid);
        
    case 2 %Enthalpy and pressure are known
        [state.T,state.rho,state.s,state.Q] =...
            RPN('HmassP_INPUTS',state.h,state.p,{'T','D','S','Q'},fluid);
        
    case 3 % Pressure and Q are known
        [state.T,state.rho,state.s,state.h] =...
            RPN('PQ_INPUTS',state.p,state.Q,{'T','D','S','H'},fluid);
        
    case 4 % Temperature and Q are known
        [state.p,state.rho,state.s,state.h] =...
            RPN('QT_INPUTS',state.Q,state.T,{'P','D','S','H'},fluid);
    
    otherwise
        error('not implemented')
end

end