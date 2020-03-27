function [state] = update_state(state,fluid,mode)

switch mode
    case 1 %Temperature and pressure are known
        [state.rho] = RP1('PT_INPUTS',state.p,state.T,'D',fluid);
        [state.h]   = RP1('PT_INPUTS',state.p,state.T,'H',fluid);
        [state.s]   = RP1('PT_INPUTS',state.p,state.T,'S',fluid);
        [state.Q]   = RP1('PT_INPUTS',state.p,state.T,'Q',fluid);
        
    case 2 %Enthalpy and pressure are known
        [state.T]   = RP1('HmassP_INPUTS',state.h,state.p,'T',fluid);
        [state.rho] = RP1('HmassP_INPUTS',state.h,state.p,'D',fluid);
        [state.s]   = RP1('HmassP_INPUTS',state.h,state.p,'S',fluid);
        [state.Q]   = RP1('HmassP_INPUTS',state.h,state.p,'Q',fluid);
        
    case 3 % Pressure and Q are known
        [state.T]   = RP1('PQ_INPUTS',state.p,state.Q,'T',fluid);
        [state.rho] = RP1('PQ_INPUTS',state.p,state.Q,'D',fluid);
        [state.s]   = RP1('PQ_INPUTS',state.p,state.Q,'S',fluid);
        [state.h]   = RP1('PQ_INPUTS',state.p,state.Q,'H',fluid);
        
    case 4 % Temperature and Q are known
        [state.p]   = RP1('QT_INPUTS',state.Q,state.T,'P',fluid);
        [state.rho] = RP1('QT_INPUTS',state.Q,state.T,'D',fluid);
        [state.s]   = RP1('QT_INPUTS',state.Q,state.T,'S',fluid);
        [state.h]   = RP1('QT_INPUTS',state.Q,state.T,'H',fluid);
    otherwise
        error('not implemented')
end

end