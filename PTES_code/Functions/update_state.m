function [state] = update_state (state,handle,read,TAB,mode)

switch read
    case 'CP'
        if mode == 1 %Temperature and pressure are known
                [state.rho,state.h,state.s,state.Q,~] = CP5('PT_INPUTS',state.p,state.T,'D','H','S','Q','T',handle);
        elseif mode == 2 %Enthalpy and pressure are known
                [state.T,state.rho,state.s,state.Q,~] = CP5('HmassP_INPUTS',state.h,state.p,'T','D','S','Q','P',handle);
        elseif mode == 3 % Pressure and Q are known
                [state.T,state.rho,state.s,state.h,~] = CP5('PQ_INPUTS',state.p,state.Q,'T','D','S','H','P',handle);                
        elseif mode == 4 % Temperature and Q are known
                [state.p,state.rho,state.s,state.h,~] = CP5('QT_INPUTS',state.Q,state.T,'P','D','S','H','T',handle);  
        else
            error('not implemented')
        end
        
    case 'TAB'        
        if mode == 1 %Temperature and pressure are known            
            state.h   = rtab1(TAB(:,1),TAB(:,2),state.T,0);
            state.rho = 1./rtab1(TAB(:,1),TAB(:,3),state.T,0);
            state.s   = rtab1(TAB(:,1),TAB(:,4),state.T,0);
        elseif mode == 2 %Enthalpy and pressure are known
            state.T   = rtab1(TAB(:,2),TAB(:,1),state.h,1);
            state.rho = 1./rtab1(TAB(:,2),TAB(:,3),state.h,1);
            state.s   = rtab1(TAB(:,2),TAB(:,4),state.h,1);
        else
            error('not implemented')
        end
        
    otherwise
        error('not implemented')
end

end