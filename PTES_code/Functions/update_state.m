function [state] = update_state (state,handle,read,A,mode)

switch read
    case 'CP'
        if mode == 1 %Temperature and pressure are known
            %fprintf(1,'\np = %.4g [bar], T = %.4g [K]\n',fluid.p/1e5,fluid.T)
%             try
                %state.rho = CP1('PT_INPUTS',state.p,state.T,'D',handle);
                %state.h   = CP1('PT_INPUTS',state.p,state.T,'H',handle);
                %state.s   = CP1('PT_INPUTS',state.p,state.T,'S',handle);
                %state.Q   = CP1('PT_INPUTS',state.p,state.T,'Q',handle);
                [state.rho,state.h,state.s,state.Q,~] = CP5('PT_INPUTS',state.p,state.T,'D','H','S','Q','T',handle);
%             catch
%                 warning('***averaging properties in update subroutine***')
%                 T1 = state.T*1.01;
%                 T2 = state.T*0.99;
%                 p1 = state.p*1.01;
%                 p2 = state.p*0.99;
%                 rho11 = CP1('PT_INPUTS',p1,T1,'D',handle);
%                 rho12 = CP1('PT_INPUTS',p1,T2,'D',handle);
%                 rho21 = CP1('PT_INPUTS',p2,T1,'D',handle);
%                 rho22 = CP1('PT_INPUTS',p2,T2,'D',handle);
%                 state.rho = 0.25*(rho11+rho12+rho21+rho22);
%                 h11   = CP1('PT_INPUTS',p1,T1,'H',handle);
%                 h12   = CP1('PT_INPUTS',p1,T2,'H',handle);
%                 h21   = CP1('PT_INPUTS',p2,T1,'H',handle);
%                 h22   = CP1('PT_INPUTS',p2,T2,'H',handle);
%                 state.h = 0.25*(h11+h12+h21+h22);
%                 s11   = CP1('PT_INPUTS',p1,T1,'S',handle);
%                 s12   = CP1('PT_INPUTS',p1,T2,'S',handle);
%                 s21   = CP1('PT_INPUTS',p2,T1,'S',handle);
%                 s22   = CP1('PT_INPUTS',p2,T2,'S',handle);
%                 state.s = 0.25*(s11+s12+s21+s22);
%             end
        elseif mode == 2 %Enthalpy and pressure are known
            %fprintf(1,'\np = %.4g [bar], h = %.4g [J/kg]\n',fluid.p/1e5,fluid.h/1e6)
%             try
                %state.T   = CP1('HmassP_INPUTS',state.h,state.p,'T',handle);
                %state.rho = CP1('HmassP_INPUTS',state.h,state.p,'D',handle);
                %state.s   = CP1('HmassP_INPUTS',state.h,state.p,'S',handle);
                %state.Q   = CP1('HmassP_INPUTS',state.h,state.p,'Q',handle);
                [state.T,state.rho,state.s,state.Q,~] = CP5('HmassP_INPUTS',state.h,state.p,'T','D','S','Q','P',handle);
%             catch
%                 warning('***averaging properties in update subroutine***')
%                 h1 = state.h*1.01;
%                 h2 = state.h*0.99;
%                 p1 = state.p*1.01;
%                 p2 = state.p*0.99;
%                 T11   = CP1('HmassP_INPUTS',h1,p1,'T',handle);
%                 T12   = CP1('HmassP_INPUTS',h1,p2,'T',handle);
%                 T21   = CP1('HmassP_INPUTS',h2,p1,'T',handle);
%                 T22   = CP1('HmassP_INPUTS',h2,p2,'T',handle);
%                 state.T = 0.25*(T11+T12+T21+T22);
%                 rho11 = CP1('HmassP_INPUTS',h1,p1,'D',handle);
%                 rho12 = CP1('HmassP_INPUTS',h1,p2,'D',handle);
%                 rho21 = CP1('HmassP_INPUTS',h2,p1,'D',handle);
%                 rho22 = CP1('HmassP_INPUTS',h2,p2,'D',handle);
%                 state.rho = 0.25*(rho11+rho12+rho21+rho22);
%                 s11   = CP1('HmassP_INPUTS',h1,p1,'S',handle);
%                 s12   = CP1('HmassP_INPUTS',h1,p2,'S',handle);
%                 s21   = CP1('HmassP_INPUTS',h2,p1,'S',handle);
%                 s22   = CP1('HmassP_INPUTS',h2,p2,'S',handle);
%                 state.s = 0.25*(s11+s12+s21+s22);
%             end
        else
            error('not implemented')
        end
        
    case 'TAB'        
        if mode == 1 %Temperature and pressure are known            
            state.h   = rtab1(A(:,1),A(:,2),state.T,0);
            state.rho = 1./rtab1(A(:,1),A(:,3),state.T,0);
            state.s   = rtab1(A(:,1),A(:,4),state.T,0);
        elseif mode == 2 %Enthalpy and pressure are known
            state.T   = rtab1(A(:,2),A(:,1),state.h,1);
            state.rho = 1./rtab1(A(:,2),A(:,3),state.h,1);
            state.s   = rtab1(A(:,2),A(:,4),state.h,1);
        else
            error('not implemented')
        end
        
    otherwise
        error('not implemented')
end

end