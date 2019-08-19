function [fluid] = update (fluid,ind,mode)

% Obtain linear index of interest and extract fluid.state
idx   = sub2ind(size(fluid.state),ind(1),ind(2));
state = fluid.state(idx);

state = update_state(state,fluid.handle,fluid.read,fluid.A,mode);

% switch fluid.read
%     case 'CP'
%         if mode == 1 %Temperature and pressure are known
%             %fprintf(1,'\np = %.4g [bar], T = %.4g [K]\n',fluid.p/1e5,fluid.T)
%             try
%                 state.rho = CoolProp.PropsSI('D','T',state.T,'P',state.p,fluid.name);
%                 state.h   = CoolProp.PropsSI('H','T',state.T,'P',state.p,fluid.name);
%                 state.s   = CoolProp.PropsSI('S','T',state.T,'P',state.p,fluid.name);
%             catch
%                 warning('***averaging properties in update subroutine***')
%                 T1 = state.T*1.01;
%                 T2 = state.T*0.99;
%                 p1 = state.p*1.01;
%                 p2 = state.p*0.99;
%                 rho11 = CoolProp.PropsSI('D','T',T1,'P',p1,fluid.name);
%                 rho12 = CoolProp.PropsSI('D','T',T1,'P',p2,fluid.name);
%                 rho21 = CoolProp.PropsSI('D','T',T2,'P',p1,fluid.name);
%                 rho22 = CoolProp.PropsSI('D','T',T2,'P',p2,fluid.name);
%                 state.rho = 0.25*(rho11+rho12+rho21+rho22);
%                 h11 = CoolProp.PropsSI('H','T',T1,'P',p1,fluid.name);
%                 h12 = CoolProp.PropsSI('H','T',T1,'P',p2,fluid.name);
%                 h21 = CoolProp.PropsSI('H','T',T2,'P',p1,fluid.name);
%                 h22 = CoolProp.PropsSI('H','T',T2,'P',p2,fluid.name);
%                 state.h = 0.25*(h11+h12+h21+h22);
%                 s11 = CoolProp.PropsSI('S','T',T1,'P',p1,fluid.name);
%                 s12 = CoolProp.PropsSI('S','T',T1,'P',p2,fluid.name);
%                 s21 = CoolProp.PropsSI('S','T',T2,'P',p1,fluid.name);
%                 s22 = CoolProp.PropsSI('S','T',T2,'P',p2,fluid.name);
%                 state.s = 0.25*(s11+s12+s21+s22);
%             end
%         elseif mode == 2 %Enthalpy and pressure are known
%             %fprintf(1,'\np = %.4g [bar], h = %.4g [J/kg]\n',fluid.p/1e5,fluid.h/1e6)
%             try
%                 state.T   = CoolProp.PropsSI('T','H',state.h,'P',state.p,fluid.name);
%                 state.rho = CoolProp.PropsSI('D','H',state.h,'P',state.p,fluid.name);
%                 state.s   = CoolProp.PropsSI('S','H',state.h,'P',state.p,fluid.name);
%             catch
%                 error('not implemented')
%             end
%         else
%             error('not implemented')
%         end
%         
%     case 'TAB'        
%         if mode == 1 %Temperature and pressure are known            
%             state.h   = rtab1(fluid.A(:,1),fluid.A(:,2),state.T,0);
%             state.rho = 1./rtab1(fluid.A(:,1),fluid.A(:,3),state.T,0);
%             state.s   = rtab1(fluid.A(:,1),fluid.A(:,4),state.T,0);
%         elseif mode == 2 %Enthalpy and pressure are known
%             state.T   = rtab1(fluid.A(:,2),fluid.A(:,1),state.T,1);
%             state.rho = 1./rtab1(fluid.A(:,2),fluid.A(:,3),state.T,1);
%             state.s   = rtab1(fluid.A(:,2),fluid.A(:,4),state.T,1);
%         else
%             error('not implemented')
%         end
%         
%     otherwise
%         error('not implemented')
% end

fluid.state(idx) = state;

end