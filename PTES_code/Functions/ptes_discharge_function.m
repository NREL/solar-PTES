function [ineff] = ptes_discharge_function(gas, fluidH, fluidC, HT, CT, environ, T0, T1, pbot, PRr, PRch, mdot, Nc_dis, Ne_dis, eta, eff, ploss, t_ch, mode) %#ok<*INUSL,*INUSD>

PTES_DISCHARGE

if mode == 0 %PTES
    chi = RTeff(gas,t_ch,t_dis);
elseif mode == 2 %Heat engine only
    chi = HEeff(gas,HT.B(2).B - HT.A(1).B,t_dis);
else
    error('not implemented');
end

ineff = 1 - chi;

%fprintf(1,'\nchi = %8.3f, PRr = %8.3f\n',chi, PRr);

end