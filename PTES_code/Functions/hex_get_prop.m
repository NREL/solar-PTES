function [ Tv, Cpv, hv ] = hex_get_prop( fluid, TC1, TH2, pressure, n )
% Obtain the Tv, Cpv and hv arrays of a given fluid for the hex
% subroutines

Tv  = linspace(TC1,TH2,n)'; %temperature array between TC1 and TH2

%Set Cp and hv arrays
if strcmp(fluid.read,'CP') %read from CoolProp
    
    [Cpv,hv,~,~,~] = CP5('PT_INPUTS',linspace(pressure,pressure,n)',Tv,'C','H','T','T','T',fluid.handle);
    
    %Fill in zero-gaps by interpolation (sometimes created at/around critical temperature)
    for i0=1:n
        if Cpv(i0) == 0 && hv(i0) == 0 % identify gap
            if i0 == 1 % gap at begginning of array
                i_up    = find(Cpv(find(~Cpv,1,'first')+1:end),2,'first')+find(~Cpv,1,'first'); % first two non-zero elements after i0
                Cpv(i0) = Cpv(i_up(1)) - (Cpv(i_up(2)) - Cpv(i_up(1)))/(i_up(2) - i_up(1))*(i_up(1) - i0);
                hv(i0)  = hv(i_up(1))  - (hv(i_up(2))  - hv(i_up(1)))/ (i_up(2) - i_up(1))*(i_up(1) - i0);
            elseif i0 == n || all(~Cpv(i0:n)) % gap at end of array
                i_down  = find(Cpv(1:find(~Cpv,1,'first')),2,'last'); % first two non-zero elements before i0
                Cpv(i0) = Cpv(i_down(1)) + (Cpv(i_down(2)) - Cpv(i_down(1)))/(i_down(2) - i_down(1))*(i0 - i_down(1));
                hv(i0)  = hv(i_down(1))  + (hv(i_down(2))  - hv(i_down(1)))/ (i_down(2) - i_down(1))*(i0 - i_down(1));
            else
                i_up   = find(Cpv(find(~Cpv,1,'first')+1:end),1,'first')+find(~Cpv,1,'first'); % first non-zero element after i0
                i_down = find(Cpv(1:find(~Cpv,1,'first')),1,'last'); % first non-zero element before i0
                Cpv(i0) = Cpv(i_down) + (Cpv(i_up) - Cpv(i_down))/(i_up - i_down)*(i0 - i_down);
                hv(i0)  = hv(i_down)  + (hv(i_up)  - hv(i_down)) /(i_up - i_down)*(i0 - i_down);
            end
        end
    end
    
elseif strcmp(fluid.read,'TAB') %read from table
    Tx  = fluid.A(:,1);
    Cpy = fluid.A(:,5);
    hy  = fluid.A(:,2);
    Cpv = rtab1(Tx,Cpy,Tv,0);
    hv  = rtab1(Tx,hy,Tv,0);
else
    error('not implemented')
end

end

