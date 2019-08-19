function [ Mdot ] = total_Mdot( fluid_streams, ind )

% Compute total mass flow rate of the fluid_streams
Mdot = 0;
n = numel(fluid_streams);
for i=1:n
    Mdot   = Mdot + fluid_streams(i).state(ind(1),ind(2)).mdot;
end

end

