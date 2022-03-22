% This script estimates the cost of the steamPCM system and also prints out
% geometric information


Vsteel = Npipe * Lp * pi * dp * 0.5e-3 ;
Vpcm   = Lp * pi * (do^2 - dp^2) / 4 ;

fprintf('\n---------------------------------\n');
fprintf('GEOMETRIC AND ECONOMIC RESULTS:\n');
fprintf('---------------------------------\n');
fprintf('\nGeometry:\n')

fprintf('Pipe length:                 %8.2f m\n',Lp)
fprintf('Internal pipe diameter:      %8.2f m\n',dp)
fprintf('PCM thickness:               %8.2f m\n',0.5*(do-dp))
fprintf('Total tube diameter:         %8.2f m\n',do)
fprintf('Number of tubes per tank:    %8.2f \n',Npipe)
fprintf('Number of tanks:             %8.2f \n',Ntank)
fprintf('PCM volume per tank:         %8.2f m3\n',Vpcm)
fprintf('Steel volume per tank:       %8.2f m3\n',Vsteel)
fprintf('Total PCM volume:            %8.2f m3\n',Vpcm*Ntank)
fprintf('Total steel volume:          %8.2f m3\n',Vsteel*Ntank)

fprintf('\nEconomics:\n')
fprintf('PCM cost:                    %8.2f $\n',Vpcm*Ntank*Cpcm)
fprintf('Steel cost:                  %8.2f $\n',Vsteel*Ntank*Csteel)
fprintf('Total cost:                  %8.2f $\n',Vpcm*Ntank*Cpcm + Vsteel*Ntank*Csteel)

