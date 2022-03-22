% This script plots out profiles along the length of the storage 
% It only works if there are two phases - i.e. one charge and one discharge
% Otherwise there is too much data and the user should really decide what
% information should be plotted.

if Nload > 2
   fprintf('There are more than two charge/discharge phases. Profiles will not be plotted.\n'); 
   
else
    
    % PCM melted fraction profiles
    figure(1)
    plot(XPsave(:,:,1))
    xlabel('Gridstep along PCM length')
    ylabel('PCM melted fraction')
    title('PCM melted fraction during phase 1 (charge)')
    
    figure(2)
    plot(XPsave(:,:,2))
    xlabel('Gridstep along PCM length')
    ylabel('PCM melted fraction')
    title('PCM melted fraction during phase 2 (discharge)')
    
    % PCM temperature profiles
    figure(3)
    plot(TPsave(:,:,1)-273.15)
    xlabel('Gridstep along PCM length')
    ylabel('PCM temperature $$^\circ$$C')
    title('PCM temperature during phase 1 (charge)')
    
    figure(4)
    plot(TPsave(:,:,2)-273.15)
    xlabel('Gridstep along PCM length')
    ylabel('PCM temperature $$^\circ$$C')
    title('PCM temperature during phase 2 (discharge)')
    
    % Steam quality profiles
    figure(5)
    plot(XSsave(:,:,1))
    xlabel('Gridstep along PCM length')
    ylabel('Steam quality')
    title('Steam quality during phase 1 (charge)')
    
    figure(6)
    plot(XSsave(:,:,2))
    xlabel('Gridstep along PCM length')
    ylabel('Steam quality')
    title('Steam quality during phase 2 (discharge)')
    
    % Steam temperature profiles
    figure(7)
    plot(TSsave(:,:,1)-273.15)
    xlabel('Gridstep along PCM length')
    ylabel('Steam temperature $$^\circ$$C')
    title('Steam temperature during phase 1 (charge)')
    
    figure(8)
    plot(TSsave(:,:,2)-273.15)
    xlabel('Gridstep along PCM length')
    ylabel('Steam temperature $$^\circ$$C')
    title('Steam temperature during phase 2 (discharge)')
    
    % Steam enthalpy profiles
    figure(9)
    plot(HSsave(:,:,1)/1e3)
    xlabel('Gridstep along PCM length')
    ylabel('Steam enthalpy, kJ/kg')
    title('Steam enthalpy during phase 1 (charge)')
    
    figure(10)
    plot(HSsave(:,:,2)/1e3)
    xlabel('Gridstep along PCM length')
    ylabel('Steam enthalpy, kJ/kg')
    title('Steam enthalpy during phase 2 (discharge)')
   
   
end