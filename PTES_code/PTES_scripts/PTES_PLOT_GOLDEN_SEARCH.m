figure(3)
plot(xv(1:iter).*PRch,1-yv(1:iter),'^k')
title('Golden search')
xlabel('Discharge pressure ratio')
ylabel('Efficiency')
xlim([PRr_min PRr_max].*PRch)