function [fluidhybrid] = update_hybrid(fluidhybrid,ind,mode)
% Update state of a fluid.

% Obtain linear index of interest and extract fluid.state
idx   = sub2ind(size(fluidhybrid.statehybrid),ind(1),ind(2));
statehybrid = fluidhybrid.statehybrid(idx);

% Update state
statehybrid = update_state_hybrid(statehybrid,fluidhybrid,mode);

% Save back into fluid.state arrays
fluidhybrid.statehybrid(idx) = statehybrid;

end