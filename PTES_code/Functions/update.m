function [fluid] = update(fluid,ind,mode)
% Update state of a fluid.

% Obtain linear index of interest and extract fluid.state
idx   = sub2ind(size(fluid.state),ind(1),ind(2));
state = fluid.state(idx);

% Update state
state = update_state(state,fluid,mode);

% Save back into fluid.state arrays
fluid.state(idx) = state;

end