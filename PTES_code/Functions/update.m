function [fluid] = update (fluid,ind,mode)

% Obtain linear index of interest and extract fluid.state
idx   = sub2ind(size(fluid.state),ind(1),ind(2));
state = fluid.state(idx);

state = update_state(state,fluid.handle,fluid.read,fluid.TAB,fluid.IDL,mode);

fluid.state(idx) = state;

end