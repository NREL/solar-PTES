function [fluid,environ,iF,iE] = hex_set(fluid,indF,environ,indE,Tset,eff,ploss)

% Import fluid.state and fluid.stage
state = fluid.state(indF(1),indF(2));
stage = fluid.stage(indF(1),indF(2));

% Extract initial conditions
T_st    = state.T;
p_st    = state.p;
h_st    = state.h;
s_st    = state.s;
T0      = environ.T0;

% Set final conditions and update
state.T = T_st + eff*(Tset-T_st);
state.p = p_st*(1-ploss);
state   = update_state(state,fluid.handle,fluid.read,fluid.TAB,1);

% Compute stage energy flows
stage.Dh  = state.h - h_st;
stage.type= 'hex_reject';
if all([T_st >= Tset, Tset >= T0]) || all([T_st <= Tset, Tset <= T0]) %cooling/heating versus environment
    %fprintf(1,'\nCooling/heating versus environment!\n')
    stage.w    = 0;
    stage.q    = stage.Dh + stage.w;
    stage.sirr = state.s - s_st - stage.Dh/T0;
elseif all([T_st <= Tset, Tset > T0]) %joule heating
    %fprintf(1,'\nJoule heating!\n')
    stage.q    = 0;
    stage.w    = stage.q - stage.Dh;
    stage.sirr = state.s - s_st;
else
    error('Selected conditions do not allow cooling/heating versus ambient or Joule heating');
end
environ.sink(indE(1),indE(2)).DHdot = - stage.q*state.mdot;

% Reassing values to fluid.state and fluid.stage
fluid.state(indF(1),indF(2)+1) = state;
fluid.stage(indF(1),indF(2))   = stage;

%Increase stage counter
iF  = indF(2) + 1;
iE  = indE(2) + 1;
end

