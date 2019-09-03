function [tank_state] = update_tank (tank_state, fluid, T0, mode)

tank_state = update_state(tank_state,fluid.handle,fluid.read,fluid.TAB,mode);

tank_state.V = tank_state.M/tank_state.rho;
tank_state.H = tank_state.M*tank_state.h;
tank_state.S = tank_state.M*tank_state.s;
tank_state.B = tank_state.H - T0*tank_state.S;

end