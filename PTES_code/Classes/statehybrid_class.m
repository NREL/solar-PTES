classdef statehybrid_class
    properties
        T    = 0 % temperature, K
        p    = 0 % pressure, Pa
        mdot = 0 % mass flow rate, kg/s
        rho  = 0 % density, kg/m3
        h = 0 % specific enthalpy, J/kg
        s = 0 % specific entropy, J/kg/K
        Q = 0 % vapour quality
        Cp = 0 %specific heat
        mu = 0 %dynamic viscosity
        k_therm = 0 %thermal conductivity
    end    
end
