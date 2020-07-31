classdef tank_state_class
    properties
        % Intensive/specific properties
        T = 0 % temperature, K
        p = 0 % pressure, Pa
        rho = 0 % density, kg/m3
        h = 0 % specific enthalpy, J/kg
        s = 0 % specific entropy, J/kg/K
        Q = 0 % vapour quality
        
        % Extensive properties
        M = 0 % mass, kg
        V = 0 % volume, m3
        H = 0 % enthalpy, J
        S = 0 % entropy, J/K
        B = 0 % exergy, J
    end
end