classdef double_tank_class
    properties
        A = tank_state_class; % source tank (full at start of charge cycle)
        B = tank_state_class; % sink tank (full at start of discharge cycle)
    end
    methods
         function obj = double_tank_class(num)
            obj.A(1:num) = tank_state_class;
            obj.B(1:num) = tank_state_class;
         end
         function obj = reset_tanks(obj,T,p)
            obj.A(:) = tank_state_class;
            obj.B(:) = tank_state_class;
            obj.A(1).T = T;
            obj.A(1).p = p;
        end
    end
end