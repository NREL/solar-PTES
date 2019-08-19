classdef environment_class
    properties
        T0 = 0
        p0 = 0
        sink = sink_class
    end
    methods
        function obj = environment_class(T0,p0,num)
            obj.T0 = T0;
            obj.p0 = p0;
            obj.sink(1:2,1:num) = sink_class; % first row for charge, second row for discharge
        end
    end
end

