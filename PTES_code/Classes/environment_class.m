classdef environment_class
    properties
        T0 = 0
        p0 = 0
        sink = sink_class
    end
    methods
        function obj = environment_class(T0,p0,numPeriods,numStates)
            obj.T0 = T0;
            obj.p0 = p0;
            obj.sink(1:numPeriods,1:numStates) = sink_class; % first row for charge, second row for discharge
        end
    end
end

