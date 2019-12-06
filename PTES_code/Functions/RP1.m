function [output1] = RP1(input_pair,input1,input2,out1,fluid)

if strcmp(fluid.read,'CP')
    output1 = CP1(input_pair,input1,input2,out1,fluid.handle);
    
elseif strcmp(fluid.read,'TAB')
    
    switch out1
        case 'T'
            v  = fluid.TAB(:,1); %temperature
        case 'H'
            v  = fluid.TAB(:,2); %enthalpy
        otherwise
            error('not implemented')
    end
    
    switch input_pair
        case 'HmassP_INPUTS'
            % Tabulated data is pressure independent at the moment. Only
            % input1 (enthalpy) is used. Since enthalpy data is ordered in
            % regular intervals, the rtab1 function in mode 0 is used
            % instead of interp1 (which is slower).
            x  = fluid.TAB(:,2); %enthalpy
            output1 = rtab1(x,v,input1,0);
            
        case 'PT_INPUTS'
            % Tabulated data is pressure independent at the moment. Only
            % input2 (temperature) is used. Since temperature data is NOT ordered in regular
            % intervals, the interp1 function is used.
            x  = fluid.TAB(:,1); %temperature
            output1 = interp1(x,v,input2);
            
        otherwise
            error('not implemented')
    end
end

end