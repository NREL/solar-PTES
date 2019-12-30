function [output1] = RP1(input_pair,input1,input2,out1,fluid)

if strcmp(fluid.read,'CP')
    output1 = CP1(input_pair,input1,input2,out1,fluid.handle);
    
elseif strcmp(fluid.read,'TAB')
    
    switch out1
        case 'T'
            v  = fluid.TAB(:,1); %temperature
        case 'H'
            v  = fluid.TAB(:,2); %enthalpy
        case 'S'
            v  = fluid.TAB(:,4) ; % entropy
        otherwise
            error('not implemented')
    end
    
    switch input_pair
        case 'HmassP_INPUTS'
            % Tabulated data is pressure independent at the moment. Only
            % input1 (enthalpy) is used.
            x  = fluid.TAB(:,2); %enthalpy
            output1 = interp1(x,v,input1);
            
        case 'PT_INPUTS'
            % Tabulated data is pressure independent at the moment. Only
            % input2 (temperature) is used.
            x  = fluid.TAB(:,1); %temperature
            output1 = interp1(x,v,input2);
            
        otherwise
            error('not implemented')
    end

elseif strcmp(fluid.read,'IDL')
    % Calculate properties of an ideal gas
    
    switch input_pair
        case 'HmassP_INPUTS'
            switch out1
                case 'T'
                    output1 = input1 / fluid.IDL.cp ;
                case 'S'
                    T       = input1 / fluid.IDL.cp ; % Find T first, then S
                    output1 = fluid.IDL.cp * log(T) - fluid.IDL.R * log(input2) - fluid.IDL.s0 ;
                case 'D'
                    T       = input1 / fluid.IDL.cp ; % Find T first, then density
                    output1 = input2 / (fluid.IDL.R * T) ;
            end
            
        case 'PT_INPUTS'
            switch out1
                case 'H'
                    output1 = fluid.IDL.cp  * input2 - fluid.IDL.h0 ;
                case 'S'
                    output1 = fluid.IDL.cp * log(input2) - fluid.IDL.R * log(input1) - fluid.IDL.s0 ;
            end
        case 'PSmass_INPUTS'
            switch out1
                case 'T'
                    output1 = exp( (input2 + fluid.IDL.R * input1)/fluid.IDL.cp) ;
                case 'H'
                    T       = exp( (input2 + fluid.IDL.R * input1)/fluid.IDL.cp) ;
                    output1 = fluid.IDL.cp * T - h0 ;
            end
    end
    
end

end