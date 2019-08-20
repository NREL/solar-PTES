function [output1] = CP1(input_pair,input1,input2,out1,handle)
% Wrapper for the low-level CoolProp function "AbstractState_keyed_output"

% Usage example (trivial answer):
% Tcrit = CP1(0,0,0,'Tcrit',gas.handle)

% Usage example (non-trivial answer):
% rho = CP1('PT_INPUTS',1e5,300,'D',gas.handle)
% rho = CP1('PT_INPUTS',[1e5,1e5,1e5],[300,400,500],'D',gas.handle)
% For enthalpy and pressure use: 'HmassP_INPUTS'
% For pressure and entropy use:  'PSmass_INPUTS'
% Others: 'QT_INPUTS', 'PQ_INPUTS'...
% or visit documentation in CoolProp Doxygen:
% http://www.coolprop.org/_static/doxygen/html/namespace_cool_prop.html#a58e7d98861406dedb48e07f551a61efb

% Declaring error variables
ierr = 0; buffer_size = 10;
herr= char((1:1:buffer_size));

% Obtaining index of parameter to compute
output(1,1) = calllib('coolprop','get_param_index',out1);

if input_pair==0 % Trivial answer
    output1 = calllib('coolprop','AbstractState_keyed_output',handle,output,ierr,herr,buffer_size);
    
else % Non-trivial answer
    % Obtain input pair index
    INPUTS = calllib('coolprop','get_input_pair_index',input_pair);
    
    % Find size of input arrays
    n1 = length(input1);
    n2 = length(input2);
    if (n1==n2)
        n = n1;
    else
        error('arrays input1 and input2 must have the same size')
    end
    
    if n==1 %single element input
        % Update AbstractState and extract
        calllib('coolprop','AbstractState_update',handle,INPUTS,input1,input2,ierr,herr,buffer_size)
        output1 = calllib('coolprop','AbstractState_keyed_output',handle,output,ierr,herr,buffer_size);
        
    else
        input1Ptr = libpointer('doublePtr',input1);
        input2Ptr = libpointer('doublePtr',input2);
        
        %Choosing parameter to compute
        outputs = calllib('coolprop','get_param_index',out1);
        
        %Creating ouput pointers
        out1Ptr = libpointer('doublePtr',zeros(n,1));
        
        calllib('coolprop','AbstractState_update_and_1_out',handle,INPUTS,input1Ptr,input2Ptr,n,outputs,out1Ptr,ierr,herr,buffer_size);
        
        %Saving computed values to array
        output1=get(out1Ptr,'Value');
        
    end
end

end

