function [output1,output2,output3,output4,output5] = CP5(input_pair,input1,input2,out1,out2,out3,out4,out5,handle)
% Wrapper for the low-level CoolProp function "AbstractState update and 5
% out"

% Usage example
% CP5(gas.handle,'PT_INPUTS',linspace(1e5,2e5,10)',linspace(300,400,10)','D','C','H','T','P')

% For enthalpy and pressure use: 'HmassP_INPUTS'

% Declaring error variables
ierr = 0; buffer_size = 10;
herr= char((1:1:buffer_size));
% Obtain input pair index
INPUTS = calllib('coolprop','get_input_pair_index',input_pair);

n = length(input1);
input1Ptr = libpointer('doublePtr',input1);
input2Ptr = libpointer('doublePtr',input2);

outputs=zeros(5,1);
%Choosing parameters to compute
outputs(1,1) = calllib('coolprop','get_param_index',out1);
outputs(2,1) = calllib('coolprop','get_param_index',out2);
outputs(3,1) = calllib('coolprop','get_param_index',out3);
outputs(4,1) = calllib('coolprop','get_param_index',out4);
outputs(5,1) = calllib('coolprop','get_param_index',out5);

%Creating ouput pointers
out1Ptr = libpointer('doublePtr',zeros(n,1));
out2Ptr = libpointer('doublePtr',zeros(n,1));
out3Ptr = libpointer('doublePtr',zeros(n,1));
out4Ptr = libpointer('doublePtr',zeros(n,1));
out5Ptr = libpointer('doublePtr',zeros(n,1));

calllib('coolprop','AbstractState_update_and_5_out',handle,INPUTS,input1Ptr,input2Ptr,n,outputs,out1Ptr,out2Ptr,out3Ptr,out4Ptr,out5Ptr,ierr,herr,buffer_size);

%Saving computed values to array
output1=get(out1Ptr,'Value');
output2=get(out2Ptr,'Value');
output3=get(out3Ptr,'Value');
output4=get(out4Ptr,'Value');
output5=get(out5Ptr,'Value');
end

