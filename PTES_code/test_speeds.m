% test_speeds

N = 1e5 ;
cp = 1000 ;
tic
for i = 1:N
    %h1 = RPN('PT_INPUTS',1e5,400,'H',fluidH) ; % Interpolation of tables
    %h1 = RPN('PT_INPUTS',10e5,400,'H',gas) ; % Calls coolprop
    h1 = 400 * cp ; % Calculates enthalpy directly
end
toc
