% Exploring the possibility of creating our own property tables...
% To run this script, one must first run the first few lines of the PTES.m
% script (to add the paths and load CoolProp).


% TODO:
% (0) improve documentation
% (1) try to speed up table creation by using PARFOR loops.
% (2) create T-s diagram and plot errors there (including saturation curve)
% (3) find way to obtain library file from fluid's name (perhaps just name
% the variable file after the fluid?)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables
clear;

% Set fluid
name = 'Water';
fluid = fluid_class(name,'WF','CP','HEOS',1,1);

switch name
    case 'Water'
        % Set minimum and maximum temperatures, pressures and enthalpies
        Tmin  = CP1(0,0,0,'Tmin',fluid.handle);
        Tmax  = 1500;
        pmin  = CP1(0,0,0,'pmin',fluid.handle);
        pmax  = 1000*1e5;
        hmin  = CP1('PT_INPUTS',pmax,Tmin,'H',fluid.handle);
        hmax  = CP1('PT_INPUTS',pmin,Tmax,'H',fluid.handle);
    otherwise
        error('not implemented')
end

% Set number of rows and number of columns
num = 300;

% Set arrays
p_array = logspace(log10(pmin),log10(pmax),num);
h_array = linspace(hmin,hmax,num);
T_array = linspace(Tmin,Tmax,num);

% Set list of properties to extract
hp_list = {'T','S','D','C','CONDUCTIVITY','VISCOSITY','PRANDTL','Q'};
pT_list = {'H','S','D','C','CONDUCTIVITY','VISCOSITY','PRANDTL','Q'};

%%% FILL TABLE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HP = makePROPS(h_array,p_array,'HmassP_INPUTS',hp_list,fluid.handle);
PT = makePROPS(p_array,T_array,'PT_INPUTS',pT_list,fluid.handle);
LIB.name = name;
LIB.HP = HP;
LIB.PT = PT;

% Save library into a .mat file
save('./Other/LIB.mat','LIB')

%% OBTAIN PROPERTIES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load library
load('./Other/LIB.mat')

% Interpolation arrays
nq = 33;
p_arr = logspace(log10(pmin),log10(pmax),nq);
h_arr = linspace(hmin,hmax,nq);
[H_arr,P_arr] = meshgrid(h_arr,p_arr);

D_cool = zeros(size(H_arr));
for ic=1:size(H_arr,2)
    D_cool(:,ic) = CP1('HmassP_INPUTS',H_arr(:,ic),P_arr(:,ic),'D',fluid.handle);
end

D_tab = RTAB(LIB.HP, 'D', H_arr, P_arr, 'makima');
err =  abs((D_cool - D_tab)./D_cool)*100;
err_max = max(err,[],'all');

figure(1)
surf(H_arr,P_arr,D_tab)
xlabel('Enthalpy')
ylabel('Pressure')
zlabel('Density')
title([sprintf('Max error = %.3f ',err_max),'$$\%$$']);

figure(2)
imagesc(err,[0.1,10])
colorbar




%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PROPS = makePROPS(x_array,y_array,inputs,list,handle)

[X,Y] = meshgrid(x_array,y_array);
TAB    = zeros([size(X),length(list)]);
for ip=1:length(list)
    for ic=1:size(X,2) %fill column by column
        TAB(:,ic,ip)   = CP1(inputs,X(:,ic),Y(:,ic),list{ip},handle);
    end
end
PROPS.TAB  = TAB;
PROPS.X    = X;
PROPS.Y    = Y;
PROPS.list = list;

end

function result = RTAB (PROPS, out, xq, yq, method)

% Obtain index corresponding to property 'out'
ind = strcmp(out,PROPS.list);
if ~any(ind)
    error('invalid property selected')
end

result = interp2(PROPS.X,PROPS.Y,PROPS.TAB(:,:,ind),xq,yq,method);

end

