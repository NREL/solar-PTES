% Exploring the possibility of creating our own property tables...
% To run this script, one must first run the first few lines of the PTES.m
% script (to add the paths and load CoolProp).


% TO DO:
% - Work on implementation. How to move beyond h-p? Find out why TTSE was
% failing on PTES code!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables
clear;

% Enter debugging mode if an error occurs
dbstop if error
%dbclear all

%%% CREATE TABLE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
create = 0;
switch create
    case 1
        tic
        Water.HP = LIB_THERM('Water','HmassP_INPUTS',100);
        Water.name = Water.HP.name;
        
        % Save library into a .mat file
        save('./Other/LIB/Water.mat','Water')
        toc
        
        % Set list of properties to extract
        hp_list = {'T','S','D','U','Q','C','CONDUCTIVITY','VISCOSITY','PRANDTL'};
        list = hp_list(1:5);
        
        tic
        % Create denser table for comparisson of interpolated values.
        % Obtain exact values from Coolprop HEOS.
        X  = Water.HP.X;
        Y  = Water.HP.Y;
        x  = X(1,:);
        y  = Y(:,1)';
        xmin = min(x);
        xmax = max(x);
        ymin = min(y);
        ymax = max(y);
        nq   = 3*length(x)-2;
        xq   = linspace(xmin,xmax,nq);
        yq   = logspace(log10(ymin),log10(ymax),nq);
        [Xq,Yq] = meshgrid(xq,yq);
        
        %figure(4)
        %semilogy(X,Y,'k.'); hold on;
        %semilogy(Xq,Yq,'rs'); hold off;
        %keyboard
        
        fluid = fluid_class(Water.name,'WF','CP','HEOS',1,1);
        zcool = zeros([size(Xq),length(list)]);
        for ip=1:length(list)
            for ic=1:size(Xq,2)
                zcool(:,ic,ip)  = CP1('HmassP_INPUTS',Xq(:,ic),Yq(:,ic),list{ip},fluid.handle);
            end
        end
        
        % Save table into a .mat file
        save('./Other/LIB/Water.mat','Water','zcool','list','Xq','Yq')
        toc
end

%%% OBTAIN PROPERTIES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables
clear;

% Load library
load('./Other/LIB/Water.mat')

tic
for i0=1:1
% Obtain interpolated tables
[Water.HP, ztab] = RLIB(Water.HP, list, Xq, Yq);
end
toc

% Compare against results from CoolProp using TTSE
fluidT = fluid_class(Water.name,'WF','CP','TTSE',1,1);
tic
for i0=1:1
zcoolT = zeros([size(Xq),length(list)]);
for ip=1:length(list)
    for ic=1:size(Xq,2)
        zcoolT(:,ic,ip) = CP1('HmassP_INPUTS',Xq(:,ic),Yq(:,ic),list{ip},fluidT.handle);
    end
end
QT = zcoolT(:,:,5);
QT(QT < 0) = -1.0;
zcoolT(:,:,5) = QT;
end
toc

% Compute and print errors
err  = abs((zcool - ztab)./zcool)*100;
err2 = abs((zcool - zcoolT)./zcool)*100;
fprintf(1,'CUSTOM interpolation:\n')
fprintf(1,'Max error    = %8.3f   %%\n',max(err,[],'all'))
fprintf(1,'Mean error   = %10.1e %%\n',mean(err,'all'))
fprintf(1,'Median error = %10.1e %%\n',median(err,'all'))
fprintf(1,'\n')
fprintf(1,'TTSE interpolation:\n')
fprintf(1,'Max error    = %8.3f   %%\n',max(err2,[],'all'))
fprintf(1,'Mean error   = %10.1e %%\n',mean(err2,'all'))
fprintf(1,'Median error = %10.1e %%\n',median(err2,'all'))
fprintf(1,'\n\n')

% Plot grids
X   = Water.HP.X;
Y   = Water.HP.Y;
XL  = Water.HP.XL;
YL  = Water.HP.YL;
XG  = Water.HP.XG;
YG  = Water.HP.YG;
XT  = Water.HP.XT;
YT  = Water.HP.YT;
xL  = Water.HP.xL;
xG  = Water.HP.xG;
yL  = Water.HP.yL;
yG  = Water.HP.yG;
dxL = Water.HP.dxL;
dxG = Water.HP.dxG;
figure(1)
semilogy(xL,yL,'b',xG,yG,'b'); hold on;
semilogy(xL-dxL,yL,xG+dxG,yG)
semilogy(X,Y,'k.');
semilogy(XL,YL,'k.');
semilogy(XG,YG,'k.');
semilogy(XT,YT,'k.');
botm = 0.1;
mark = max(err(:,:,:),[],3) > botm;
semilogy(Xq(mark),Yq(mark),'rs');
hold off;
xlabel('Enthalpy, J/kg/K')
ylabel('Pressure, bar')
xlim([min(X,[],'all'),max(X,[],'all')])
ylim([min(Y,[],'all'),max(Y,[],'all')])

figure(2)
semilogy(xL,yL,xG,yG,xL-dxL,yL,xG+dxG,yG); hold on;
semilogy(X,Y,'k.',Xq,Yq,'r.');
hold off;
xlabel('Enthalpy, J/kg/K')
ylabel('Pressure, bar')

% Invert vertical axis for plotting
err(:,:,:)  = err(end:-1:1,:,:);
for i=1:length(list)
    plot_error(err(:,:,i),[botm,10*botm],'Enthalpy','Pressure',10+i,list{i})
end


%%% SUPPORT FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_error (ERR,range,xname,yname,fignum,property)

figure(fignum)
clf(fignum)
imagesc(ERR,range)
colorbar
xlabel(xname)
ylabel(yname)
str = {'Errors:',...
    [sprintf('Max     = %.2g',max(ERR,[],'all')), ' %'],...
    [sprintf('Mean   = %.2g',mean(ERR,'all')), ' %']...
    [sprintf('Median= %.2g',median(ERR,'all')) ' %']};
dim = [0.52 0.1 0.4 0.4]; %[x y w h]
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','white');

switch property
    case 'T'
        tit = 'Temperature [T]';
    case 'H'
        tit = 'Enthalpy [H]';
    case 'S'
        tit = 'Entropy [S]';
    case 'D'
        tit = 'Density [D]';
    case 'U'
        tit = 'Internal Energy [U]';
    case 'C'
        tit = 'Isobaric heat capacity [Cp]';
    case 'CONDUCTIVITY'
        tit = 'Conductivity [k]';
    case 'VISCOSITY'
        tit = 'Viscosity [mu]';
    case 'PRANDTL'
        tit = 'Prandtl [Pr]';
    case 'Q'
        tit = 'Vapour quality [Q]';
    otherwise
        error('not implemented')
end
title(tit)

end
