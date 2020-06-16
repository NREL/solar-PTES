% Exploring the possibility of creating our own property tables...
% To run this script, one must first run the first few lines of the PTES.m
% script (to add the paths and load CoolProp).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables
clear;

% Enter debugging mode if an error occurs
dbstop if error
%dbclear all

inputs = 'HmassP_INPUTS';
%inputs = 'PQ_INPUTS';

%%% CREATE TABLE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
create = 0;
create_only = 0;

switch create_only
    case 1
        tic
        LIB.Water.HmassP_INPUTS = LIB_THERM('Water','HmassP_INPUTS',100);
        toc
        tic
        LIB.Water.PQ_INPUTS     = LIB_THERM('Water','PQ_INPUTS',1e3);
        toc
        LIB.Water.name = 'Water';
        save('./LIB/LIB/LIB.mat','LIB')
        return
end

switch create
    case 1
        switch inputs
            case 'HmassP_INPUTS'
                tic
                % Create table
                LIB.Water.(inputs) = LIB_THERM('Water','HmassP_INPUTS',50);
                LIB.Water.name     = LIB.Water.(inputs).name;
                
                % Save library into a .mat file
                save('./LIB/LIB/LIB.mat','LIB')
                toc
                
                % Set list of properties to extract
                hp_list = {'T','S','D','U','Q','C','CONDUCTIVITY','VISCOSITY','PRANDTL'};
                list = hp_list(1:5);
                
                tic
                % Create denser table for comparisson of interpolated values.
                % Obtain exact values from Coolprop HEOS.
                X  = LIB.Water.(inputs).X;
                Y  = LIB.Water.(inputs).Y;
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
                
                fluid = fluid_class(LIB.Water.name,'WF','CP','HEOS',1,1);
                zcool = zeros([size(Xq),length(list)]);
                for ip=1:length(list)
                    for ic=1:size(Xq,2)
                        zcool(:,ic,ip)  = CP1('HmassP_INPUTS',Xq(:,ic),Yq(:,ic),list{ip},fluid.handle);
                    end
                end
                toc
                
            case 'PQ_INPUTS'
                tic
                % Create table
                LIB.Water.(inputs) = LIB_THERM('Water','PQ_INPUTS',1e3);
                LIB.Water.name = LIB.Water.(inputs).name;
                
                % Save library into a .mat file
                save('./LIB/LIB/LIB.mat','LIB')
                toc
                
                % Set list of properties to extract
                pq_list = {'H','T','S','D','U','C','CONDUCTIVITY','VISCOSITY','PRANDTL'};
                list = pq_list;%(1:5);
                
                tic
                % Create denser table for comparisson of interpolated values.
                % Obtain exact values from Coolprop HEOS.
                xL   = LIB.Water.(inputs).xL;
                num  = round(length(xL)/20);
                xmin = min(xL);
                xmax = max(xL);
                nq   = 3*num-2;
                xq   = logspace(log10(xmin),log10(xmax),nq);
                yq   = linspace(0,1,nq);
                [Xq,Yq] = meshgrid(xq,yq);
                
                fluid = fluid_class(LIB.Water.name,'WF','CP','HEOS',1,1);
                zcool = zeros([size(Xq),length(list)]);
                for ip=1:length(list)
                    for ic=1:size(Xq,2)
                        zcool(:,ic,ip)  = CP1('PQ_INPUTS',Xq(:,ic),Yq(:,ic),list{ip},fluid.handle);
                    end
                end
                toc
                
            otherwise
                error('not implemented');
        end
        
        % Save query table into a .mat file
        save('./LIB/LIB/Query.mat','zcool','list','Xq','Yq','inputs');
end

%%% OBTAIN PROPERTIES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables
clear;

% Load library
load('./LIB/LIB/LIB.mat')
load('./LIB/LIB/Query.mat')

% Restrict list of variables to be compared
%list  = list(1:5);
%zcool = zcool(:,:,1:5);

tic
for i0=1:1
% Obtain interpolated tables
[LIB.Water.(inputs), ztab] = RLIB(LIB.Water.(inputs), inputs, Xq, Yq, list);
end
toc

% Compare against results from CoolProp using TTSE
fluidT = fluid_class(LIB.Water.name,'WF','CP','TTSE',1,1);
tic
for i0=1:1
zcoolT = zeros([size(Xq),length(list)]);
for ip=1:length(list)
    for ic=1:size(Xq,2)
        switch inputs
            case 'HP'
                zcoolT(:,ic,ip) = CP1('HmassP_INPUTS',Xq(:,ic),Yq(:,ic),list{ip},fluidT.handle);
            case 'PQ'
                zcoolT(:,ic,ip) = CP1('PQ_INPUTS',Xq(:,ic),Yq(:,ic),list{ip},fluidT.handle);
        end
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
fprintf(1,'Max error    = %8.3f   %%\n',max(err,[],'all','omitnan'))
fprintf(1,'Mean error   = %10.1e %%\n',mean(err,'all','omitnan'))
fprintf(1,'Median error = %10.1e %%\n',median(err,'all','omitnan'))
fprintf(1,'\n')
fprintf(1,'TTSE interpolation:\n')
fprintf(1,'Max error    = %8.3f   %%\n',max(err2,[],'all','omitnan'))
fprintf(1,'Mean error   = %10.1e %%\n',mean(err2,'all','omitnan'))
fprintf(1,'Median error = %10.1e %%\n',median(err2,'all','omitnan'))
fprintf(1,'\n\n')

% Plot grids
switch inputs
    case 'HmassP_INPUTS'
        TAB = LIB.Water.(inputs);
        X   = TAB.X;
        Y   = TAB.Y;
        XL  = TAB.XL;
        YL  = TAB.YL;
        XG  = TAB.XG;
        YG  = TAB.YG;
        XT  = TAB.XT;
        YT  = TAB.YT;
        xL  = TAB.xL;
        xG  = TAB.xG;
        yL  = TAB.yL;
        yG  = TAB.yG;
        dxL = TAB.dxL;
        dxG = TAB.dxG;
        figure(1)
        semilogy(xL,yL,'b',xG,yG,'b'); hold on;
        semilogy(xL-dxL,yL,xG+dxG,yG)
        semilogy(X,Y,'k.');
        semilogy(XL,YL,'k.');
        semilogy(XG,YG,'k.');
        semilogy(XT,YT,'k.');
        botm = 0.01;
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
    case 'PQ_INPUTS'
        TAB  = LIB.Water.(inputs);
        xL   = TAB.xL;
        xq   = Xq(1,:);
        botm = 0.01;
        mark = max(err(:,:,:),[],3) > botm;
        figure(1)
        loglog(xL,xL,'ks'); hold on;
        loglog(xq,xq,'ro'); hold off;
        figure(2)
        semilogx(Xq,Yq,'k.'); hold on;
        semilogx(Xq(mark),Yq(mark),'rs');
        hold off;
        for i=1:length(list)
            plot_error(err(:,:,i),[botm,10*botm],'Pressure','Vapour quality',10+i,list{i})
        end
    otherwise
        error('not implemented')
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
