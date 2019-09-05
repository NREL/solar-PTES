% The following are extracts of outdated features of various subroutines of
% the PTES code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZATION OF HEX T-Q DIAGRAM
if WM==1 && fignum~=0
    % Set Cp and enthalpy arrays
    [TvH,CpvH,hvH] = hex_get_prop(fluidH,TC1,TH2,n);
    [TvC,CpvC,hvC] = hex_get_prop(fluidC,TC1,TH2,n);
    % Find DT distribution in heat space
    TC = zeros(1,n);
    TH = zeros(1,n);
    dQ = QT/(n-1);
    QS  = zeros(1,n);
    TC(1)  = TC1;
    TH(n) = TH2;
    for i=1:n-1
        TC(i+1) = TC(i) + dQ/(mC*rtab(TC(i),TvC,CpvC));
        TC(i+1) = TC(i) + dQ/(mC*rtab(0.5*(TC(i)+TC(i+1)),TvC,CpvC));
        QS(i+1) = QS(i) + dQ;
    end
    for i=n:-1:2
        TH(i-1) = TH(i) - dQ/(mH*rtab(TH(i),TvH,CpvH));
        TH(i-1) = TH(i) - dQ/(mH*rtab(0.5*(TH(i)+TH(i-1)),TvH,CpvH));
    end
    
    figure(fignum)
    fCname = fluidC.name; if strcmp(fCname,'mineral_oil'), fCname = 'Mineral Oil'; end
    fHname = fluidH.name; if strcmp(fHname,'mineral_oil'), fHname = 'Mineral Oil'; end
    if strcmp(fCname,'solar_salt'), fCname = 'Solar Salt'; end
    if strcmp(fHname,'solar_salt'), fHname = 'Solar Salt'; end
    LC = sprintf('%s, %.1f bar',fCname,pC_st/1e5);
    LH = sprintf('%s, %.1f bar',fHname,pH_st/1e5);
    if mod(fignum,2)==0
        plot(QS./QS(n),TC,QS./QS(n),TH)
        legend(LC,LH,'Location','SouthEast')
    else
        plot(QS./QS(n),TH,QS./QS(n),TC)
        legend(LH,LC,'Location','SouthEast')
    end
    ylabel('Temperature [K]')
    xlabel('Normalised heat transfer')
    %ylim([50,300])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display gas states
disp([[gas1(1,1:i).T]'-273.15,[gas1(1,1:i).p]'/1e5,(1:i)']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain linear index of interest and import fluid.state
idx   = sub2ind(size(fluid.state),ind(1),ind(2));
state = fluid.state(idx);
stage = fluid.stage(idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separately write and plot charge and discharge cycles
% Write charge cycle
if WM==1
    %WRITE_CHARGE
    %PLOT_CHARGE
end

% Write discharge cycle
if WM==1
    %WRITE_DISCHARGE
    %PLOT_DISCHARGE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find optimal pressure ratio for discharge cycle, PR_dis
% On main script:
optimise  = 1;  % run discharge cycle several times and find optimal PR_dis
tol = 1e-3;
switch optimise
    case 1
        %PTES_OPT_DISCHARGE; % run PTES discharge to find optimal PRdis
        save('./Outputs/WS.mat');
        [ PR_dis, ~, ~ ] = golden_search( @discharge_function, PRfactMin*PR, PRfactMax*PR, tol*PR, 'Min');
    case 0
        PR_dis = PRfact*PR;
    otherwise
        error('not implemented')
end
if any([PR_dis/PR < (PRfactMin + 2*tol),PR_dis/PR > (PRfactMax - 2*tol)])
    err_PR_opt = 1;
    warning('PTES_OPT_DISCHARGE didnt find optimal point!!!')
else
    err_PR_opt = 0;
end

% On ENERGY_BALANCE subroutine:
% Manage warnings
switch optimise
    case 1
        if err_PR_opt == 1
            warning('PTES_OPT_DISCHARGE didnt find optimal point!!!')
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute freezing curve of Ethylene Glycol as a function of concentration
T_freeze = zeros(1,57);
fluidCname = cell(1,57);
for i0 = 1:57
    if i0 < 10
        fluidCname{i0} = strcat('INCOMP::MEG2[0.0',num2str(i0-1),']');
    else
        fluidCname{i0} = strcat('INCOMP::MEG2[0.',num2str(i0-1),']');
    end
    T_freeze(i0) = CoolProp.Props1SI('T_freeze',fluidCname{i0}) - 273.15;
end
disp(fluidCname{i0})
disp(T_freeze)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading input variables from .txt file

% Read input variables from .txt file
fileID = fopen('inputs.txt','r');
fgetl(fileID); fgetl(fileID);
in_num(1:13)  = fscanf(fileID,'%*s %f %*s %*s',[1 Inf])';
fgetl(fileID); fgetl(fileID);
in_str{1} = fscanf(fileID,'%*s %s %*s',[1 1]); fgetl(fileID);
in_str{2} = fscanf(fileID,'%*s %s %*s',[1 1]); fgetl(fileID);
in_str{3} = fscanf(fileID,'%*s %s %*s',[1 1]); fgetl(fileID);
fgetl(fileID); fgetl(fileID); fgetl(fileID);
in_num(14:16) = fscanf(fileID,'%*s %f %*s %*s',[1 3])';

% Atmospheric conditions and cycle parameters
T0      = in_num(1) + 273.15;   % ambient temp, K
p0      = in_num(2) * 1e5;      % ambient pressure, Pa
TC_0    = in_num(3) + 273.15;   % temperature of discharged cold fluid
TH_0    = in_num(4) + 273.15;   % temperature of discharged hot fluid
mdot    = in_num(5);            % working fluid mass flow rate, kg/s
t_ch    = in_num(6) * 3600;     % charge time, s
pmax    = in_num(7) * 1e5;      % top pressure, Pa
PR      = in_num(8);            % charge pressure ratio
PR_fact = in_num(9);            % PR_dis = PR*PR_fact.
Nc_ch   = in_num(10);           % number of compressions during charge
Ne_ch   = in_num(11);           % number of expansions during charge
setTmax = in_num(12);           % set Tmax? (substitutes PR).
Tmax    = in_num(13) + 273.15;  % maximum temp at compressor outlet
gasName    = in_str{1}; % working fluid
fluidHname = in_str{2}; % hot storage fluid
fluidCname = in_str{3}; % cold storage fluid
Nc_dis = Ne_ch; % compressions during discharge
Ne_dis = Nc_ch; % expansions during discharge

% Set component parameters
eta    = in_num(14);  % polytropic efficiency
eff    = in_num(15);  % heat exchanger effectiveness
ploss  = in_num(16);  % pressure loss in HEXs

% Set operation modes
mode       = in_num(17); % cycle mode: 0=PTES, 1=Heat pump, 2=Heat engine
multi_run  = in_num(18); % run cycle several times with different parameters?
make_plots = in_num(19); % make plots?
savefig    = in_num(20); % save figures at the end?

% Variables to run cycle multiple times and plot curves. The variables must
% have been defined in the PTES_SET_MULTI_RUN script
if multi_run==1
    Vpnt = 'eta';   % variable along curve
    Npnt = 5;       % points on curve    
    pnt1 = 0.85;    % min value
    pnt2 = 0.95;    % max value
    Apnt = linspace(pnt1,pnt2,Npnt); % array
    Vcrv = 'eff';  % variable between curves
    Ncrv = 3;      % number of curves    
    crv1 = 0.95;   % min value
    crv2 = 0.99;   % max value
    Acrv = linspace(crv1,crv2,Ncrv); % array
    SET_MULTI_VAR;
else
    Npnt=1; Ncrv=1;
    % Start new logfile
    delete ./Outputs/log.txt
    diary  ./Outputs/log.txt
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading from CoolProp (low interface)

% Declaring error variables
ierr = 0; buffer_size = 10;
herr= char((1:1:buffer_size));
% Selecting input par values
PT_INPUTS = calllib('coolprop','get_input_pair_index','PT_INPUTS'); % CoolProp.get_input_pair_index('PT_INPUTS')

input1 = linspace(pressure,pressure,n)';
input2 = Tv;
input1Ptr = libpointer('doublePtr',input1);
input2Ptr = libpointer('doublePtr',input2);

outputs=zeros(5,1);
%Choosing parameters to compute
outputs(1,1) = calllib('coolprop','get_param_index','C');
outputs(2,1) = calllib('coolprop','get_param_index','H');
outputs(3,1) = calllib('coolprop','get_param_index','T');
outputs(4,1) = calllib('coolprop','get_param_index','T');
outputs(5,1) = calllib('coolprop','get_param_index','T');

%Creating ouput pointers
out1Ptr = libpointer('doublePtr',zeros(n,1));
out2Ptr = libpointer('doublePtr',zeros(n,1));
out3Ptr = libpointer('doublePtr',zeros(n,1));
out4Ptr = libpointer('doublePtr',zeros(n,1));
out5Ptr = libpointer('doublePtr',zeros(n,1));

calllib('coolprop','AbstractState_update_and_5_out',fluid.handle,PT_INPUTS,input1Ptr,input2Ptr,n,outputs,out1Ptr,out2Ptr,out3Ptr,out4Ptr,out5Ptr,ierr,herr,buffer_size);

%Saving computed values to array
Cpv=get(out1Ptr,'Value');
hv =get(out2Ptr,'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polytropic compression/expansion 
% Conventional approach
state.p = var;
pv = logspace(log10(p1),log10(state.p),num);
hv = zeros(1,num);
rhov = zeros(1,num);
hv(1) = h1;
rhov(1) = rho1;
for i0 = 1:(num-1)
    if (pv(i0+1)>0.99*Pcrit && pv(i0+1)<1.00*Pcrit), pv(i0+1) = Pcrit*0.99; end
    Dp = (pv(i0+1) - pv(i0));
    Dh = Dp/(rhov(i0)*eta^n);
    hv(i0+1)  = hv(i0) + Dh;
    rhov(i0+1) = CP1('HmassP_INPUTS',hv(i0+1),pv(i0+1),'D',fluid.handle);
    rhoAV = 0.5*(rhov(i0+1) + rhov(i0));
    Dh = Dp/(rhoAV*eta^n);
    hv(i0+1)  = hv(i0) + Dh;
end
state.h = hv(num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop inside PTES.m scrip, separating the PTES mode, the Heat Pump
% only mode, and the Heat Engine only mode.
tic %start timer

% Reinitialise arrays (gas, fluids and tanks) to zero and do other
% preliminary tasks
PTES_INITIALISE

for iL = 1:Load.num
    switch mode
        case 0 % PTES
            % Run charge cycle, compute state of storage tanks and run
            % discharge cycle
            switch Load.type(iL)
                case 'chg'
                    PTES_CHARGE
                    
                case 'str'
                    PTES_TANKS_STORAGE % no storage loss considered at the moment
                    
                case 'dis'
                    if optimise % obtain optimal PRr
                        f = @(PRr) ptes_discharge_function(gas, fluidH, fluidC, HT, CT, environ,...
                            T0, T1, pbot, PRr, PRch, Load.mdot(iL), Nc_dis, Ne_dis,...
                            eta, eff, ploss, Load.time(iL), mode);
                        
                        [PRr,ineff,xv,yv,iter] = golden_search(f,PRr_min,PRr_max,0.005,'Min',100);
                    end
                    PTES_DISCHARGE
            end
        case 1 % Heat pump only
            PTES_CHARGE
            
        case 2 % Heat engine only
            PTES_SOLAR_TANKS
            if optimise % obtain optimal PRr
                f = @(PRr) ptes_discharge_function(gas, fluidH, fluidC, HT, CT, environ,...
                    T0, T1, pbot, PRr, PRch, Load.mdot(iL), Nc_dis, Ne_dis,...
                    eta, eff, ploss, Load.time(iL), mode);
                
                [PRr,ineff,xv,yv,iter] = golden_search(f,PRr_min,PRr_max,0.005,'Min',100);
            end
            PTES_DISCHARGE
    end
end

% Compute energy balance
PTES_ENERGY_BALANCE

toc %stop timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

