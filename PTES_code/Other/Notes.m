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
hv =get(out2Ptr,'Value'); %#ok<*NASGU>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polytropic compression/expansion 
% Conventional approach
state.p = var; %#ok<*LTARG>
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

% Update flow area of stream 1 considering discrete number of tubes, and
% use number of tubes and tube thickness to compute heat transfer area of
% stream 2 (shell-side)

% Compute flow area
S1.Af  = S1.mdot/S1.G; % total flow area (stream 1)
Af1one = pi/4*S1.D^2;    % flow area of one channel
N1     = round(S1.Af/Af1one); % number of channels
S1.Af  = N1*Af1one;    % update total flow area

% Update G, Re, Cf and St
S1.G  = S1.mdot/S1.Af;
S1.Re = S1.G*S1.D/S1.mu;
[S1.Cf, S1.St] = developed_flow(S1.Re,S1.Pr,'circular');

% Set t1 (tube thickness)
t1_hoop = S1.p*S1.D/(2*sigma); %thickness for which sigma = hoop stress.
t1 = max([t1_min, 2*t1_hoop]);

% Compute A2
S2.A = L*N1*pi()*(S1.D+t1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EQUAL_CP ALGORITHM FOR HEX FUNCTION
switch algorithm
    case 'DT_min'
        
        % Compute hH1 for which DTmin = 0, using golden search method
        f1   = @(hH1) DTmin(mH,mC,hH2,hC1,hvH,TvH,hvC,TvC,n,'hH1',hH1);
        TolX = (hH2-hH1_min)/1e4; %tolerance
        options = optimset('TolX',TolX);%,'Display','iter');
        hH1  = fminbnd(f1,hH1_min,hH2,options);
        
        % Compute QMAX
        QMAX = mH*(hH2 - hH1);
        
    case 'Equal_Cp'
        
        % Compute QMAX using the fact that mH*CpH = mC*CpC for any intermediate
        % pinch point
        [~, CpH] = get_Cp( fluidH, TC1, TH2, pH2, n );
        [Tv,CpC] = get_Cp( fluidC, TC1, TH2, pC1, n );
        CH  = mH*CpH;
        CC  = mC*CpC;
        fx  = CH - CC;
        fx0 = fx(1);
        np  = 0; %number of pinch points found
        ip  = zeros(n-1,1); %pinch point index
        for i0=2:n-1
            if fx(i0)*fx0<0
                np = np + 1;
                ip(np) = i0; %pinch point index
                fx0 = fx(i0);
            end
        end
        Tp = zeros(np,1);
        Qp = zeros(size(Tp));
        for i0=1:np
            Tp(i0) = Tv(ip(i0));
            QC = mC*(rtab1(TvC,hvC,Tp(i0),1)-hC1); % mC*DhC from TC1 to Tp
            QH = mH*(hH2-rtab1(TvH,hvH,Tp(i0),1)); % mH*DhH from Tp to TH2
            Qp(i0) = QC + QH;
        end
        % QMAX has to be the minimum between Qp, QC and QH:
        QMAX = min([Qp;QMAX0]);
end

function [ Tv, Cpv ] = get_Cp( fluid, T1, T2, pressure, n )
% Obtain the Tv and Cpv arrays of a given fluid for the hex subroutines.
% Data ordered in regular intervals of T. Cp as a function of T.

if strcmp(fluid.read,'CP') %read from CoolProp
    
    Tv   = linspace(T1,T2,n)'; % Temperature array between TC1 and TH2
    pv  = ones(size(Tv)).*pressure; % pressure array
    dT   = (T2-T1)/(n*2);
    h_up = CP1('PT_INPUTS',pv,Tv+dT,'H',fluid.handle);
    h_lo = CP1('PT_INPUTS',pv,Tv-dT,'H',fluid.handle);
    Cpv  = (h_up - h_lo)/(2*dT);
    
elseif strcmp(fluid.read,'TAB') %read from table
    
    T   = fluid.TAB(:,1);
    Cp  = fluid.TAB(:,5);    
    Tv  = linspace(T1,T2,n)'; % Temperature array between TC1 and TH2
    Cpv = rtab_1D_1out(T,Cp,Tv,0);
else
    error('not implemented')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FITTING USING LINEAR TRANSFORMATIONS

% Obtain fitting coefficients
fit_mode = 1; %#ok<*FISST>
x1 = interp_coeffs(h1,T1,fit_mode);
x2 = interp_coeffs(h2,T2,fit_mode);
x3 = interp_coeffs(h3,T3,fit_mode);

% Obtain fitting arrays
y1  = linspace(min(h1),max(h1),10)';
Tf1 = interp_poly(x1,y1,fit_mode);
y2  = linspace(min(h2),max(h2),10)';
Tf2 = interp_poly(x2,y2,fit_mode);
y3  = linspace(min(h3),max(h3),10)';
Tf3 = interp_poly(x3,y3,fit_mode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Logarithmically spaced arrays
PM   = logspace(log10(Pmin),log10(Pcrit+dP),length(Psat))';
% Adapt the pressure array to make sure that it has one point at Pcrit.
ip   = find(PM>Pcrit,1,'first');
Pr   = (Pcrit/Pmin)^(1/(ip-1));
PM   = Pmin*Pr.^(0:length(Psat));
% Adapt the pressure array to make sure that it has one point at Pcrit.
ip   = find(PM>Pcrit,1,'first');
PM1  = logspace(log10(Pmin),log10(Pcrit),ip);
Pr   = PM1(2)/PM1(1);
PM2  = Pcrit*Pr.^(1:(length(Psat)-ip));
PM   = [PM1,PM2]';
% xM grid, containing one vertical line at the point separating HsatL and
% HsatG
a  = linspace(xL(iM),HsatL(end),nsx+1);
da = a(end) - a(end-1);
b  = linspace(HsatL(end)+da,xR(iM),nsx-1);
xM(iM,:) = [a,b];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check difference between CoolProp handles
output1 = CP1(input_pair,input1,input2,out1,handle);
output2 = CP1(input_pair,input1,input2,out1,fluid.HEOS);
if any(abs(output1./output2 - 1) > 1e-4)
    warning('Original CoolProp handle proving innacurate')
    keyboard
    output1 = output2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%