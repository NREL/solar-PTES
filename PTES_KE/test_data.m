%close all
global  DEXP1 DEXP2 DEXP3 DEXP4 DCOND DCMP1 DCMP2 DCMP3 DCMP4 FWH1 FWH2 ...
        PREH BOILER SUPH REH CYC data %#ok<NUSED>
    
%% Cycle Metrics
cyc_dat.name    = "Steam Cycle";
CYC             = CYCLE(cyc_dat);
% Work in (consider split steam flow)
DCMP1.W         = DCMP1.w*DCMP1.in.mdot / MEG;                          % Work input of pump1 (MW)
DCMP2.W         = DCMP2.w*DCMP2.in.mdot / MEG;                          % Work input of pump2 (MW)
DCMP3.W         = DCMP3.w*DCMP3.in.mdot / MEG;                          % Work input of pump3 (MW)
CYC.WinD        = (DCMP1.W + DCMP2.W + DCMP3.W);                        % Total work input (MW)
%Work out
DEXP3.W         = DEXP3.w*DEXP3.in.mdot / MEG;                          % Work output of expander3 (MW)
DEXP2.W         = DEXP2.w*DEXP2.in.mdot / MEG;                          % Work output of expander2 (MW)
DEXP1.W         = DEXP1.w*DEXP1.in.mdot / MEG;                          % Work output of expander1 (MW)
CYC.WoutD       = (DEXP3.W + DEXP2.W + DEXP3.W);                        % Total work output (MW)
% Net work (power)
CYC.WnetD       = (CYC.WinD + CYC.WoutD);                               % Net power generated by cycle (MW)
% solar inputs (reheater only reheats split after HP turbine)
REH.Qc          = REH.qc * REH.Cin.mdot / MEG;                          % reheater heat input (MW)    
PREH.Qc         = PREH.qc * PREH.Cin.mdot / MEG;                        % preheater heat input (MW)   
BOILER.Qc      	= BOILER.qc * BOILER.Cin.mdot / MEG;                    % boiler heat input (MW)
SUPH.Qc         = SUPH.qc * SUPH.Cin.mdot / MEG;                        % superheat heat input (MW)
CYC.Qin         = (PREH.Qc + BOILER.Qc + SUPH.Qc + REH.Qc);             % Total heat input (MW)
CYC.Qrej        = DCOND.qh*DCOND.Hin.mdot / MEG;                        % Rejected heat by condenser(MW)
% cycle efficiency
CYC.eta         = 100. * CYC.WnetD / CYC.Qin;            	% Cycle Efficinecy (%)
%% Create array for the cycle data points
DIS(1,1:9)  = ["Point","Component","T, K","P, bar","h, kJ/kg","s, kJ/kg.K","Mass flow, kg/s","q","Density, kg/m3"];
% HP Expander
DIS(2,1:9)  = [1, "EXP1 in",    DEXP1.in.T - degC, DEXP1.in.P/BAR, DEXP1.in.h/KIL, DEXP1.in.s/KIL, DEXP1.in.mdot, DEXP1.in.q, DEXP1.in.rho] ;
DIS(3,1:9)  = [2, "EXP1 out",	DEXP1.out.T - degC, DEXP1.out.P/BAR, DEXP1.out.h/KIL, DEXP1.out.s/KIL, DEXP1.out.mdot, DEXP1.out.q, DEXP1.out.rho] ;
% Solar Reheat
DIS(4,1:9)  = [3, "REH in",	REH.Cin.T - degC, REH.Cin.P/BAR, REH.Cin.h/KIL, REH.Cin.s/KIL, REH.Cin.mdot, REH.Cin.q, REH.Cin.rho] ;
DIS(5,1:9)  = [4, "REH out",REH.Cout.T - degC, REH.Cout.P/BAR, REH.Cout.h/KIL, REH.Cout.s/KIL, REH.Cout.mdot, REH.Cout.q, REH.Cout.rho] ;
% IP Expander
DIS(6,1:9)  = [5, "EXP2 in",    DEXP2.in.T - degC, DEXP2.in.P/BAR, DEXP2.in.h/KIL, DEXP2.in.s/KIL, DEXP2.in.mdot, DEXP2.in.q, DEXP2.in.rho] ;
DIS(7,1:9)  = [6, "EXP2 out",   DEXP2.out.T - degC, DEXP2.out.P/BAR, DEXP2.out.h/KIL, DEXP2.out.s/KIL, DEXP2.out.mdot, DEXP2.out.q, DEXP2.out.rho] ;
% LP Expander
DIS(8,1:9)  = [7, "EXP3 in",    DEXP3.in.T - degC, DEXP3.in.P/BAR, DEXP3.in.h/KIL, DEXP3.in.s/KIL, DEXP3.in.mdot, DEXP3.in.q, DEXP3.in.rho] ;
DIS(9,1:9)  = [8, "EXP3 out",   DEXP3.out.T - degC, DEXP3.out.P/BAR, DEXP3.out.h/KIL, DEXP3.out.s/KIL, DEXP3.out.mdot, DEXP3.out.q, DEXP3.out.rho] ;
% Condenser
DIS(10,1:9) = [9,  "COND in",   DCOND.Hin.T - degC, DCOND.Hin.P/BAR, DCOND.Hin.h/KIL, DCOND.Hin.s/KIL, DCOND.Hin.mdot, DCOND.Hin.q, DCOND.Hin.rho] ;
DIS(11,1:9) = [10, "COND out",  DCOND.Hout.T - degC, DCOND.Hout.P/BAR, DCOND.Hout.h/KIL, DCOND.Hout.s/KIL, DCOND.Hout.mdot, DCOND.Hout.q, DCOND.Hout.rho] ;
% Compressor (Pump) 1
DIS(12,1:9) = [11, "CMP1 in",   DCMP1.in.T - degC, DCMP1.in.P/BAR, DCMP1.in.h/KIL, DCMP1.in.s/KIL, DCMP1.in.mdot, DCMP1.in.q, DCMP1.in.rho] ;
DIS(13,1:9) = [12, "CMP1 out",  DCMP1.out.T - degC, DCMP1.out.P/BAR, DCMP1.out.h/KIL, DCMP1.out.s/KIL, DCMP1.out.mdot, DCMP1.out.q, DCMP1.out.rho] ;
% FWH 1
DIS(14,1:9) = [13, "FWH1 Cin",  FWH1.Cin.T - degC, FWH1.Cin.P/BAR, FWH1.Cin.h/KIL, FWH1.Cin.s/KIL, FWH1.Cin.mdot, FWH1.Cin.q, FWH1.Cin.rho] ; 
DIS(15,1:9) = [14, "FWH1 Hin",  FWH1.Hin.T - degC, FWH1.Hin.P/BAR, FWH1.Hin.h/KIL, FWH1.Hin.s/KIL, FWH1.Hin.mdot, FWH1.Hin.q, FWH1.Hin.rho] ; 
DIS(16,1:9) = [15, "FWH1 out",  FWH1.Fout.T - degC, FWH1.Fout.P/BAR, FWH1.Fout.h/KIL, FWH1.Fout.s/KIL, FWH1.Fout.mdot, FWH1.Fout.q, FWH1.Fout.rho] ; 
% Compressor (Pump) 2
DIS(17,1:9) = [16, "CMP2 in",   DCMP2.in.T - degC, DCMP2.in.P/BAR, DCMP2.in.h/KIL, DCMP2.in.s/KIL, DCMP2.in.mdot, DCMP2.in.q, DCMP2.in.rho] ;
DIS(18,1:9) = [17, "CMP2 out",  DCMP2.out.T - degC, DCMP2.out.P/BAR, DCMP2.out.h/KIL, DCMP2.out.s/KIL, DCMP2.out.mdot, DCMP2.out.q, DCMP2.out.rho] ;
% FWH 2
DIS(19,1:9) = [18, "FWH2 Cin",  FWH2.Cin.T - degC, FWH2.Cin.P/BAR, FWH2.Cin.h/KIL, FWH2.Cin.s/KIL, FWH2.Cin.mdot, FWH2.Cin.q, FWH2.Cin.rho] ; 
DIS(20,1:9) = [19, "FWH2 Hin",  FWH2.Hin.T - degC, FWH2.Hin.P/BAR, FWH2.Hin.h/KIL, FWH2.Hin.s/KIL, FWH2.Hin.mdot, FWH2.Hin.q, FWH2.Hin.rho] ; 
DIS(21,1:9) = [20, "FWH2 out",  FWH2.Fout.T - degC, FWH2.Fout.P/BAR, FWH2.Fout.h/KIL, FWH2.Fout.s/KIL, FWH2.Fout.mdot, FWH2.Fout.q, FWH2.Fout.rho] ; 
% Compressor (Pump) 3
DIS(22,1:9) = [21, "CMP3 in",   DCMP3.in.T - degC, DCMP3.in.P/BAR, DCMP3.in.h/KIL, DCMP3.in.s/KIL, DCMP3.in.mdot, DCMP3.in.q, DCMP3.in.rho] ;
DIS(23,1:9) = [22, "CMP3 out",  DCMP3.out.T - degC, DCMP3.out.P/BAR, DCMP3.out.h/KIL, DCMP3.out.s/KIL, DCMP3.out.mdot, DCMP3.out.q, DCMP3.out.rho] ;
% Solar Preheat
DIS(24,1:9) = [23, "PREH in"	 PREH.Cin.T - degC, PREH.Cin.P/BAR, PREH.Cin.h/KIL, PREH.Cin.s/KIL, PREH.Cin.mdot, PREH.Cin.q, PREH.Cin.rho] ;
DIS(25,1:9) = [24, "PREH out", PREH.Cout.T - degC, PREH.Cout.P/BAR, PREH.Cout.h/KIL, PREH.Cout.s/KIL, PREH.Cout.mdot, PREH.Cout.q, PREH.Cout.rho] ;
% Solar generation
DIS(26,1:9) = [25, "BOILER in", BOILER.Cin.T - degC, BOILER.Cin.P/BAR, BOILER.Cin.h/KIL, BOILER.Cin.s/KIL, BOILER.Cin.mdot, BOILER.Cin.q, BOILER.Cin.rho] ;
DIS(27,1:9) = [26, "BOILER out",BOILER.Cout.T - degC, BOILER.Cout.P/BAR, BOILER.Cout.h/KIL, BOILER.Cout.s/KIL, BOILER.Cout.mdot, BOILER.Cout.q, BOILER.Cout.rho] ;
% Solar superheat
DIS(28,1:9) = [27, "SUPH in", SUPH.Cin.T - degC, SUPH.Cin.P/BAR, SUPH.Cin.h/KIL, SUPH.Cin.s/KIL, SUPH.Cin.mdot, SUPH.Cin.q, SUPH.Cin.rho] ;
DIS(29,1:9) = [28, "SUPH out",SUPH.Cout.T - degC, SUPH.Cout.P/BAR, SUPH.Cout.h/KIL, SUPH.Cout.s/KIL, SUPH.Cout.mdot, SUPH.Cout.q, SUPH.Cout.rho] ;
%Completed cycle
DIS(30,1:9) = [1, "EXP1 in",    DEXP1.in.T - degC, DEXP1.in.P/BAR, DEXP1.in.h/KIL, DEXP1.in.s/KIL, DEXP1.in.mdot, DEXP1.in.q, DEXP1.in.rho] ;

DIS(35,1:3) = ["Component" , "Work in/out (MW)", "Heat in (MW)"];
DIS(36,1:3) = ["EXP1", DEXP1.w * DEXP1.in.mdot / MEG, 0]; 
DIS(37,1:3) = ["EXP2", DEXP2.w * DEXP2.in.mdot / MEG, 0]; 
DIS(38,1:3) = ["EXP3", DEXP3.w * DEXP3.in.mdot / MEG, 0]; 
DIS(39,1:3) = ["COND", 0, CYC.Qrej / MEG];
DIS(40,1:3) = ["CMP1", DCMP1.w * DCMP1.in.mdot / MEG, 0]; 
DIS(41,1:3) = ["CMP2", DCMP2.w * DCMP2.in.mdot / MEG, 0]; 
DIS(42,1:3) = ["CMP3", DCMP3.w * DCMP3.in.mdot / MEG, 0]; 
DIS(43,1:3) = ["PREH", 0 , PREH.qc * PREH.Cout.mdot / MEG]; 
DIS(44,1:3) = ["BOILER", 0 , BOILER.qc * BOILER.Cout.mdot / MEG];
DIS(45,1:3) = ["SUPH", 0 , SUPH.qc * SUPH.Cout.mdot / MEG]; 
DIS(46,1:3) = ["REH", 0 , REH.qc * REH.Cout.mdot / MEG];  
%% TS Diagram
% retrieve data
S           = str2double(DIS([2:2:14,17,19,22:2:30],6));            % Extracting T data but skipping the hot inlets of FWH's and overlapped points (points 14, 19)
T           = str2double(DIS([2:2:14,17,19,22:2:30],3));            % Extracting S data but skipping the hot inlets of FWH's and overlapped points (points 14, 19)
H           = str2double(DIS([2:2:14,17,19,22:2:30],5));            % Extracting S data but skipping the hot inlets of FWH's and overlapped points (points 14, 19)

FWH1_T      = linspace(str2double(DIS(16,3)), str2double(DIS(15,3)), 50)  + degC;         % Isobaric line from hot inlet to outlet of FWH1
FWH2_T      = linspace(str2double(DIS(21,3)), str2double(DIS(20,3)), 50)  + degC;         % Isobaric line from hot inlet to outlet of FWH2
FWH1_S      = zeros(1,length(FWH1_T)); 
FWH2_S      = zeros(1,length(FWH2_T));
for i = 1:50     
    FWH1_S(i)   = CoolProp.PropsSI('S','P',FWH1.Cin.P,'T',FWH1_T(i),data.name) / KIL;   % Isobaric line from hot inlet to outlet of FWH1
    FWH2_S(i)   = CoolProp.PropsSI('S','P',FWH2.Cin.P,'T',FWH2_T(i),data.name) / KIL;   % Isobaric line from hot inlet to outlet of FWH2
end
P_sat_plot  = 0.1:0.1:220.6;
T_sat_L     = zeros(1,length(P_sat_plot));
T_sat_V     = zeros(1,length(P_sat_plot));
S_sat_L     = zeros(1,length(P_sat_plot));
S_sat_V     = zeros(1,length(P_sat_plot));
for i = 1:length(P_sat_plot)
    T_sat_L(i)      = CoolProp.PropsSI('T','P',P_sat_plot(i)*1e5,'Q',0,data.name);
    T_sat_V(i)     	= CoolProp.PropsSI('T','P',P_sat_plot(i)*1e5,'Q',1,data.name);
    S_sat_L(i)      = CoolProp.PropsSI('S','P',P_sat_plot(i)*1e5,'T',T_sat_L(i)-0.01,data.name) / KIL;
    S_sat_V(i)      = CoolProp.PropsSI('S','P',P_sat_plot(i)*1e5,'T',T_sat_V(i)+0.01,data.name) / KIL;
end
%{
% plot
figure(1)
fontsize = 14;
hold on
plot(S,T,'r-o','Linewidth',1.5)
% point labels if wanted
for i = [1:5 13]
    text(S(i)+0.1,T(i),num2str(i),'FontSize',fontsize)
end
text(S(6)-0.5,T(6)+15,'6,7','Fontsize',fontsize)
text(S(8)-0.75,T(8)+10,'8,9','Fontsize',fontsize)
text(S(10)-1,T(10)+10,'10,11','Fontsize',fontsize)
text(S(12)-0.5,T(12)+10,'12','Fontsize',fontsize)
plot(FWH1_S,FWH1_T - degC,'r-','Linewidth',1.5)
plot(FWH2_S,FWH2_T - degC,'r-','Linewidth',1.5)
plot(S_sat_L,T_sat_L - degC,'k-','Linewidth',1.5)
plot(S_sat_V,T_sat_V - degC,'k-','Linewidth',1.5)
ylabel('Temperature, C','fontweight','bold','fontsize',fontsize)
xlabel('Specific Entropy, kJ/kg/K','fontweight','bold','fontsize',fontsize)
set(gca,'fontweight','bold','fontsize',fontsize);
grid()

%}