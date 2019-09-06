close all
%% Joule Cycle
CHG(1,1:9) = ["Point","Component","T, C","P, bar","h, kJ/kg","s, kJ/kg.K","q","Density, kg/m3", "Flow, m^3/s"];
CHG(2,1:9) = [1, CCMP1.type, CCMP1.in.T-degC, CCMP1.in.P/BAR, CCMP1.in.h/KIL, CCMP1.in.s/KIL,CCMP1.in.q, CCMP1.in.rho, -m_joule/CCMP1.in.rho] ;
CHG(3,1:9) = [2, CCMP1.type, CCMP1.out.T-degC, CCMP1.out.P/BAR, CCMP1.out.h/KIL, CCMP1.out.s/KIL,CCMP1.out.q, CCMP1.out.rho, -m_joule/CCMP1.out.rho] ;

CHG(4,1:9) = [2, cHS.type, cHS.Hin.T-degC, cHS.Hin.P/BAR, cHS.Hin.h/KIL, cHS.Hin.s/KIL, cHS.Hin.q, cHS.Hin.rho, -m_joule/cHS.Hin.rho] ;
CHG(5,1:9) = [3, cHS.type, cHS.Hout.T-degC, cHS.Hout.P/BAR, cHS.Hout.h/KIL, cHS.Hout.s/KIL, cHS.Hout.q, cHS.Hout.rho, -m_joule/cHS.Hout.rho] ;

CHG(6,1:9) = [3, cRCP.type, cRCP.Hin.T-degC, cRCP.Hin.P/BAR, cRCP.Hin.h/KIL, cRCP.Hin.s/KIL, cRCP.Hin.q, cRCP.Hin.rho, -m_joule/cRCP.Hin.rho] ;
CHG(7,1:9) = [4, cRCP.type, cRCP.Hout.T-degC, cRCP.Hout.P/BAR, cRCP.Hout.h/KIL, cRCP.Hout.s/KIL, cRCP.Hout.q, cRCP.Hout.rho,-m_joule/cRCP.Hout.rho] ;

CHG(8,1:9) = [4, CEXP1.type, CEXP1.in.T-degC, CEXP1.in.P/BAR, CEXP1.in.h/KIL, CEXP1.in.s/KIL, CEXP1.in.q, CEXP1.in.rho, -m_joule/CEXP1.in.rho] ;
CHG(9,1:9) = [5, CEXP1.type, CEXP1.out.T-degC, CEXP1.out.P/BAR, CEXP1.out.h/KIL, CEXP1.out.s/KIL, CEXP1.out.q, CEXP1.out.rho, -m_joule/CEXP1.out.rho] ;

CHG(10,1:9) = [5, cCS1.type, cCS1.Cin.T-degC, cCS1.Cin.P/BAR, cCS1.Cin.h/KIL, cCS1.Cin.s/KIL, cCS1.Cin.q, cCS1.Cin.rho, -m_joule/cCS1.Cin.rho] ;
CHG(11,1:9) = [6, cCS1.type, cCS1.Cout.T-degC, cCS1.Cout.P/BAR, cCS1.Cout.h/KIL, cCS1.Cout.s/KIL, cCS1.Cout.q, cCS1.Cout.rho, -m_joule/cCS1.Cout.rho] ;

CHG(12,1:9) = [7, CEXP2.type, CEXP2.in.T-degC, CEXP2.in.P/BAR, CEXP2.in.h/KIL, CEXP2.in.s/KIL, CEXP2.in.q, CEXP2.in.rho, -m_joule/CEXP2.in.rho] ;
CHG(13,1:9) = [8, CEXP2.type, CEXP2.out.T-degC, CEXP2.out.P/BAR, CEXP2.out.h/KIL, CEXP2.out.s/KIL, CEXP2.out.q, CEXP2.out.rho, -m_joule/CEXP2.out.rho] ;

CHG(14,1:9) = [9, cCS2.type, cCS2.Cin.T-degC, cCS2.Cin.P/BAR, cCS2.Cin.h/KIL, cCS2.Cin.s/KIL, cCS2.Cin.q, cCS2.Cin.rho, -m_joule/cCS2.Cin.rho] ;
CHG(15,1:9) = [10, cCS2.type, cCS2.Cout.T-degC, cCS2.Cout.P/BAR, cCS2.Cout.h/KIL, cCS2.Cout.s/KIL, cCS2.Cout.q, cCS2.Cout.rho, -m_joule/cCS2.Cout.rho] ;

CHG(16,1:9) = [7, CEXP3.type, CEXP3.in.T-degC, CEXP3.in.P/BAR, CEXP3.in.h/KIL, CEXP3.in.s/KIL, CEXP3.in.q, CEXP3.in.rho, -m_joule/CEXP3.in.rho] ;
CHG(17,1:9) = [8, CEXP3.type, CEXP3.out.T-degC, CEXP3.out.P/BAR, CEXP3.out.h/KIL, CEXP3.out.s/KIL, CEXP3.out.q, CEXP3.out.rho, -m_joule/CEXP3.out.rho] ;

CHG(18,1:9) = [9, cCS3.type, cCS3.Cin.T-degC, cCS3.Cin.P/BAR, cCS3.Cin.h/KIL, cCS3.Cin.s/KIL, cCS3.Cin.q, cCS3.Cin.rho, -m_joule/cCS3.Cin.rho] ;
CHG(19,1:9) = [10, cCS3.type, cCS3.Cout.T-degC, cCS3.Cout.P/BAR, cCS3.Cout.h/KIL, cCS3.Cout.s/KIL, cCS3.Cout.q, cCS3.Cout.rho, -m_joule/cCS3.Cout.rho] ;

CHG(20,1:9) = [11, cRCP.type, cRCP.Cin.T-degC, cRCP.Cin.P/BAR, cRCP.Cin.h/KIL, cRCP.Cin.s/KIL, cRCP.Cin.q, cRCP.Cin.rho, -m_joule/cRCP.Cin.rho] ;
CHG(21,1:9) = [12, cRCP.type, cRCP.Cout.T-degC, cRCP.Cout.P/BAR, cRCP.Cout.h/KIL, cRCP.Cout.s/KIL, cRCP.Cout.q, cRCP.Cout.rho, -m_joule/cRCP.Cout.rho] ;
%% HS
CHG(22,1:9) = ["Point","Hot Storage","T, C","P, bar","h, kJ/kg","s, kJ/kg.K","q","Density, kg/m3", "Flow, kg/s"];
CHG(23,1:9) = [3, cHS.type, cHS.Cin.T-degC, cHS.Cin.P/BAR, cHS.Cin.h/KIL, cHS.Hout.s/KIL, cHS.Cin.q, cHS.Cin.rho, cHS.Cin.mdot] ;
CHG(24,1:9) = [2, cHS.type, cHS.Cout.T-degC, cHS.Cout.P/BAR, cHS.Cout.h/KIL, cHS.Hin.s/KIL, cHS.Cout.q, cHS.Cout.rho, cHS.Cout.mdot] ;
%% CS1
CHG(25,1:9) = ["Point","Cold Storage","T, C","P, bar","h, kJ/kg","s, kJ/kg.K","q","Density, kg/m3", "Flow, kg/s"];
CHG(26,1:9) = [6, cCS1.type, cCS1.Hin.T-degC, cCS1.Hin.P/BAR, cCS1.Hin.h/KIL, cCS1.Cout.s/KIL, cCS1.Hin.q, cCS1.Hin.rho, cCS1.Hin.mdot] ;
CHG(27,1:9) = [5, cCS1.type, cCS1.Hout.T-degC, cCS1.Hout.P/BAR, cCS1.Hout.h/KIL, cCS1.Cin.s/KIL, cCS1.Hout.q, cCS1.Hout.rho, cCS1.Hout.mdot] ;
CHG(28,1:9) = [10, cCS2.type, cCS2.Hin.T-degC, cCS2.Hin.P/BAR, cCS2.Hin.h/KIL, cCS2.Cout.s/KIL, cCS2.Hin.q, cCS2.Hin.rho, cCS2.Hin.mdot] ;
CHG(29,1:9) = [9, cCS2.type, cCS2.Hout.T-degC, cCS2.Hout.P/BAR, cCS2.Hout.h/KIL, cCS2.Cin.s/KIL, cCS2.Hout.q, cCS2.Hout.rho, cCS2.Hout.mdot] ;
%% Steam cycle
DIS(1,1:9)  = ["Point","Component","T, K","P, bar","h, kJ/kg","s, kJ/kg.K","Mass flow, kg/s","q","Density, kg/m3"];
CYC_Rank.eta_net    = CYC_Rank.WnetD / CYC_Joule.WnetC *100;    % Net cycle efficiency 
% HP Expander
DIS(2,1:9)  = [1, "EXP1 in",    DEXP1.in.T - degC, DEXP1.in.P/BAR, DEXP1.in.h/KIL, DEXP1.in.s/KIL, DEXP1.in.mdot, DEXP1.in.q, DEXP1.in.rho] ;
DIS(3,1:9)  = [2, "EXP1 out",	DEXP1.out.T - degC, DEXP1.out.P/BAR, DEXP1.out.h/KIL, DEXP1.out.s/KIL, DEXP1.out.mdot, DEXP1.out.q, DEXP1.out.rho] ;
% Solar Reheat
DIS(4,1:9)  = [3, "SOL_REH in",	SOL_reh.Cin.T - degC, SOL_reh.Cin.P/BAR, SOL_reh.Cin.h/KIL, SOL_reh.Cin.s/KIL, SOL_reh.Cin.mdot, SOL_reh.Cin.q, SOL_reh.Cin.rho] ;
DIS(5,1:9)  = [4, "SOL_REH out",SOL_reh.Cout.T - degC, SOL_reh.Cout.P/BAR, SOL_reh.Cout.h/KIL, SOL_reh.Cout.s/KIL, SOL_reh.Cout.mdot, SOL_reh.Cout.q, SOL_reh.Cout.rho] ;
% IP Expander
DIS(6,1:9)  = [5, "EXP2 in",    DEXP2.in.T - degC, DEXP2.in.P/BAR, DEXP2.in.h/KIL, DEXP2.in.s/KIL, DEXP2.in.mdot, DEXP2.in.q, DEXP2.in.rho] ;
DIS(7,1:9)  = [6, "EXP2 out",   DEXP2.out.T - degC, DEXP2.out.P/BAR, DEXP2.out.h/KIL, DEXP2.out.s/KIL, DEXP2.out.mdot, DEXP2.out.q, DEXP2.out.rho] ;
% LP Expander
DIS(8,1:9)  = [7, "EXP3 in",    DEXP3.in.T - degC, DEXP3.in.P/BAR, DEXP3.in.h/KIL, DEXP3.in.s/KIL, DEXP3.in.mdot, DEXP3.in.q, DEXP3.in.rho] ;
DIS(9,1:9)  = [8, "EXP3 out",   DEXP3.out.T - degC, DEXP3.out.P/BAR, DEXP3.out.h/KIL, DEXP3.out.s/KIL, DEXP3.out.mdot, DEXP3.out.q, DEXP3.out.rho] ;
% Condenser
DIS(10,1:9) = [9,  "COND in",   dREJ.Hin.T - degC, dREJ.Hin.P/BAR, dREJ.Hin.h/KIL, dREJ.Hin.s/KIL, dREJ.Hin.mdot, dREJ.Hin.q, dREJ.Hin.rho] ;
DIS(11,1:9) = [10, "COND out",  dREJ.Hout.T - degC, dREJ.Hout.P/BAR, dREJ.Hout.h/KIL, dREJ.Hout.s/KIL, dREJ.Hout.mdot, dREJ.Hout.q, dREJ.Hout.rho] ;
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
DIS(24,1:9) = [23, "SOL_PRE in"	 SOL_pre.Cin.T - degC, SOL_pre.Cin.P/BAR, SOL_pre.Cin.h/KIL, SOL_pre.Cin.s/KIL, SOL_pre.Cin.mdot, SOL_pre.Cin.q, SOL_pre.Cin.rho] ;
DIS(25,1:9) = [24, "SOL_PRE out", SOL_pre.Cout.T - degC, SOL_pre.Cout.P/BAR, SOL_pre.Cout.h/KIL, SOL_pre.Cout.s/KIL, SOL_pre.Cout.mdot, SOL_pre.Cout.q, SOL_pre.Cout.rho] ;
% Solar generation
DIS(26,1:9) = [25, "SOL_GEN in", SOL_gen.Cin.T - degC, SOL_gen.Cin.P/BAR, SOL_gen.Cin.h/KIL, SOL_gen.Cin.s/KIL, SOL_gen.Cin.mdot, SOL_gen.Cin.q, SOL_gen.Cin.rho] ;
DIS(27,1:9) = [26, "SOL_GEN out",SOL_gen.Cout.T - degC, SOL_gen.Cout.P/BAR, SOL_gen.Cout.h/KIL, SOL_gen.Cout.s/KIL, SOL_gen.Cout.mdot, SOL_gen.Cout.q, SOL_gen.Cout.rho] ;
% Solar superheat
DIS(28,1:9) = [27, "SOL_SUP in", SOL_sup.Cin.T - degC, SOL_sup.Cin.P/BAR, SOL_sup.Cin.h/KIL, SOL_sup.Cin.s/KIL, SOL_sup.Cin.mdot, SOL_sup.Cin.q, SOL_sup.Cin.rho] ;
DIS(29,1:9) = [28, "SOL_SUP out",SOL_sup.Cout.T - degC, SOL_sup.Cout.P/BAR, SOL_sup.Cout.h/KIL, SOL_sup.Cout.s/KIL, SOL_sup.Cout.mdot, SOL_sup.Cout.q, SOL_sup.Cout.rho] ;
%Completed cycle
DIS(30,1:9) = [1, "EXP1 in",    DEXP1.in.T - degC, DEXP1.in.P/BAR, DEXP1.in.h/KIL, DEXP1.in.s/KIL, DEXP1.in.mdot, DEXP1.in.q, DEXP1.in.rho] ;

DIS(35,1:3) = ["Component" , "Work in/out (MW)", "Heat in (MW)"];
DIS(36,1:3) = ["EXP1", DEXP1.w * DEXP1.in.mdot / MEG, 0]; 
DIS(37,1:3) = ["EXP2", DEXP2.w * DEXP2.in.mdot / MEG, 0]; 
DIS(38,1:3) = ["EXP3", DEXP3.w * DEXP3.in.mdot / MEG, 0]; 
DIS(39,1:3) = ["COND", 0, CYC_Rank.Qrej / MEG];
DIS(40,1:3) = ["CMP1", DCMP1.w * DCMP1.in.mdot / MEG, 0]; 
DIS(41,1:3) = ["CMP2", DCMP2.w * DCMP2.in.mdot / MEG, 0]; 
DIS(42,1:3) = ["CMP3", DCMP3.w * DCMP3.in.mdot / MEG, 0]; 
DIS(43,1:3) = ["SOL_PRE", 0 , SOL_pre.qc * SOL_pre.Cout.mdot / MEG]; 
DIS(44,1:3) = ["SOL_GEN", 0 , SOL_gen.qc * SOL_gen.Cout.mdot / MEG];
DIS(45,1:3) = ["SOL_SUP", 0 , SOL_sup.qc * SOL_sup.Cout.mdot / MEG]; 
DIS(46,1:3) = ["SOL_reh", 0 , SOL_reh.qc * SOL_reh.Cout.mdot / MEG];  
%% Rankine Diagram
S           = str2double(DIS([2:2:14,17,19,22:2:30],6));            % Extracting T data but skipping the hot inlets of FWH's and overlapped points (points 14, 19)
T           = str2double(DIS([2:2:14,17,19,22:2:30],3));            % Extracting S data but skipping the hot inlets of FWH's and overlapped points (points 14, 19)
H           = str2double(DIS([2:2:14,17,19,22:2:30],5));            % Extracting S data but skipping the hot inlets of FWH's and overlapped points (points 14, 19)

% labels      = {'1','2','3','4','5','6','7','8','9','10','11','12','13',''};
FWH1_T      = linspace(str2double(DIS(16,3)), str2double(DIS(15,3)), 50)  + degC;         % Isobaric line from hot inlet to outlet of FWH1
FWH2_T      = linspace(str2double(DIS(21,3)), str2double(DIS(20,3)), 50)  + degC;         % Isobaric line from hot inlet to outlet of FWH2
FWH1_S      = zeros(1,length(FWH1_T)); 
FWH2_S      = zeros(1,length(FWH2_T));
for i = 1:50     
    FWH1_S(i)   = CoolProp.PropsSI('S','P',FWH1.Cin.P,'T',FWH1_T(i),rank_data.name) / KIL;   % Isobaric line from hot inlet to outlet of FWH1
    FWH2_S(i)   = CoolProp.PropsSI('S','P',FWH2.Cin.P,'T',FWH2_T(i),rank_data.name) / KIL;   % Isobaric line from hot inlet to outlet of FWH2
end
P_sat_plot  = 0.1:0.1:220.6;
T_sat_L     = zeros(1,length(P_sat_plot));
T_sat_V     = zeros(1,length(P_sat_plot));
S_sat_L     = zeros(1,length(P_sat_plot));
S_sat_V     = zeros(1,length(P_sat_plot));
for i = 1:length(P_sat_plot)
    T_sat_L(i)      = CoolProp.PropsSI('T','P',P_sat_plot(i)*1e5,'Q',0,rank_data.name);
    T_sat_V(i)     	= CoolProp.PropsSI('T','P',P_sat_plot(i)*1e5,'Q',1,rank_data.name);
    S_sat_L(i)      = CoolProp.PropsSI('S','P',P_sat_plot(i)*1e5,'T',T_sat_L(i)-0.01,rank_data.name) / KIL;
    S_sat_V(i)      = CoolProp.PropsSI('S','P',P_sat_plot(i)*1e5,'T',T_sat_V(i)+0.01,rank_data.name) / KIL;
end
%% Steam TS and HS diagrams
figure(1)
% subplot(2,1,1)
fontsize = 14;
hold on
plot(S,T,'r-o','Linewidth',1.5)
% for i = [1:5 13]
%     text(S(i)+0.1,T(i),num2str(i),'FontSize',14)
% end
% text(S(6)-0.5,T(6)+15,'6,7','Fontsize',14)
% text(S(8)-0.75,T(8)+10,'8,9','Fontsize',14)
% text(S(10)-1,T(10)+10,'10,11','Fontsize',14)
% text(S(12)-0.5,T(12)+10,'12','Fontsize',14)
plot(FWH1_S,FWH1_T - degC,'r-','Linewidth',1.5)
plot(FWH2_S,FWH2_T - degC,'r-','Linewidth',1.5)
plot(S_sat_L,T_sat_L - degC,'k-','Linewidth',1.5)
plot(S_sat_V,T_sat_V - degC,'k-','Linewidth',1.5)
ylabel('Temperature, C','fontweight','bold','fontsize',fontsize)
xlabel('Specific Entropy, kJ/kg/K','fontweight','bold','fontsize',fontsize)
set(gca,'fontweight','bold','fontsize',12);
grid()
%% Joule TS diagram
S_chg       = str2double(CHG([2:2:20,21],6));
T_chg       = str2double(CHG([2:2:20,21],3));

S_HS        = str2double(CHG(23:24,6));
T_HS        = str2double(CHG(23:24,3));

S_CS        = str2double(CHG(26:27,6));
T_CS        = str2double(CHG(26:27,3));

S_CS2     	= str2double(CHG(28:29,6));
T_CS2     	= str2double(CHG(28:29,3));
figure(2)
hold on
plot(S_chg,T_chg,'r-x','Linewidth',1.5)
ylim([-100 600])
for i = [1 2 10]
    text(S_chg(i) + 0.01, T_chg(i),num2str(i),'FontSize',10)
end
for i = [3 4]
    text(S_chg(i) - 0.025, T_chg(i) + 5,num2str(i),'FontSize',10)
end
for i = [5 7 9]
    text(S_chg(i), T_chg(i) - 10,num2str(i),'FontSize',10)
end
for i = [6 8]
    text(S_chg(i), T_chg(i) + 10,num2str(i),'FontSize',10)
end
plot(S_HS,T_HS, 'b:','Linewidth',1.5)
plot(S_CS,T_CS, 'b:','Linewidth',1.5)
plot(S_CS2,T_CS2, 'b:','Linewidth',1.5)
ylabel('Temperature, �C','fontweight','bold','fontsize',fontsize)
xlabel('Specific Entropy, kJ/kgK','fontweight','bold','fontsize',fontsize)
set(gca,'fontweight','bold','fontsize',12);
grid()


