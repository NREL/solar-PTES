function [mdot,reh_Tin,pre_Tin] = Steam_fxn(P_cond)
% This function runs the Rankine steam cycle at variable condenser
% pressure. The input file defines certain global variables pertinent to
% cycle design.
global BAR ...
    data FLD mdot_des mdot_tot DEN DE DCN HP_sat LP_sat HP_Tin IP_Tin  ...
    DE_eta DC_eta eff_ex_mode TTD xf_des dP ...
    DEXP1  DEXP2 DEXP3 DEXP4 DCOND DCMP1 DCMP2 DCMP3 DCMP4 FWH1 FWH2 ...
    PREH BOILER SUPH REH 
        
% Set up the temperature and pressure profiles along the steam cycle
HT_sat      = CoolProp.PropsSI('T','P',HP_sat*BAR,'Q',1,data.name);     %high saturation temperature (K)
LT_sat      = CoolProp.PropsSI('T','P',LP_sat*BAR,'Q',1,data.name);     %low saturation temperature (K)
T_sat       = zeros(1,DEN+1);                                           %construct linear profile for saturation temperature 
P_sat       = zeros(1,DEN+1);                                           %construct array for design saturation pressure
T_ex        = zeros(1,DEN+1);                                           %construct array for design exit temperatures

% Define design flow rates
m_refs      = mdot_des * [1 (1 - xf_des(1)) (1 - xf_des(1)) * (1 - xf_des(2))];
xf_guess    = xf_des;
err_xf      = 1; 
while err_xf > 1e-3    %While loop to find extractions. The cycle needs to asssume extractions to that all systems upstream of FWH's can be modeled
    % Guess flow rates
    m_guess     = mdot_tot * [1 (1 - xf_guess(1)) (1 - xf_guess(1)) * (1 - xf_guess(2))];

% Construct Expanders
    DE	= struct('Pin',cell(1,DEN),'Tin',cell(1,DEN),'beta',cell(1,DEN),'eta',cell(1,DEN),'type',cell(1,DEN));
    for i = 1:DEN+1
        T_sat(i)    = HT_sat - (HT_sat - LT_sat)*((i-1)/DEN);                   % set up temperature profile at design conditions
        P_sat(i)    = CoolProp.PropsSI('P','T', T_sat(i),'Q',1,data.name);    % "" pressure profile at design conditions
        T_ex(i)     = T_sat(i);
        % The last T,P coincides with condensor inlet; the inlet temperature of HP (and IP) expanders is superheated
    end 
    % Override some variables
    P_sat(DEN + 1)          = P_cond * BAR; % redefine cond. P/LP expander outlet
    T_sat(DEN + 1)          = CoolProp.PropsSI('T','P', P_sat(DEN + 1),'Q',1,data.name); % redefine LP expander outlet T
    T_ex(DEN + 1)           = T_sat(DEN + 1);
    T_ex(1)                 = HP_Tin;           % override HP expander temperature to superheat, C
    if DEN > 2
        T_ex(2:DEN-1)       = IP_Tin;          	% override IP expander temperature to reheat, C
    end    
    for i = 1:DEN
        DE(i).beta  = P_sat(i)/P_sat(i+1);      % define pressure ratio of each turbine based upon saturation profile
        DE(i).Pin   = P_sat(i);                 % turbine inlet pressures are Psat(1-3)
        DE(i).Tin   = T_ex(i);
        DE(i).eta   = eta_exp(DE_eta,m_refs(i),m_guess(i),eff_ex_mode); % expander effiicnecy updated if flowrate is off-design
        DE(i).type  = "EXP";
        switch i
            case 1
                DEXP1 = MACHINE(DE(i)) ;
            case 2 
                DEXP2 = MACHINE(DE(i)) ;
            case 3
                DEXP3 = MACHINE(DE(i)) ;
            case 4
                DEXP4 = MACHINE(DE(i)) ;
            case 5
               error('Code not written for case where there are > 5 discharging expanders');
        end 
    end

% Construct Compressors (Pumps)
    DC	= struct('Pin',cell(1,DCN),'Tin',cell(1,DCN),'beta',cell(1,DCN),'eta',cell(1,DCN),'type',cell(1,DCN));
    for i = 1:DCN
        % these conditions of the compressors are overridden later in the code
        % but these are design conditions
        DC(i).beta  = DE(DCN+1-i).beta;                 % pressure ratio of compressors is opposite order of expanders
        DC(i).Pin   = P_sat(DCN+2-i);                   % compressors are in opposite order as expanders
        DC(i).Tin   = T_sat(DCN+2-i) - TTD;             % TTD ensures no steam enters compressors 
        DC(i).eta   = eta_cmp(DC_eta,m_refs(DCN+1-i),m_guess(DCN+1-i)); % compressor effiicnecy updated if flowrate is off-design
        DC(i).type  = "CMP";
        switch i
            case 1
                DCMP1 = MACHINE(DC(i)) ;
            case 2 
                DCMP2 = MACHINE(DC(i)) ;
            case 3
                DCMP3 = MACHINE(DC(i)) ;
            case 4
                DCMP4 = MACHINE(DC(i)) ;
            case 5
                error('Code not written for case where there are > 5 discharging compressors');
        end
    end

% Rankine cycle - start at HP expander inlet
    % HP Expander
    DEXP1.in        = flow_state(DEXP1.in,'PT',FLD);
    DEXP1.in.mdot   = m_guess(1);
    DEXP1.out.mdot  = DEXP1.in.mdot;
    DEXP1           = calc_machine(DEXP1,FLD);

    % HP Expander --> Reheater
    hx_dat.type     = "Reheater" ;
    REH          = HEATEX(hx_dat);
    REH.Cdp      = dP; 
    REH.Cin      = DEXP1.out;
    REH.Cin.mdot = m_guess(2);
    REH.Cout.mdot= REH.Cin.mdot;
    REH.Cout.T   = DEXP2.in.T;
    REH          = calc_hx(REH,"cold",FLD);

    % Reheater --> IP Expander    
    DEXP2.in        = REH.Cout;
    DEXP2.out.mdot  = DEXP2.in.mdot;
    DEXP2           = calc_machine(DEXP2,FLD);
    
    % IP Expander --> LP Expander
    DEXP3.in        = DEXP2.out;
    DEXP3.in.mdot   = m_guess(3);
    DEXP3.out.mdot  = DEXP3.in.mdot;
    DEXP3           = calc_machine(DEXP3,FLD);
    
    % LP Expander --> Condenser
    hx_dat.type      = "Condenser";
    DCOND            = HEATEX(hx_dat);
    DCOND.Hin        = DEXP3.out;
    DCOND.Hdp        = dP;
    DCOND.Hout.T     = CoolProp.PropsSI('T','P',DCOND.Hin.P*(1 - DCOND.Hdp),'Q',0,data.name) - TTD;  % ensures exit stream is subcooled at exit pressure    
    DCOND.Hout.mdot  = DCOND.Hin.mdot;
    DCOND            = calc_hx(DCOND,"hot",FLD) ;
    
    % Condenser --> CMP1
    DCMP1.in        = DCOND.Hout;
    DCMP1.beta      = DEXP3.in.P/DCMP1.in.P; % redefines pressure ratio to account for pressure drop in condesner
    DCMP1.out.mdot  = DCMP1.in.mdot;
    DCMP1           = calc_machine(DCMP1,FLD);
   
    % IP Expander/CMP1 --> FWH1
    fwh_dat.type    = "OFWH";
    FWH1            = HEATEX(fwh_dat);
    FWH1.Hin        = DEXP2.out;
    FWH1.Hin.mdot   = m_guess(2) - m_guess(3);
    FWH1.Cin        = DCMP1.out;
    FWH1.Fout.T     = CoolProp.PropsSI('T','P',FWH1.Cin.P,'Q',0,data.name) - TTD; % exit steam is subcooled TTD below saturation
    FWH1.Fout.mdot  = FWH1.Hin.mdot + FWH1.Cin.mdot;
    FWH1            = calc_hx(FWH1,"mix",FLD);
   
    % FWH1 --> CMP2
    DCMP2.in        = FWH1.Fout;
    DCMP2.beta      = DEXP2.in.P/DCMP2.in.P; % redefines pressure ratio to account for pressure drop in condesner
    DCMP2.out.mdot  = DCMP2.in.mdot;
    DCMP2           = calc_machine(DCMP2,FLD);

    % HP Expander/CMP2 --> FWH2
    fwh_dat.type    = "OFWH";
    FWH2            = HEATEX(fwh_dat);
    FWH2.Hin        = DEXP1.out;
    FWH2.Hin.mdot   = m_guess(1) - m_guess(2);
    FWH2.Cin        = DCMP2.out;
    FWH2.Fout.T     = CoolProp.PropsSI('T','P',FWH2.Cin.P,'Q',0,data.name) - TTD;% exit steam is subcooled TTD below saturation
    FWH2.Fout.mdot  = FWH2.Hin.mdot + FWH2.Cin.mdot;
    FWH2            = calc_hx(FWH2,"mix",FLD);

    % FWH2 --> CMP3
    DCMP3.in        = FWH2.Fout;
    DCMP3.beta      = DEXP1.in.P/(1 - dP)^2/DCMP3.in.P; %compresses beyond EXP1 design to account for dP in heaters.
    DCMP3.out.mdot  = DCMP3.in.mdot;
    DCMP3           = calc_machine(DCMP3,FLD);

    % HP --> Solar preheat
    hx_dat.type     = "Preheater" ;
    PREH            = HEATEX(hx_dat);
    PREH.Cdp        = dP; 
    PREH.Cin        = DCMP3.out; 
    PREH.Cout.T     = CoolProp.PropsSI('T','P',DCMP3.out.P*(1 - PREH.Cdp),'Q',0,data.name) - 0.001; % (-0.001 K)saturation temperature at exit pressure of preheater
    PREH.Cout.mdot  = PREH.Cin.mdot;
    PREH            = calc_hx(PREH,"cold",FLD);

    % Solar preheat --> Steam Generation
    hx_dat.type     = "Boiler" ;
    BOILER          = HEATEX(hx_dat);
    BOILER.Cdp      = 0; % assume 0 pressure loss since the steam is not really flowing in a boiler...
    BOILER.Cin      = PREH.Cout;
    BOILER.Cout.T   = CoolProp.PropsSI('T','P',BOILER.Cin.P,'Q',1,data.name) + 0.001; % (+0.001 K) ensures steam undergoes phase change at pressure conditions
    BOILER.Cout.mdot = BOILER.Cin.mdot;
    BOILER          = calc_hx(BOILER,"cold",FLD);

    % Solar superheat
    hx_dat.type     = "Superheater" ;
    SUPH            = HEATEX(hx_dat);
    SUPH.Cdp        = dP; 
    SUPH.Cin        = BOILER.Cout;
    SUPH.Cout.T     = DEXP1.in.T;
    SUPH.Cout.mdot  = SUPH.Cin.mdot;
    SUPH            = calc_hx(SUPH,"cold",FLD);

    %overide extraction guesses
    xf_new          = [FWH2.x_f FWH1.x_f];
    err_xf          = norm(xf_new - xf_guess);
    xf_guess        = xf_new;
end
% outputs
    mdot            = mdot_tot * [1 (1 - xf_guess(1)) (1 - xf_guess(1)) * (1 - xf_guess(2))];
    reh_Tin         = DEXP1.out.T;
    pre_Tin         = DCMP1.out.T;

end
%% Off-Design Steam Functions
% expander (turbine) efficiency 
function eta_OD = eta_exp(eta_ref,m_ref,m_dot,mode)
% Patnode
    if mode == 1
        eta_red     =  1 - (0.191 - 0.409*(m_dot/m_ref) + 0.218*(m_dot/m_ref)^2);
%Atlas Copco
    elseif mode == 2
        eta_red     = (53.145*(m_dot/m_ref)^3 - 214.73*(m_dot/m_ref)^2 + 267.18*(m_dot/m_ref) - 5.3216) / 100;
    else %patnode by default
        eta_red     = 1 - (0.191 - 0.409*(m_dot/m_ref) + 0.218*(m_dot/m_ref)^2);
    end
    eta_OD          = eta_ref * eta_red;
end
% compressor (pump) efficiency 
function eta_OD = eta_cmp(eta_ref,m_ref,m_dot)
    eta_OD      = (2 * (m_dot/m_ref) - (m_dot/m_ref)^2) * eta_ref;
end