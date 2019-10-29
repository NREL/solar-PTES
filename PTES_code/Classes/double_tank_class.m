classdef double_tank_class
    properties
        name   % substance (storage medium)
        job    % 'WF' (working fluid) or 'SM' (storage media)
        read   % 'CP' (CoolProp) or 'TAB' (table)
        handle % integer to identify CoolProp AbstractState
        TAB    % matrix of thermophysical properties
        A = tank_state_class; % source tank (full at start of charge cycle)
        B = tank_state_class; % sink tank (full at start of discharge cycle)
        WL_chg % exergetic loss during charge
        WL_str % exergetic loss during storage
        WL_dis % exergetic loss during discharge
    end
    methods        
         function obj = double_tank_class(fluid,TA,pA,MA,TB,pB,MB,T0,num)
             obj.name   = fluid.name;
             obj.job    = fluid.job;
             obj.read   = fluid.read;
             obj.handle = fluid.handle;
             obj.TAB    = fluid.TAB;
             obj.A(1:num) = tank_state_class;
             obj.B(1:num) = tank_state_class;
             obj.A(1).T = TA;
             obj.A(1).p = pA;
             obj.A(1).M = MA;
             obj.A(1)   = update_tank_state(obj,obj.A(1),T0,1);
             obj.B(1).T = TB;
             obj.B(1).p = pB;
             obj.B(1).M = MB;
             obj.B(1)   = update_tank_state(obj,obj.B(1),T0,1);
             obj.WL_chg = 0;
             obj.WL_str = 0;
             obj.WL_dis = 0;
         end
         
         function obj = reset_tanks(obj,TA,pA,MA,TB,pB,MB,T0)
             obj.A(:) = tank_state_class;
             obj.B(:) = tank_state_class;
             obj.A(1).T = TA;
             obj.A(1).p = pA;
             obj.A(1).M = MA;
             obj.A(1)   = update_tank_state(obj,obj.A(1),T0,1);
             obj.B(1).T = TB;
             obj.B(1).p = pB;
             obj.B(1).M = MB;
             obj.B(1)   = update_tank_state(obj,obj.B(1),T0,1);
             obj.WL_chg = 0;
             obj.WL_str = 0;
             obj.WL_dis = 0;
         end
         
         function tank_state = update_tank_state(obj, tank_state, T0, mode)
             tank_state   = update_state(tank_state,obj.handle,obj.read,obj.TAB,mode);             
             tank_state.V = tank_state.M/tank_state.rho;
             tank_state.H = tank_state.M*tank_state.h;
             tank_state.S = tank_state.M*tank_state.s;
             tank_state.B = tank_state.H - T0*tank_state.S;             
         end
         
         function [obj] = run_tanks(obj,fluid_streams,iL,Load,T0)
             % Compute total mass, enthalpy and entropy flows of the fluid_streams
             Mdot = 0;
             Hdot_in  = 0; % Enthalpy flow into sink tank
             Hdot_out = 0; % Enthalpy flow out of source tank
             Sdot_in  = 0; % Entropy  flow into sink tank
             n = numel(fluid_streams);
             for i=1:n
                 Mdot     = Mdot     + fluid_streams(i).state(iL,2).mdot;
                 Hdot_in  = Hdot_in  + fluid_streams(i).state(iL,2).h*fluid_streams(i).state(iL,2).mdot;
                 Hdot_out = Hdot_out + fluid_streams(i).state(iL,1).h*fluid_streams(i).state(iL,1).mdot;
                 Sdot_in  = Sdot_in  + fluid_streams(i).state(iL,2).s*fluid_streams(i).state(iL,2).mdot;
             end
             
             % Set conditions of combined stream entering sink tank
             mix_stream   = fluid_streams(1).state(iL,2);             
             mix_stream.h = Hdot_in/Mdot;
             mix_stream.p = fluid_streams(1).state(iL,2).p;
             mix_stream.mdot = Mdot;
             mix_stream   = update_state(mix_stream,fluid_streams(1).handle,fluid_streams(1).read,fluid_streams(1).TAB,2);
             s_mix = mix_stream.s;
             
             % Select sink tank and source tank depending on operation mode
             if any(strcmp(Load.type(iL),{'chg','chgCO2'}))
                 SO1   = obj.A(iL);   %source tank
                 SO2   = obj.A(iL+1); %source tank
                 SI1   = obj.B(iL);   %sink tank
                 SI2   = obj.B(iL+1); %sink tank
             elseif any(strcmp(Load.type(iL),{'dis','disCO2'}))
                 SI1   = obj.A(iL);   %sink tank
                 SI2   = obj.A(iL+1); %sink tank
                 SO1   = obj.B(iL);   %source tank
                 SO2   = obj.B(iL+1); %source tank
             end
             t = Load.time(iL);
             
             % Update end conditions of source tank A
             SO2.M = SO1.M - Mdot*t;
             SO2.H = SO1.H - Hdot_out*t;
             SO2.h = SO2.H/SO2.M;
             SO2.p = SO1.p;
             SO2   = update_tank_state(obj,SO2,T0,2);
             
             % Update end conditions of sink tank B
             SI2.M = SI1.M + Mdot*t;
             SI2.H = SI1.H + Hdot_in*t;
             SI2.h = SI2.H/SI2.M;
             SI2.p = SI1.p;
             SI2   = update_tank_state(obj,SI2,T0,2);             
             
             % Compute entropy generation of mixing
             S_irr1 = (s_mix*Mdot - Sdot_in)*t; %mixing of streams before sink tank inlet
             S_irr2 = SI2.S - (SI1.S + s_mix*Mdot*t); %mixing of streams with fluid inside sink tank
             S_irr  = S_irr1 + S_irr2;
             if any(strcmp(Load.type(iL),{'chg','chgCO2'})) 
                 obj.A(iL+1) = SO2;
                 obj.B(iL+1) = SI2;
                 obj.WL_chg = obj.WL_chg + T0*S_irr;
             elseif any(strcmp(Load.type(iL),{'dis','disCO2'}))
                 obj.A(iL+1) = SI2;
                 obj.B(iL+1) = SO2;
                 obj.WL_dis = obj.WL_dis + T0*S_irr;
             end
         end
    end
end