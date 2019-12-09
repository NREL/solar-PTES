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
        
        % Tank volumes
        fluid_mass = 0;
        fluid_volA = 0;
        fluid_volB = 0;
        tank_volA  = 0;
        tank_volB  = 0;
        
        % Costs
        tankA_cost = econ_class(0,0,0,0) ;
        tankB_cost = econ_class(0,0,0,0) ;
        fluid_cost = econ_class(0,0,0,0) ;
        
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
             
             obj.tankA_cost = econ_class(1, 0.2, 5, 0.2) ;
             obj.tankB_cost = econ_class(1, 0.2, 5, 0.2) ;
             obj.fluid_cost = econ_class(1, 0.2, 5, 0.2) ;
             
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
             elseif any(strcmp(Load.type(iL),{'dis','ran','disCO2'}))
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
             elseif any(strcmp(Load.type(iL),{'dis','ran','disCO2'}))
                 obj.A(iL+1) = SI2;
                 obj.B(iL+1) = SO2;
                 obj.WL_dis = obj.WL_dis + T0*S_irr;
             end
         end
         
         % Calculate some stats for the tank including max volume of fluid,
         % max. mass of fluid
         
         % This is necessary because if there are consecutive charge cycles
         % it may not be clear what the total change in fluid mass is
         function [obj] = tank_stats(obj)
             
             min_volA = 1e11;
             min_masA = 1e11;
             
             max_volA = 0;
             max_masA = 0;
             
             min_volB = 1e11;
             min_masB = 1e11;
             
             max_volB = 0;
             max_masB = 0;
             
             N = numel(obj.A) ;
             
             for i = 1 : N
                
                 min_volA = min(obj.A(i).V,min_volA);
                 min_volB = min(obj.B(i).V,min_volB);
                 
                 min_masA = min(obj.A(i).M,min_masA);
                 min_masB = min(obj.B(i).M,min_masB);
                 
                 max_volA = max(obj.A(i).V,max_volA);
                 max_volB = max(obj.B(i).V,max_volB);
                 
                 max_masA = max(obj.A(i).M,max_masA);
                 max_masB = max(obj.B(i).M,max_masB);
                 
             end
             
             
             obj.fluid_mass = max_masA - min_masA ;
             obj.fluid_volA = max_volA - min_volA ;
             obj.fluid_volB = max_volB - min_volB ;
             
             % Add on a bit of extra volume to calculate the tank volume
             fact = 1.1 ;
             obj.tank_volA = obj.fluid_volA * fact ;
             obj.tank_volB = obj.fluid_volB * fact ;
             
         end
         
         
         % Calculate the cost of the tank
         function [obj] = tank_cost(obj)
             
             mode = obj.tankA_cost.cost_mode ;
             % Numerous cost correlations available
             switch mode
                 case 1
                     % Tank cost based upon Peters + Timmerhaus
                     % See solar-PTES Q2 report, equation PV2 - updated
                     % with CEindex for 2019
                     Acost = 6629.4 * obj.tank_volA^0.557 ;
                     Bcost = 6629.4 * obj.tank_volB^0.557 ;
                     
                     % Increase cost if pressurized
                     p = obj.A(1).p / 1e5 ;
                     if p > 7
                         Cfact = 0.922 + 0.0335*p - 0.0003*p^2 +1e-6*p^3 ;
                         Acost = Acost * Cfact ;
                         Bcost = Bcost * Cfact ;
                     end
                     
                 case 2
                     error('Mode not available for tank costs')
             end
             
             obj.tankA_cost.COST = Acost ;
             obj.tankB_cost.COST = Bcost ;
                
         end
         
         % Calculate the cost of the fluid. cost_kg is the cost per kg of
         % fluid
         function [obj] = fld_cost(obj,cost_kg)
             
             obj.fluid_cost.COST = obj.fluid_mass * cost_kg ;
                
         end
         
         
    end
end