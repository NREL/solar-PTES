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
        ins_volA   = 0;
        ins_volB   = 0;
        
        % Costs
        costdat % A structure containing information for tank costing
        tankA_cost = econ_class(0,0,0,0) ; % Containment cost
        tankB_cost = econ_class(0,0,0,0) ; % Containment cost
        insA_cost  = econ_class(0,0,0,0) ; % Insulation cost
        insB_cost  = econ_class(0,0,0,0) ; % Insulation cost
        fluid_cost = econ_class(1,0,0,0) ; % Fluid cost
        
    end
    methods
         function obj = double_tank_class(fluid,TA,pA,MA,TB,pB,MB,T0,costdat,num)
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
             
             obj.costdat = costdat ;
             
             obj.tankA_cost = econ_class(costdat.tankmode, 0.2, 5, 0.2) ;
             obj.tankB_cost = econ_class(costdat.tankmode, 0.2, 5, 0.2) ;
             obj.insA_cost  = econ_class(1, 0.2, 5, 0.2) ;
             obj.insB_cost  = econ_class(1, 0.2, 5, 0.2) ;
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
             tank_state   = update_state(tank_state,obj,mode);             
             tank_state.V = tank_state.M/tank_state.rho;
             tank_state.H = tank_state.M*tank_state.h;
             tank_state.S = tank_state.M*tank_state.s;
             tank_state.B = tank_state.H - T0*tank_state.S;             
         end
         
         function [obj] = run_tanks(obj,iL,fluid,i_out,i_in,Load,T0)
             % RUN_TANKS Compute effect of various fluid streams entering
             % and leaving one double storage tank
             %
             %      USAGE: RUN_TANKS(obj, iL, fluid, i_out, i_in, Load, T0)
             %
             %      e.g.  run_tanks(HT, iL, fluidH, [1,3], [2,4], Load, T0)
             %
             %      OBJ is an instance of the double_tank_class.
             %
             %      FLUID is an instance of the fluid_class.
             %
             %      IL is the load index. I_OUT and I_IN are arrays
             %      containing the indices of the streams leaving the
             %      source tank and entering the sink tank, respectively.
             %
             %      This method computes the total mass, enthalpy and
             %      entropy flows of the fluid streams and adds/substracts
             %      them to/from the sink/source tanks.
             
             
             if any([i_out==0, i_in==0,isempty(i_in),isempty(i_out)])
                 warning('run_tanks might be using "empty" fluid streams')
                 obj.A(iL+1) = obj.A(iL);
                 obj.B(iL+1) = obj.B(iL);
                 return
             end
             
             % Compute flow rates into and out of the tanks
             Mdot_in  = 0; % Mass      flow rate into sink tank
             Hdot_in  = 0; % Enthalpy  flow rate into sink tank
             Sdot_in  = 0; % Entropy   flow rate into sink tank (before mixing)
             Mdot_out = 0; % Mass      flow rate out of source tank
             Hdot_out = 0; % Enthalpy  flow rate out of source tank
             Sdot_out = 0; % Entropy   flow rate out of source tank
             for i = i_out
                 Mdot_out = Mdot_out + fluid.state(iL,i).mdot;
                 Hdot_out = Hdot_out + fluid.state(iL,i).h.*fluid.state(iL,i).mdot;
                 Sdot_out = Sdot_out + fluid.state(iL,i).s.*fluid.state(iL,i).mdot;
             end
             for i = i_in
                 Mdot_in  = Mdot_in  + fluid.state(iL,i).mdot;
                 Hdot_in  = Hdot_in  + fluid.state(iL,i).h.*fluid.state(iL,i).mdot;
                 Sdot_in  = Sdot_in  + fluid.state(iL,i).s.*fluid.state(iL,i).mdot;
             end
             if abs(Mdot_out/Mdot_in - 1) > 1e-6
                 error('mass flow rates do not match!')
             else
                 Mdot = Mdot_out;
             end
             
             % Set conditions of combined stream entering sink tank
             mix_state   = fluid.state(iL,i_in(1));
             mix_state.h = Hdot_in/Mdot;
             mix_state.p = fluid.state(iL,i_in(1)).p;
             mix_state.mdot = Mdot;
             mix_state   = update_state(mix_state,fluid,2);
             s_mix = mix_state.s; % entropy flow into sink tank (after mixing)
             
             t = Load.time(iL);
             
             switch obj.job
                 case 'SF'
                     % Select sink tank and source tank depending on operation mode
                     if any(strcmp(Load.type(iL),{'chg','chgCO2','chgTSCO2'}))
                         SO1   = obj.A(iL);   %source tank
                         SO2   = obj.A(iL+1); %source tank
                         SI1   = obj.B(iL);   %sink tank
                         SI2   = obj.B(iL+1); %sink tank
                     elseif any(strcmp(Load.type(iL),{'dis','ran','disCO2','rcmpCO2','disTSCO2'}))
                         SI1   = obj.A(iL);   %sink tank
                         SI2   = obj.A(iL+1); %sink tank
                         SO1   = obj.B(iL);   %source tank
                         SO2   = obj.B(iL+1); %source tank
                     end
                     
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
                     if any(strcmp(Load.type(iL),{'chg','chgCO2','chgTSCO2'}))
                         obj.A(iL+1) = SO2;
                         obj.B(iL+1) = SI2;
                     elseif any(strcmp(Load.type(iL),{'dis','ran','disCO2','rcmpCO2','disTSCO2'}))
                         obj.A(iL+1) = SI2;
                         obj.B(iL+1) = SO2;
                     end
                     
                 case 'ENV'
                     
                     % Atmospheric tanks remain unmodified, as they are
                     % considered to be infinite
                     obj.A(iL+1) = obj.A(iL);
                     obj.B(iL+1) = obj.B(iL);
                     
                     % Compute entropy generation due to bringing the warm
                     % air streams (which absorbed rejected heat) back to
                     % ambient conditions: sirr = Ds - Dh/T0
                     S_irr = (Sdot_out -  Sdot_in - (Hdot_out - Hdot_in)/T0)*t;
                     
                 otherwise
                     error('not implemented')
             end
             
             if any(strcmp(Load.type(iL),{'chg','chgCO2','chgTSCO2'}))
                 obj.WL_chg = obj.WL_chg + T0*S_irr;
             elseif any(strcmp(Load.type(iL),{'dis','ran','disCO2','rcmpCO2','disTSCO2'}))
                 obj.WL_dis = obj.WL_dis + T0*S_irr;
             end
             
             if any([obj.WL_chg<0,obj.WL_dis<0])
                 %keyboard
             end
         end
         
         % Calculate some stats for the tank including max volume of fluid,
         % max. mass of fluid
         % This is necessary because if there are consecutive charge cycles
         % it may not be clear what the total change in fluid mass is
         function [obj] = tank_stats(obj)
             
             min_volA = 1e21;
             min_masA = 1e21;
             
             max_volA = 0;
             max_masA = 0;
             
             min_volB = 1e21;
             min_masB = 1e21;
             
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
             obj.tank_volA = obj.fluid_volA * obj.costdat.over_fac ;
             obj.tank_volB = obj.fluid_volB * obj.costdat.over_fac ;
             
             % TANK A: length, diameter, surface area
             AD = (4. * obj.tank_volA / (pi * obj.costdat.AR)) ^ (1./3.) ;
             AR = 0.5 * AD ;
             AL = AD * obj.costdat.AR ;
             AA = pi * AD * (AL + AD/4.) ;
             AA = 8 * AA / 5;
             
             CF = RP1('PT_INPUTS',obj.A(1).p,obj.A(1).T,'CPMASS',obj) ;
             tau = obj.costdat.tau * 24 * 3600 ; % convert from days to seconds
             UA = obj.tank_volA * obj.A(1).rho * CF / (tau * AA) ; % Overall heat transfer coef
             
             % Insulation thickness and volume
             tinsA = AR * (exp(obj.costdat.ins_k / (AR * UA)) - 1) ;
             obj.ins_volA = pi * tinsA * (2 + tinsA) * AL ; % Side wall insulation
             obj.ins_volA = obj.ins_volA + 2 * pi * tinsA * (AR + tinsA)^2 ; % Top and bottom insulation
             
             % New tank length, diameter, and volume
             AL = AL + 2 * tinsA ;
             AD = AD + 2 * tinsA ;
             obj.tank_volA = 0.25 * pi * AL * AD^2 ;
             
             % TANK B: length, diameter, surface area
             BD = (4. * obj.tank_volB / (pi * obj.costdat.AR)) ^ (1./3.) ;
             BR = BD * 0.5 ;
             BL = BD * obj.costdat.AR ;
             BA = pi * BD * (BL + BD/4.) ;
             
             CF = RP1('PT_INPUTS',obj.B(1).p,obj.B(1).T,'CPMASS',obj) ;
             tau = obj.costdat.tau * 24 * 3600 ; % convert from days to seconds
             UB = obj.tank_volB * obj.B(1).rho * CF / (tau * BA) ; % Overall heat transfer coef
             
             % Insulation thickness and volume
             tinsB = BR * (exp(obj.costdat.ins_k / (BR * UB)) - 1) ;
             obj.ins_volB = pi * tinsB * (2 + tinsB) * BL ; % Side wall insulation
             obj.ins_volB = obj.ins_volB + 2 * pi * tinsB * (BR + tinsB)^2 ; % Top and bottom insulation
             
             % New tank length, diameter, and volume
             BL = BL + 2 * tinsB ;
             BD = BD + 2 * tinsB ;
             obj.tank_volB = 0.25 * pi * BL * BD^2 ;
             
             
         end
         
         
         % Calculate the cost of the tank
         function [obj] = tank_cost(obj, CEind)
             
             curr = 2019 ; % Current year
             mode = obj.tankA_cost.cost_mode ;
             
             % Numerous cost correlations available
             switch mode
                 case 0 
                     
                     Acost = 0.01 ;
                     Bcost = 0.01 ;
                 
                 case 1
                     % Tank cost based upon Peters + Timmerhaus
                     % See solar-PTES Q2 report, equation PV2 - updated
                     % with CEindex for 2019
                     
                     % Maximum tank volume is 50e3 m3. 
                     maxV  = 50e3 ;
                     
                     if obj.tank_volA > maxV
                         nA    = floor(obj.tank_volA / maxV); % Number of maxV tanks
                         vA    = mod(obj.tank_volA / maxV,1)  ; % Remaining volume
                         Acost = nA * 3829 * maxV^0.557 ;
                         
                         Acost = Acost + 3829 * vA^0.557 ;
                     else
                         Acost = 3829 * obj.tank_volA^0.557 ;
                     end
                     
                     if obj.tank_volB > maxV
                         nB    = floor(obj.tank_volB / maxV); % Number of maxV tanks
                         vB    = mod(obj.tank_volB / maxV,1)  ; % Remaining volume
                         Bcost = nB * 3829 * maxV^0.557 ;
                         
                         Bcost = Bcost + 3829 * vB^0.557 ;
                     else
                         Bcost = 3829 * obj.tank_volB^0.557 ;
                     end
                     
                     % Increase cost if pressurized
                     p = obj.A(1).p / 1e5 ;
                     if p > 7
                         Cfact = 0.922 + 0.0335*p - 0.0003*p^2 +1e-6*p^3 ;
                         Acost = Acost * Cfact ;
                         Bcost = Bcost * Cfact ;
                     end
                     
                     Acost = Acost * CEind(curr) / CEind(1990) ;
                     Bcost = Bcost * CEind(curr) / CEind(1990) ;
                     
                 case 2
                     % NETL, carbon steel fixed roof, eq. PV1.2
                     Acost = 40288 + 50.9 * obj.tank_volA ;
                     Bcost = 40288 + 50.9 * obj.tank_volB ;
                     
                     Acost = Acost * CEind(curr) / CEind(1998) ;
                     Bcost = Bcost * CEind(curr) / CEind(1998) ;
                     
                 case 3
                     % FIxed roof storage tank, eq. PV5
                     Acost = 8e3 + 600 * obj.tank_volA^0.78 ;
                     Bcost = 8e3 + 600 * obj.tank_volB^0.78 ;
                     
                     Acost = Acost * CEind(curr) / CEind(2009) ;
                     Bcost = Bcost * CEind(curr) / CEind(2009) ;
                     
                 case 4 
                     % Tank up to 10 bar, eq. PV 6
                     Acost = 10e3 + 4.5e3 * obj.tank_volA^0.6 ;
                     Bcost = 10e3 + 4.5e3 * obj.tank_volB^0.6 ;
                     
                     Acost = Acost * CEind(curr) / CEind(1998) ;
                     Bcost = Bcost * CEind(curr) / CEind(1998) ;
                     
                 case 5
                     % Tank cost based upon Peters + Timmerhaus
                     % But updated to match Jack plant data
                     
                     % Maximum tank volume is 50e3 m3. 
                     maxV  = 50e3 ;
                     
                     if obj.tank_volA > maxV
                         nA    = floor(obj.tank_volA / maxV); % Number of maxV tanks
                         vA    = mod(obj.tank_volA / maxV,1)  ; % Remaining volume
                         Acost = nA * 13673 * maxV^0.557 ;
                         
                         Acost = Acost + 13673 * vA^0.557 ;
                     else
                         Acost = 13673 * obj.tank_volA^0.557 ;
                     end
                     
                     if obj.tank_volB > maxV
                         nB    = floor(obj.tank_volB / maxV); % Number of maxV tanks
                         vB    = mod(obj.tank_volB / maxV,1)  ; % Remaining volume
                         Bcost = nB * 13673 * maxV^0.557 ;
                         
                         Bcost = Bcost + 13673 * vB^0.557 ;
                     else
                         Bcost = 13673 * obj.tank_volB^0.557 ;
                     end
                     
                     % Increase cost if pressurized
                     p = obj.A(1).p / 1e5 ;
                     if p > 7
                         Cfact = 0.922 + 0.0335*p - 0.0003*p^2 +1e-6*p^3 ;
                         Acost = Acost * Cfact ;
                         Bcost = Bcost * Cfact ;
                     end
                     
                     Acost = Acost * CEind(curr) / CEind(2015) ;
                     Bcost = Bcost * CEind(curr) / CEind(2015) ;
                     
                 case 6
                     % Storage tank cost based on Couper
                     % Up to 11 million gallons for field constructed tanks
                     % These costs seem ridiculously high! Looks like there
                     % is a typo in the book 1218 --> 1.218
                     
                     % Maximum tank volume is 50e3 m3. 
                     m3toGAL = 264.172 ;
                     maxV  = 50e3 * m3toGAL ;
                     CF = 2.7 ; % Cost factor for Stainless Steel 316
                     
                     volA = obj.tank_volA * m3toGAL ;
                     volB = obj.tank_volB * m3toGAL ;
                     
                     if volA > maxV
                         nA    = floor(volA / maxV); % Number of maxV tanks
                         vA    = mod(volA / maxV,1)  ; % Remaining volume
                         Acost = nA * CF * 1.218 * exp(11.662-0.6104 * log(maxV) + 0.04536 * (log(maxV)^2)) ;
                         
                         Acost = Acost + CF * 1.218 * exp(11.662-0.6104 * log(vA) + 0.04536 * (log(vA)^2)) ;
                     else
                         Acost = CF * 1.218 * exp(11.662-0.6104 * log(volA) + 0.04536 * (log(volA)^2)) ;
                     end
                     
                     if volB > maxV
                         nB    = floor(volB / maxV); % Number of maxV tanks
                         vB    = mod(volB / maxV,1)  ; % Remaining volume
                         Bcost = nB * CF * 1.218 * exp(11.662-0.6104 * log(maxV) + 0.04536 * (log(maxV)^2)) ;
                         
                         Bcost = Bcost + CF * 1.218 * exp(11.662-0.6104 * log(vB) + 0.04536 * (log(vB)^2)) ;
                     else
                         Bcost = CF * 1.218 * exp(11.662-0.6104 * log(volB) + 0.04536 * (log(volB)^2)) ;
                     end
                     
                     % Increase cost if pressurized
                     p = obj.A(1).p / 1e5 ;
                     if p > 7
                         Cfact = 0.922 + 0.0335*p - 0.0003*p^2 +1e-6*p^3 ;
                         Acost = Acost * Cfact ;
                         Bcost = Bcost * Cfact ;
                     end
                     
                     Acost = Acost * CEind(curr) / CEind(2010) ;
                     Bcost = Bcost * CEind(curr) / CEind(2010) ;
                     
                 case 10
                     
                     error('Mode not available for tank costs')
             end
             
             obj.tankA_cost.COST = Acost ;
             obj.tankB_cost.COST = Bcost ;
                
         end
         
         % Calculate the cost of the fluid. cost_kg is the cost per kg of fluid
         function [obj] = fld_cost(obj, CEind)
             
             curr = 2019 ;
             cost_kg = obj.costdat.fld_cost ;
             COST = obj.fluid_mass * cost_kg ;
             COST = COST * CEind(curr) / CEind(2019) ;
             
             obj.fluid_cost.COST = COST ;
             
         end
         
         % Calculate the cost of the insulation. cost_kg is the cost per kg of fluid
         function [obj] = ins_cost(obj, CEind)
             
             curr = 2019 ;
             cost_kg = obj.costdat.ins_cost ;
             
             COST = obj.ins_volA * obj.costdat.ins_rho * cost_kg ;
             COST = COST * CEind(curr) / CEind(2019) ;
             obj.insA_cost.COST = COST ;
             
             COST = obj.ins_volB * obj.costdat.ins_rho * cost_kg ;
             COST = COST * CEind(curr) / CEind(2019) ;
             obj.insB_cost.COST = COST ;

         end
         
         
    end
end