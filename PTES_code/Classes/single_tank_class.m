classdef single_tank_class
    properties
        name   % substance (storage medium)
        job    % 'WF' (working fluid) or 'SM' (storage media)
        read   % 'CP' (CoolProp) or 'TAB' (table)
        handle % integer to identify CoolProp AbstractState
        HEOS   % integer to identify CoolProp's HEOS AbstractState (if necessary)
        TAB    % matrix of thermophysical properties
        A = tank_state_class; %source/sink single tank state
        DH_chg % enthalpy change during charge
        DH_str % enthalpy change during storage
        DH_dis % enthalpy change during discharge
        WL_chg % exergetic loss during charge
        WL_str % exergetic loss during storage
        WL_dis % exergetic loss during discharge
        
        % Tank volumes
        fluid_mass = 0;
        fluid_volA = 0;
        tank_volA  = 0;
        ins_volA   = 0;
        
        % Costs
        costdat % A structure containing information for tank costing
        tankA_cost = econ_class(0,0,0,0) ; % Containment cost
        insA_cost  = econ_class(0,0,0,0) ; % Insulation cost
        fluid_cost = econ_class(1,0,0,0) ; % Fluid cost
        
    end
    methods
         function obj = single_tank_class(fluid,TA,pA,MA,T0,costdat,num)
             obj.name   = fluid.name;
             obj.job    = fluid.job;
             obj.read   = fluid.read;
             obj.handle = fluid.handle;
             obj.HEOS   = fluid.HEOS;
             obj.TAB    = fluid.TAB;
             obj.A(1:num) = tank_state_class;
             obj.A(1).T = TA;
             obj.A(1).p = pA;
             obj.A(1).M = MA;
             obj.A(1)   = update_tank_state(obj,obj.A(1),T0,1);
             obj.DH_chg = 0;
             obj.DH_str = 0;
             obj.DH_dis = 0;
             obj.WL_chg = 0;
             obj.WL_str = 0;
             obj.WL_dis = 0;
             
             obj.costdat = costdat ;
             
             obj.tankA_cost = econ_class(costdat.tankmode(1), 0.2, 5, 0.2) ;
             obj.insA_cost  = econ_class(costdat.ins_cost(1), 0.2, 5, 0.2) ;
             obj.fluid_cost = econ_class(costdat.fld_cost(1), 0.2, 5, 0.2) ;
             
         end
         
         function obj = reset_single_tank(obj,TA,pA,MA,T0)
             obj.A(:) = tank_state_class;
             obj.A(1).T = TA;
             obj.A(1).p = pA;
             obj.A(1).M = MA;
             obj.A(1)   = update_tank_state(obj,obj.A(1),T0,1);
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
         
         function [obj] = run_single_tank(obj,iL,fluid,i_out,i_in,Load,T0)
             % RUN_TANKS Compute effect of various fluid streams entering
             % and leaving one single storage tank
             %
             %      USAGE: RUN_SINGLE_TANK(obj, iL, fluid, i_out, i_in, Load, T0)
             %
             %      e.g.  run_single_tank(HT, iL, fluidH, [1,3], [2,4], Load, T0)
             %
             %      OBJ is an instance of the single_tank_class.
             %
             %      FLUID is an instance of the fluid_class.
             %
             %      IL is the load index. I_OUT and I_IN are arrays
             %      containing the indices of the streams leaving and
             %      entering the tank, respectively.
             %
             %      This method computes the total mass, enthalpy and
             %      entropy flows of the fluid streams and adds/substracts
             %      them to/from the tank.
             
             % Check that pressures are the same for all incoming and
             % outcoming fluid streams
             i_tot = [i_out,i_in];
             p = fluid.state(iL,i_tot(1)).p;
             if any(abs([fluid.state(iL,i_tot).p]/p-1)>1e-6)
                 error('Inconsistent pressure streams')
             end
             obj.A(iL).p=p;
             
             if isempty(i_in)
                 i_in = 0;
             end
             if isempty(i_out)
                 i_out = 0;
             end
             
             if all([i_out==0, i_in==0])
                 warning('run_single_tank might be using "empty" fluid streams')
                 obj.A(iL+1) = obj.A(iL);
                 return
             end
             
             % Compute flow rates into and out of the tank
             Mdot_in  = 0; % Mass      flow rate into tank
             Hdot_in  = 0; % Enthalpy  flow rate into tank
             Sdot_in  = 0; % Entropy   flow rate into tank (before mixing)
             Mdot_out = 0; % Mass      flow rate out of source tank
             Hdot_out = 0; % Enthalpy  flow rate out of source tank
             Sdot_out = 0; % Entropy   flow rate out of source tank
             if i_out~=0
                 for i = i_out
                     Mdot_out = Mdot_out + fluid.state(iL,i).mdot;
                     Hdot_out = Hdot_out + fluid.state(iL,i).h.*fluid.state(iL,i).mdot;
                     Sdot_out = Sdot_out + fluid.state(iL,i).s.*fluid.state(iL,i).mdot;
                 end
             end
             if i_in~=0
                 for i = i_in
                     Mdot_in  = Mdot_in  + fluid.state(iL,i).mdot;
                     Hdot_in  = Hdot_in  + fluid.state(iL,i).h.*fluid.state(iL,i).mdot;
                     Sdot_in  = Sdot_in  + fluid.state(iL,i).s.*fluid.state(iL,i).mdot;
                 end
             end
             
             t = Load.time(iL);
             
             % Compute change in enthalpy of the tank
             DH = (Hdot_in - Hdot_out)*t;
             
             switch obj.job
                 case {'SF','WF'}
                     % Select initial and ending states of the tank
                     A1 = obj.A(iL);
                     A2 = obj.A(iL+1);
                     
                     % Compute end conditions of the tank
                     A2.M = A1.M - Mdot_out*t + Mdot_in*t;
                     A2.H = A1.H - Hdot_out*t + Hdot_in*t;
                     A2.h = A2.H/A2.M;
                     A2.p = A1.p;
                     A2   = update_tank_state(obj,A2,T0,2);
                     
                     % Compute entropy generation of mixing
                     S_irr = A2.S - A1.S + Sdot_out*t - Sdot_in*t;
                     
                     % Set new ending state of the tank
                     obj.A(iL+1) = A2;
                     
                 case 'ENV'
                     % Select initial and ending states of the tank
                     A1 = obj.A(iL);
                     A2 = obj.A(iL+1);
                     
                     % Compute end conditions of the tank
                     A2.M = A1.M - Mdot_out*t + Mdot_in*t;
                     A2.H = A1.H - Hdot_out*t + Hdot_in*t;
                     A2.h = A2.H/A2.M;
                     A2.p = A1.p;
                     A2   = update_tank_state(obj,A2,T0,2);
                     
                     % Compute entropy generation of mixing
                     S_irr = A2.S - A1.S + Sdot_out*t - Sdot_in*t;
                     
                     % Now compute entropy generation of bringing tank back
                     % to T0
                     A3   = A2;
                     A3.T = T0;
                     A3   = update_tank_state(obj,A3,T0,1);
                     S_irr2 = (A3.S -  A2.S - (A3.H - A2.H)/T0);
                     S_irr = S_irr+S_irr2;
                     
                     % Set new ending state of the tank
                     obj.A(iL+1) = A3;
                     
                 otherwise
                     error('not implemented')
             end
             
             % Compute and save lost work due to entropy generation
             if any(strcmp(Load.type(iL),{'chg','chgCO2','chgTSCO2','chgCC','chgICC','chgICC_PC'}))
                 obj.DH_chg = obj.DH_chg + DH;
                 obj.WL_chg = obj.WL_chg + T0*S_irr;
             elseif any(strcmp(Load.type(iL),{'dis','ran','disCO2','rcmpCO2','disTSCO2','disCC','disICC','disICC_PC'}))
                 obj.DH_dis = obj.DH_dis + DH;
                 obj.WL_dis = obj.WL_dis + T0*S_irr;
             else
                 error('Unrecognised Load.type operation mode')
             end
             
         end
         
         % Calculate some stats for the tank including max volume of fluid,
         % max. mass of fluid
         % This is necessary because if there are consecutive charge cycles
         % it may not be clear what the total change in fluid mass is
         function [obj] = single_tank_stats(obj)
             
             min_volA = 1e21;
             min_masA = 1e21;
             
             max_volA = 0;
             max_masA = 0;
             
             N = numel(obj.A) ;
             
             for i = 1 : N
                
                 min_volA = min(obj.A(i).V,min_volA);
                 
                 min_masA = min(obj.A(i).M,min_masA);
                 
                 max_volA = max(obj.A(i).V,max_volA);
                 
                 max_masA = max(obj.A(i).M,max_masA);
                 
             end
             
             
             obj.fluid_mass = max_masA - min_masA ;
             obj.fluid_volA = max_volA - min_volA ;
             
             % Add on a bit of extra volume to calculate the tank volume
             obj.tank_volA = obj.fluid_volA * obj.costdat.over_fac ;
             
             % TANK A: length, diameter, surface area
             AD = (4. * obj.tank_volA / (pi * obj.costdat.AR)) ^ (1./3.) ;
             AR = 0.5 * AD ;
             AL = AD * obj.costdat.AR ;
             AA = pi * AD * (AL + AD/4.) ;
             AA = 8 * AA / 5;
             
             %CF = RPN('PT_INPUTS',obj.A(1).p,obj.A(1).T,'CPMASS',obj) ;
             CF = RPN('HmassP_INPUTS',obj.A(1).h,obj.A(1).p,'CPMASS',obj) ;
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
             
         end
         
         
         % Calculate the cost of the tank
         function [obj] = single_tank_cost(obj, CEind)
             
             curr = 2019 ; % Current year
             mode = obj.tankA_cost.cost_mode ;
             
             % Numerous cost correlations available
             switch mode
                 case 0 
                     
                     Acost = 0.01 ;
                 
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
                     
                     % Increase cost if pressurized
                     p = obj.A(1).p / 1e5 ;
                     if p > 7
                         Cfact = 0.922 + 0.0335*p - 0.0003*p^2 +1e-6*p^3 ;
                         Acost = Acost * Cfact ;
                     end
                     
                     Acost = Acost * CEind(curr) / CEind(1990) ;
                     
                 case 2
                     % Maximum tank volume is 40e3 m3. 
                     maxV = 40e3 ;
                     VA   = obj.tank_volA ;
                     nA   = 1;
                     if VA > maxV
                         nA = VA / maxV; % Number of maxV tanks
                         VA = maxV ;
                     end
                    
                     % NETL, carbon steel fixed roof, eq. PV1.2
                     Acost = nA * (40288 + 50.9 * VA) ;
                     
                     Acost = Acost * CEind(curr) / CEind(1998) ;
                     
                 case 3
                     % FIxed roof storage tank, eq. PV5
                     Acost = 8e3 + 600 * obj.tank_volA^0.78 ;
                     
                     Acost = Acost * CEind(curr) / CEind(2009) ;
                     
                 case 4 
                     % Tank up to 10 bar, eq. PV 6
                     Acost = 10e3 + 4.5e3 * obj.tank_volA^0.6 ;
                     
                     Acost = Acost * CEind(curr) / CEind(1998) ;
                     
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
                     
                     % Increase cost if pressurized
                     p = obj.A(1).p / 1e5 ;
                     if p > 7
                         Cfact = 0.922 + 0.0335*p - 0.0003*p^2 +1e-6*p^3 ;
                         Acost = Acost * Cfact ;
                     end
                     
                     Acost = Acost * CEind(curr) / CEind(2015) ;
                     
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
                     
                     if volA > maxV
                         nA    = floor(volA / maxV); % Number of maxV tanks
                         vA    = mod(volA / maxV,1)  ; % Remaining volume
                         Acost = nA * CF * 1.218 * exp(11.662-0.6104 * log(maxV) + 0.04536 * (log(maxV)^2)) ;
                         
                         Acost = Acost + CF * 1.218 * exp(11.662-0.6104 * log(vA) + 0.04536 * (log(vA)^2)) ;
                     else
                         Acost = CF * 1.218 * exp(11.662-0.6104 * log(volA) + 0.04536 * (log(volA)^2)) ;
                     end
                     
                     % Increase cost if pressurized
                     p = obj.A(1).p / 1e5 ;
                     if p > 7
                         Cfact = 0.922 + 0.0335*p - 0.0003*p^2 +1e-6*p^3 ;
                         Acost = Acost * Cfact ;
                     end
                     
                     Acost = Acost * CEind(curr) / CEind(2010) ;
                     
                 case 10
                     
                     error('Mode not available for tank costs')
             end
             
             obj.tankA_cost.COST = Acost ;
                
         end
         
         % Calculate the cost of the fluid. cost_kg is the cost per kg of fluid
         function [obj] = fld_cost(obj, CEind)
             
             curr = 2019 ;
             cost_kg = obj.fluid_cost.cost_mode ;
             COST = obj.fluid_mass * cost_kg ;
             COST = COST * CEind(curr) / CEind(2019) ;
             
             obj.fluid_cost.COST = COST ;
             
         end
         
         % Calculate the cost of the insulation. cost_kg is the cost per kg of fluid
         function [obj] = ins_cost(obj, CEind)
             
             curr = 2019 ;
             cost_kg = obj.insA_cost.cost_mode ;
             
             COST = obj.ins_volA * obj.costdat.ins_rho * cost_kg ;
             COST = COST * CEind(curr) / CEind(2019) ;
             obj.insA_cost.COST = COST ;
             
             cost_kg = obj.insB_cost.cost_mode ;
             COST = obj.ins_volB * obj.costdat.ins_rho * cost_kg ;
             COST = COST * CEind(curr) / CEind(2019) ;
             obj.insB_cost.COST = COST ;

         end
         
         
    end
end