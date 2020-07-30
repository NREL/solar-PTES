% Class for motors-generators
classdef gen_class
    properties
        type      % "gen" or "motor"
        
        pow0      % nominal power (W)
        
        pow_mech  % mechanical power during each load period (W)
        
        pow_el    % electrical power during each load period (W)
        
        E         % Electricity in/out during each load period (J)
        
        WL        % Work lost during each load period (J)
        
        part_load % describes part_load characteristics
        
        % Economics
        gen_cost = econ_class(0,0,0,0) ;
        
    end
    
    methods
        function obj = gen_class(type, cost_mode)
            
            obj.type     = type ;
            obj.gen_cost = econ_class(cost_mode, 0.2, 5, 0.2) ;
            
            % PART-LOAD PERFORMANCE
            % Data is from NASA2005 paper:
            % https://doi.org/10.1109/IEMDC.2005.195790
            % Data points were extracted from Figure 9 of the paper using a
            % separate Matlab script (see eta_PPN_extract.m inside the
            % ./Other/ folder).
            % PPN is an array of fractions of nominal power. eta is the
            % corresponding array of electrico-mechanical efficiencies.
            PPN = [...
                0.1873
                0.3722
                0.5631
                0.7520
                0.9350
                1.1159
                1.2909
                1.4699
                1.6409
                1.8060
                1.9651];
            eta = [...
                93.6448
                96.5847
                97.5027
                97.9290
                98.1475
                98.2240
                98.2787
                98.2787
                98.2240
                98.2022
                98.1257
                ]/100;
            obj.part_load.PPN = PPN;
            obj.part_load.eta = eta;
            
            % Plot part-load characteristics
            %{
            figure(10)
            x = linspace(0.05,2.5,100)';
            y = interp1(PPN,eta,x,'spline');
            plot(PPN,eta,'s',x,y); hold on;
            PPN_Patnode = linspace(0.05,1.0,100)';
            eta_Patnode = 0.9 + 0.258*PPN_Patnode - 0.3*PPN_Patnode.^2 ...
                + 0.12*PPN_Patnode.^3;
            plot(PPN_Patnode,eta_Patnode,'-.'); hold off;
            xlim([0 2.2])
            ylim([0.92 1.0])
            xlabel('Part load, P/PN')
            ylabel('Motor-generator efficiency [$\%$]')
            legend('NASA 2005','Interpolation','Padnode 2006','Location','East')
            grid on;
            keyboard
            %}
        end
        
        function obj = gen_power(obj,CCMP,CEXP,DCMP,DEXP,T)
            
            % Compute nominal motor power (W0 values in compressor and
            % expanders always defined as positive and in J/s)
            WIN0 = 0.; WOUT0 = 0.;
            for ii = 1 : length(CCMP)
                WIN0 = WIN0 + CCMP(ii).W0;
            end
            for ii = 1 : length(CEXP)
                WOUT0 = WOUT0 + CEXP(ii).W0;
            end
            MOT_pow0 = abs(WIN0 - WOUT0);
            
            % Compute nominal generator power
            WIN0 = 0.; WOUT0 = 0.;
            for ii = 1 : length(DCMP)
                WIN0 = WIN0 + DCMP(ii).W0;
            end
            for ii = 1 : length(DEXP)
                WOUT0 = WOUT0 + DEXP(ii).W0;
            end
            GEN_pow0 = abs(WIN0 - WOUT0);
            
            % Set nominal power for the motor-generator couple as the
            % maximum of the two
            obj.pow0  = max([MOT_pow0,GEN_pow0]);
            
            
            % Compute actual motor-generator power being used during each
            % load period (W arrays are defined as positive for expanders,
            % negative for compressors and in units of energy, J)
            WIN  = zeros(size(T));
            WOUT = zeros(size(T));
            for ii = 1 : length(CCMP)
                WIN  = WIN  + CCMP(ii).W;
            end
            for ii = 1 : length(DCMP)
                WIN  = WIN  + DCMP(ii).W;
            end
            for ii = 1 : length(CEXP)
                WOUT = WOUT + CEXP(ii).W;
            end
            for ii = 1 : length(DEXP)
                WOUT = WOUT + DEXP(ii).W;
            end
            MOT_GEN_pow = (WIN + WOUT)./T;
            
            % Obtain the part-load fraction for each load period
            PPN = abs(MOT_GEN_pow)/obj.pow0;
            
            % Obtain the electro-mechanical efficiency for each load period
            eta = interp1(obj.part_load.PPN,obj.part_load.eta,PPN,'makima');
            
            % Print warnings
            if any(PPN<0)
                error(['Motor-generator cannot operate under negative',...
                    'part-load conditions.'])
            end
            if any(0<PPN & PPN<0.15)
                warning(['Motor-generator operating in very low part-load',...
                    ' conditions. Predicted performance might be innacurate.'])
            end
            if any(PPN>1.2)
                warning(['Motor-generator operating above nominal power.'...
                    'Such operation might not be sustainable.'])
            end
            
            chg = MOT_GEN_pow < 0;
            dis = MOT_GEN_pow > 0;
            
            MOT_GEN_pow_el = zeros(size(MOT_GEN_pow));
            MOT_GEN_pow_el(chg) = MOT_GEN_pow(chg)./eta(chg);
            MOT_GEN_pow_el(dis) = MOT_GEN_pow(dis).*eta(dis);
            
            MOT_GEN_E  = MOT_GEN_pow_el.*T;
            MOT_GEN_WL = (MOT_GEN_pow - MOT_GEN_pow_el).*T;
            
            obj.pow_mech = MOT_GEN_pow;
            obj.pow_el   = MOT_GEN_pow_el;
            obj.E        = MOT_GEN_E;
            obj.WL       = MOT_GEN_WL;
        end
               
        % Calculate the economic cost of the motor-generator
        function obj = gen_econ(obj, CEind)
            mode = obj.gen_cost.cost_mode ;
            curr = 2019 ; % Current year
            
            % Many, many correlations exist
            % If an equation number is given, then that corresponds to an
            % equation given in the solar-PTES Q2 report
            switch mode
                case 0
                    COST = 0.01 ;
                    
                case 1
                    COST = 1.85e6 * (obj.pow0 / 1.18e7)^0.94 ; % Really not sure where this came from ...
                    
                case 2
                    % Morandin eq. MG2
                    COST = 5e3 + 110 * obj.pow0 / 1e3 ;
                    COST = COST * CEind(curr) / CEind(2009) ;
                    
                case 3 
                    % Eq. MG3
                    COST = 40 * (obj.pow0)^0.67 ;
                    COST = COST * CEind(curr) / CEind(2017) ;
                    
            end
                                    
            obj.gen_cost.COST = COST ;
        end
        
        
    end
end