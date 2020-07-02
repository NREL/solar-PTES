% Class for motors-generators
classdef gen_class
    properties
        type      % "gen" or "motor"
        
        pow
                
        % Economics
        gen_cost = econ_class(0,0,0,0) ;
        
    end
    
    methods
        function obj = gen_class(type, cost_mode)
            
            obj.type     = type ;
            obj.gen_cost = econ_class(cost_mode, 0.2, 5, 0.2) ;
            
        end
        
        function obj = gen_power(obj,CMP,EXP)
            
            WIN = 0.; WOUT = 0.;
            
            for ii = 1 : length(CMP)
                WIN = WIN + CMP(ii).W0;
            end
            for ii = 1 : length(EXP)
                WOUT = WOUT + EXP(ii).W0;
            end
            
            obj.pow  = abs(WIN - WOUT) ;
            
        end
               
        % Calculate the economic cost of the compressor or expander
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
                    COST = 1.85e6 * (obj.pow / 1.18e7)^0.94 ; % Really not sure where this came from ...
                    
                case 2
                    % Morandin eq. MG2
                    COST = 5e3 + 110 * obj.pow / 1e3 ;
                    COST = COST * CEind(curr) / CEind(2009) ;
                    
                case 3 
                    % Eq. MG3
                    COST = 40 * (obj.pow)^0.67 ;
                    COST = COST * CEind(curr) / CEind(2017) ;
                    
            end
                                    
            obj.gen_cost.COST = COST ;
        end
        
        
    end
end