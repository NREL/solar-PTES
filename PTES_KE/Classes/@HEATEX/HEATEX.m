classdef HEATEX 
    % HEATEX is a construct for heat exchangers and recuperators
    
    properties
        name    % name
        note    % Any other details
        Hin     % FLOW class - hot inlet flow properties
        Hout    % FLOW class - hot outlet flow properties
        Cin     % FLOW class - cold inlet flow properties
        Cout    % FLOW class - cold outlet flow properties
        Fout    % FLOW class - flow out of feedwater preheaters
        x_f     % Steam extraction from turbine (%) if feedwater heater
        Hdp     % Hot side pressure drop
        Cdp     % Cold side pressure drop, given as a percentage of inlet pressure
        dT      % Temperature difference between hot and cold sides
        qh      % Specific heat out of hot stream, kJ/kg
        qc      % Specific heat into cold stream, kJ/kg
        Qh      % Heat out of hot stream, kJ/kg
        Qc      % Heat into cold stream, kJ/kg
        Qmax    % Maximum heat that could be transferred, kJ/kg
        eff     % Effectiveness, Qh / Qmax = Qc / Qmax
        q       % Specific heat transfer into the system (positive), kJ/kg
        WLTH    % Work lost due to thermal irreversibilities, kJ/kg
        WLQL    % Work lost due to heat leakage, kJ/kg
        WLT     % Work lost total, kJ/kg   
        type    % 'RCP' for recuperator, 'INT' for intercooler, 'REJ' for heat rejection to environment, 'HS' for hot store, 'CS' for cold store, 'SOL' for solar heat addition, 'pre', 'sup' 'vap' for preheaters, superheaters and vaporizers, 'OFWH' for open feedwaterheater
    end
    
    % Constructor function
    methods
        function obj = HEATEX(hx_dat) 
           if nargin == 1 
                                            
               if isfield(hx_dat , 'type')
                   obj.type = hx_dat.type ;
               end
               
               if isfield(hx_dat , 'Hdp')
                   obj.Hdp = hx_dat.Hdp ;
               else
                   obj.Hdp = 0.0 ;
               end
               
               if isfield(hx_dat , 'Cdp')
                   obj.Cdp = hx_dat.Cdp ;
               else
                   obj.Cdp = 0.0 ;
               end
               
               % Set up flow classes
               flow_dat.nobv = 1 ;
               obj.Hin  = FLOW(flow_dat) ;
               obj.Hout = FLOW(flow_dat) ;
               
               obj.Cin  = FLOW(flow_dat) ;
               obj.Cout = FLOW(flow_dat) ;
               
               obj.Fout = FLOW(flow_dat);
           end
        end
    end
    
    % Other functions
    methods
        
        % Calculate heat transfer into/out of a stream
        function obj = calc_hx(obj,side,fld)
            if (obj.type == 'OFWH') 
                % in OFWH, stream properties in and out are known. Find
                % extraction from turbines needed.
                obj.Fout.P  = obj.Hin.P;
                obj.Fout    = flow_state(obj.Fout,'PT',fld);
                obj.x_f     = (obj.Fout.h - obj.Cin.h)/(obj.Hin.h - obj.Cin.h);
            else
            % calculate either hot side or cold side
                if (side == "hot")
                   obj.Hout.P   = obj.Hin.P * (1 - obj.Hdp) ;
                   obj.Hout     = flow_state(obj.Hout,'PT',fld) ;
                   obj.qh       = obj.Hout.h - obj.Hin.h ;
                elseif (side == "cold")
                   obj.Cout.P   = obj.Cin.P * (1 - obj.Cdp) ;
                   obj.Cout     = flow_state(obj.Cout,'PT',fld) ;
                   obj.qc       = obj.Cout.h - obj.Cin.h ;
                else
                   error('For heat exchangers must specify either a (hot) side or a (cold) side');
                end
            end
                        
        end
    end
    
end



        
    