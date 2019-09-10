classdef FWH 
    % FWH is a construct for open feedwater heaters for steam Rankine cycles\
    
    properties
        name    % name
        note    % any other details
        Hin     % FLOW class - hot inlet flow properties
        Cin     % FLOW class - cold inlet flow properties
        Fout    % FLOW class - mixed outlet flow properties
        x_f     % steam extraction from turbine (%)
        Hdp     % Hot side pressure drop
        Cdp     % Cold side pressure drop, given as a percentage of inlet
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
        type    % 'OFWH' for open FWH
    end   
    
    %Constructor function 
    methods
        function obj = FWH(fwh_dat)
            if nargin == 1      %number of input arguments
                
                if isfield(fwh_dat, 'type')
                    obj.type    = fwh_dat.type;
                end
                
                % Set up flow classes
                flow_dat.nobv = 1 ;
                obj.Hin     = FLOW(flow_dat);
                obj.Cin     = FLOW(flow_dat);
                obj.Fout    = FLOW(flow_dat);
            end
        end
    end
    
    % other functions
    methods
        
        %calculate the energy balance on the FWH
        function obj = calc_fwh(obj,fld)
        
            % mass balance on FWH, essentially a mixing tank
            obj.Fout.P  = obj.Hin.P;
            obj.Fout    = flow_state(obj.Fout,'PT',fld);
            obj.x_f     = (obj.Fout.h - obj.Cin.h)/(obj.Hin.h - obj.Cin.h);     %steam extraction from turbines
%             obj.Fout.h  = obj.x_f * obj.Hin.h + (1 - obj.x_f) * obj.Cin.h;
%             if obj.x_f >=0 && obj.x_f <=1
                
%             else 
%                 error('Extraction did not converge');
%             end
        end
    end
end
        
    