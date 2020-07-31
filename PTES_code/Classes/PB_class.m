classdef PB_class
   properties
       Length
       Diameter
       
       hot_time
       cold_time
       
       Nu_hotstream
       Nu_coldstream
       CM
       NTU_hotstream
       
           
       eff     % Effectiveness
       pressure_drop_hot    % Pressure drop
       pressure_drop_cold
       total_cost
   end
   methods
       function obj = PB_class(Length, Diameter, hot_time, cold_time, Nu_hotstream, Nu_coldstream, CM, NTU_hotstream,eff_re,pressure_drop_hot,pressure_drop_cold,total_cost)
           obj.Length = [];
           obj.Diameter   = [];
           obj.hot_time = [];
           obj.cold_time = [];
           obj.Nu_hotstream = [];
           obj.Nu_coldstream = [];
           obj.CM=[];
           obj.NTU_hotstream = [];
           obj.eff = [];
           obj.pressure_drop_hot = [];
           obj.pressure_drop_cold = [];
           obj.total_cost = [];
       end
   end
end