% Class for miscellaneous components that don't require lots of info
classdef env_mix_class
    properties
        type      % e.g "mix"
        
        mdot    % Mass flow rate
        
        w       % specific work transfer, J/kg
        q       % specific heat transfer, J/kg
        Dh      % specific enthalpy change, J/kg
        sirr    % specific entropy generation, J/kgK
        
        W        % work transfer, W
        Q        % heat transfer, W
        DH       % enthalpy change, W
        Sirr     % entropy generation, W/K
        
    end
    
    methods
        function obj = env_mix_class(type, numPeriods)
            obj.type = type ;
            
            obj.mdot = zeros(numPeriods,1) ;
            
            obj.w    = zeros(numPeriods,1) ;
            obj.q    = zeros(numPeriods,1) ;
            obj.Dh   = zeros(numPeriods,1) ;
            obj.sirr = zeros(numPeriods,1) ;
                        
            obj.W    = zeros(numPeriods,1) ;
            obj.Q    = zeros(numPeriods,1) ;
            obj.DH   = zeros(numPeriods,1) ;
            obj.Sirr = zeros(numPeriods,1) ;
            
        end
        
        % Calculate the energy and entropy terms for the miscellaneous  component
       function obj = env_mix_energy(obj , T)  
           % T is the duration of the load cycle in seconds
           
           % Iterate through each load cycle
           for i = 1 : numel(obj.mdot(:,1))
               
               % Only evaluate if mdot is defined
               if ~isempty(obj.mdot(i,1))
                   obj.W(i,1)    = obj.w(i,1)    * obj.mdot(i,1) * T(i) ;
                   obj.Q(i,1)    = obj.q(i,1)    * obj.mdot(i,1) * T(i) ;
                   obj.DH(i,1)   = obj.Dh(i,1)   * obj.mdot(i,1) * T(i) ;
                   obj.Sirr(i,1) = obj.sirr(i,1) * obj.mdot(i,1) * T(i) ;
               end
               
           end
       end
       
        
    end
    
end