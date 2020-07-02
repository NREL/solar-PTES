% Class for miscellaneous components that don't require lots of info
classdef misc_class
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
        function obj = misc_class(type, numPeriods)
            obj.type = type ;
            
            obj.mdot = zeros(numPeriods,2) ;
            
            obj.w    = zeros(numPeriods,2) ;
            obj.q    = zeros(numPeriods,2) ;
            obj.Dh   = zeros(numPeriods,2) ;
            obj.sirr = zeros(numPeriods,2) ;
                        
            obj.W    = zeros(numPeriods,2) ;
            obj.Q    = zeros(numPeriods,2) ;
            obj.DH   = zeros(numPeriods,2) ;
            obj.Sirr = zeros(numPeriods,2) ;
            
        end
        
        % Calculate the energy and entropy terms for the miscellaneous  component
       function obj = misc_energy(obj , T)  
           % T is the duration of the load cycle in seconds
           
           % Iterate through each load cycle
           for i = 1 : numel(obj.mdot(:,1))
               
               % Only evaluate if mdot is defined
               if ~isempty(obj.mdot(i,1))
                   % Hot streams
                   obj.W(i,1)    = obj.w(i,1)    * obj.mdot(i,1) * T(i) ;
                   obj.Q(i,1)    = obj.q(i,1)    * obj.mdot(i,1) * T(i) ;
                   obj.DH(i,1)   = obj.Dh(i,1)   * obj.mdot(i,1) * T(i) ;
                   obj.Sirr(i,1) = obj.sirr(i,1) * obj.mdot(i,1) * T(i) ;
                   
                   % Cold streams
                   obj.W(i,2)    = obj.w(i,2)    * obj.mdot(i,2) * T(i) ;
                   obj.Q(i,2)    = obj.q(i,2)    * obj.mdot(i,2) * T(i) ;
                   obj.DH(i,2)   = obj.Dh(i,2)   * obj.mdot(i,2) * T(i) ;
                   obj.Sirr(i,2) = obj.sirr(i,2) * obj.mdot(i,2) * T(i) ;
               end
               
           end
       end
       
        
    end
    
end