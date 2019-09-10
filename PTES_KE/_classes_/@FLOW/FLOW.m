classdef FLOW 
    % FLOW is a construct that contains the fluid flow properties at a
    % given point in the cycle. These properties are calculated using
    % information from the FLUID class.
    
    properties
        name    % name
        note    % Any other details
        T       % Temperature
        P       % pressure
        h       % enthalpy
        s       % entropy
        b       % exergy
        q       % Quality
        rho     % density
        mdot    % mass flow rate
        n       % Number of observations - sizes arrays
    end
    
    % Constructor function
    methods
        function obj = FLOW(flow_dat) 
                   
            if isfield(flow_dat , 'nobv')            
                obj.n = flow_dat.nobv ;               
            else                
                obj.n = 1 ;                
                %warning('FLOW class is not provided with the number of observations');               
            end           
            
            % Set remaining values to zero            
            obj.T = zeros(obj.n,1) ;            
            obj.P = zeros(obj.n,1) ;            
            obj.s = zeros(obj.n,1) ;            
            obj.h = zeros(obj.n,1) ;            
            obj.b = zeros(obj.n,1) ;
            obj.q = zeros(obj.n,1) ;            
            obj.rho = zeros(obj.n,1) ;           
            obj.mdot = zeros(obj.n,1) ;
               
        end
        
    end
    
    % Other functions
    methods
        % This function calculates the remaining fluid properties
        % p is a character array that indicates which flow properties are
        % known and therefore what needs to be calculated
        % fld is the fluid class that is passed to this function
        function obj = flow_state(obj,p,fld)
            % Change things into the correct case            
            if numel(p) > 2
                error('flow_state function called with too many property inputs, p = %f > 2',numel(p));
            end
            
            % I expect the following is extremely buggy :/
            % Have to calculate properties for each of observation
            for i = 1 : obj.n
            
            if true(fld.CoolProps)
                
              if contains(['pT','Tp'],p,'IgnoreCase',true)
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.h(i) = CoolProp(fld,'H',in) ;  
                  obj.s(i) = CoolProp(fld,'S',in) ;        
                  obj.q(i) = CoolProp(fld,'Q',in) ;  
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              elseif contains(['ph','hp'],p,'IgnoreCase',true)
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.T(i) = CoolProp(fld,'T',in) ;        
                  obj.s(i) = CoolProp(fld,'S',in) ;        
                  obj.q(i) = CoolProp(fld,'Q',in) ; 
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              elseif contains(['ps','sp'],p,'IgnoreCase',true)
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.h(i) = CoolProp(fld,'H',in) ;       
                  obj.T(i) = CoolProp(fld,'T',in) ;        
                  obj.q(i) = CoolProp(fld,'Q',in) ; 
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              elseif contains(['pQ','Qp'],p,'IgnoreCase',true) 
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.h(i) = CoolProp(fld,'H',in) ;        
                  obj.s(i) = CoolProp(fld,'S',in) ;        
                  obj.T(i) = CoolProp(fld,'T',in) ; 
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              elseif contains(['Th','hT'],p,'IgnoreCase',true)
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.P(i) = CoolProp(fld,'P',in) ;        
                  obj.s(i) = CoolProp(fld,'S',in) ;        
                  obj.q(i) = CoolProp(fld,'Q',in) ; 
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              elseif contains(['Ts','sT'],p,'IgnoreCase',true)
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.h(i) = CoolProp(fld,'H',in) ;        
                  obj.P(i) = CoolProp(fld,'P',in) ;        
                  obj.q(i) = CoolProp(fld,'Q',in) ; 
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              elseif contains(['TQ','QT'],p,'IgnoreCase',true)
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.h(i) = CoolProp(fld,'H',in) ;        
                  obj.s(i) = CoolProp(fld,'S',in) ;        
                  obj.P(i) = CoolProp(fld,'P',in) ; 
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              elseif contains(['hs','sh'],p,'IgnoreCase',true)
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.T(i) = CoolProp(fld,'T',in) ;        
                  obj.P(i) = CoolProp(fld,'p',in) ;        
                  obj.q(i) = CoolProp(fld,'Q',in) ; 
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              elseif contains(['hQ','Qh'],p,'IgnoreCase',true)
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.T(i) = CoolProp(fld,'T',in) ;        
                  obj.s(i) = CoolProp(fld,'S',in) ;        
                  obj.P(i) = CoolProp(fld,'P',in) ; 
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              elseif contains(['sQ','Qs'],p,'IgnoreCase',true)
                  in.v1 = eval(strcat('obj.',p(1),'(',int2str(i),')')) ;     % Is this the dodgiest code ever?
                  in.v2 = eval(strcat('obj.',p(2),'(',int2str(i),')')) ;
                  in.n1 = upper(p(1));
                  in.n2 = upper(p(2));
                  obj.h(i) = CoolProp(fld,'H',in) ;        
                  obj.T(i) = CoolProp(fld,'T',in) ;        
                  obj.P(i) = CoolProp(fld,'P',in) ; 
                  obj.rho(i) = CoolProp(fld,'D',in) ;
              else
                 error(['Function flow state called with unknown properties.' ... 
                        'p should be a string comprising two of the following.' ...
                        'T p h s Q. Note, this is NOT case sensitive.']);
              end
              
            else
              error('Currently fluid properties can only be calculated with CoolProp. See FLUID class.');
            end
            
            % Calculate exergy
            obj.b(i) = obj.h(i) - fld.Tref * obj.s(i) ; 
            end
        end
    end
    
end

