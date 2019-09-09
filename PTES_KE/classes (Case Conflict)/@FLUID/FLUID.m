classdef FLUID 
    % FLUID Class containing data on a specific fluid from which fluid
    % properties can be calculated
    
    properties
        name
        note
        Tref    % The reference temperature
        href    % The reference enthalpy = 0 at Tref,pref
        pref    % The reference pressure
        sref    % The reference entropy = 0 at Tref,pref
        Cost
        CoolProps
    end
    
    methods
        %   Constructor function
        function obj = FLUID(f_data)
            if nargin == 1 
                
                if isfield(f_data,'name')
                    obj.name = f_data.name ; % Assign the name of the fluid
                end
                
                if isfield(f_data,'CoolProps')
                    obj.CoolProps = f_data.CoolProps ;
                else
                    obj.CoolProps = true ;      % Default value
                end
                
                % Assign reference values
                if isfield(f_data,'Tref')
                    obj.Tref = f_data.Tref ;
                else
                    obj.Tref = 273.15 + 25.0 ;
                end
                
                if isfield(f_data,'pref')
                    obj.pref = f_data.pref ;
                else
                    obj.pref = 1.0e5 ;
                end
                
                if isfield(f_data,'href')
                    obj.href = f_data.href ;
                else
                    obj.href = 0.0 ;
                end
                
                if isfield(f_data,'sref')
                    obj.sref = f_data.sref ;
                else
                    obj.sref = 0.0 ;
                end
                
                %   Determine whether fluid properties will be calculated
                %   with CoolProps or not
                if true(obj.CoolProps)
                    obj.note = 'Fluid properties are calculated with Coolprops' ;                    
                else
                    obj.note = 'Fluid properties are calculated somehow ...' ;
                    error('Currently fluid properties can only be calculated with CoolProp');
                end
                                   
                if isfield(f_data,'Cost')
                    obj.Cost = f_data.Cost ;
                else
                    obj.Cost = 0.0 ;
                end
            else
               error('Error in FLUID class constructor. Fluid data not provided') 
            end
        end
        
        % Function that calculates fluid properties using CoolProps
        function f = CoolProp(obj,out,in)
           % Evaluate coolprop function
           f = CoolProp.PropsSI(out,in.n1,in.v1,in.n2,in.v2,obj.name) ; 
        end
        
        % Function that decides which function to call to calculate fluid
        % properties
        function x = fprops(obj,out,in)
           if true(obj.CoolProps)  
               x = CoolProp(obj,out,in) ;
           else
               error('Currently fluid properties can only be calculated with CoolProp');
           end
        end
    end
    
end

