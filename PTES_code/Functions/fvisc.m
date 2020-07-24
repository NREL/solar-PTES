        % Calculate the viscosity of an ideal gas using Sutherland's law
        function mu = fvisc(fluid,T)
            
            mu0   = fluid.IDL.mu0 ;
            TVref = fluid.IDL.TVref ;
            S     = fluid.IDL.S ;
            
            mu    = mu0 .* ((TVref + S) ./ (T + S)) .* (T./TVref).^1.5 ;
            
        end
