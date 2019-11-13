classdef stream
    %STREAN CLASS Contains all the properties and geometrical details of a
    %fluid stream in the HEX
    
    properties
        % General variables
        name
        handle
        Tin
        Tout
        pin
        pout
        mdot
        Cp_mean        
        % Lookup table matrices
        TX
        pX
        CpX
        hX
        sX
        rhoX
        muX
        kX
        PrX
        % Arrays in HEX sections (sides)
        T
        p
        Cp
        h
        s
        rho
        v
        mu
        k
        Pr
        Re
        Cf
        St
        ht %heat transfer coeff
        % Averaged arrays (in sections centers)
        T_AV
        rho_AV
        Cp_AV
        % Geometry
        A
        Af
        Ax
        D
        G
    end
    
    methods
    end
    
end

