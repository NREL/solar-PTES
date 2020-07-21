classdef stream
    %STREAM CLASS Contains all the properties and geometrical details of a
    %fluid stream in the HEX
    
    properties
        % General variables
        name
        read
        handle
        HEOS
        IDL
        Tin     {mustBePositive}
        Tout    {mustBePositive}
        pin     {mustBePositive}
        pout    {mustBePositive}
        hin
        hout
        mdot
        Cp_mean {mustBePositive}
        % Arrays in HEX section corners
        T
        p
        x  % vapour quality
        Cp      
        h       
        s       
        rho     
        v       
        mu      
        k       
        Pr      %{mustBePositive}
        Re      
        Cf      
        St      
        ht %heat transfer coeff
        
        %%% Arrays for the two-phase region
        % Saturated liquid
        rhoL
        kL
        muL
        PrL
        % Saturated vapour
        rhoG
        kG
        muG
        PrG
        % Latent heat of vaporisation
        hLG
        % Heat flux
        q
        %%%
        
        % Averaged arrays (in section centers)
        T_AV    {mustBePositive}
        rho_AV  {mustBePositive}
        Cp_AV   {mustBePositive}
        dpdL
        % Geometry
        shape
        A       {mustBePositive}
        Af      {mustBePositive}
        Ax      {mustBePositive}
        D       {mustBePositive}
        G       {mustBePositive}
        L       {mustBePositive}
    end
    
    methods
    end
    
end

