classdef stream
    %STREAM CLASS Contains all the properties and geometrical details of a
    %fluid stream in the HEX
    
    properties
        % General variables
        name
        read
        handle
        HEOS
        Tin     {mustBePositive}
        Tout    {mustBePositive}
        pin     {mustBePositive}
        pout    {mustBePositive}
        mdot
        Cp_mean {mustBePositive}
        % Arrays in HEX sections (sides)
        T       {mustBePositive}
        p       {mustBePositive}
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
        % Arrays for the two-phase region
        hLG  % latent heat of vaporisation
        rhoL % density of saturated liquid
        rhoG % density of saturated vapour
        kL   % conductivity of saturated liquid
        muL  % viscosity of saturated liquid
        PrL  % Prandtl of saturated liquid
        q    % heat flux
        % Averaged arrays (in sections centers)
        T_AV    {mustBePositive}
        rho_AV  {mustBePositive}
        Cp_AV   {mustBePositive}
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

