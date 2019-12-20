classdef stream
    %STREAM CLASS Contains all the properties and geometrical details of a
    %fluid stream in the HEX
    
    properties
        % General variables
        name
        handle
        Tin     {mustBePositive}
        Tout    {mustBePositive}
        pin     {mustBePositive}
        pout    {mustBePositive}
        mdot    {mustBePositive}
        Cp_mean {mustBePositive}
        % Arrays in HEX sections (sides)
        T       {mustBePositive}
        p       {mustBePositive}
        Cp      {mustBePositive}
        h       {mustBePositive}
        s       {mustBePositive}
        rho     {mustBePositive}
        v       {mustBePositive}
        mu      {mustBePositive}
        k       {mustBePositive}
        Pr      {mustBePositive}
        Re      {mustBePositive}
        Cf      {mustBePositive}
        St      {mustBePositive}
        ht %heat transfer coeff
        % Averaged arrays (in sections centers)
        T_AV    {mustBePositive}
        rho_AV  {mustBePositive}
        Cp_AV   {mustBePositive}
        % Geometry
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

