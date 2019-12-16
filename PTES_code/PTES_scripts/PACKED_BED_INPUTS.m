%%% SET PACKED BED INPUT VARIABLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nhot = 1;
% Set up the classes
for ii = 1 : Nhot
    pbH(ii) = packbed_class('hot') ;
end

%for ii = 1 : Ncld
%    pbC(ii) = packbed_class('cold') ;
%end



% *** HOT PACKED BEDS *** %
for ii = 1 : Nhot
   
    % GEOMETRY
    pbH(ii).dp  = 0.02 ;    % Particle diameter
    pbH(ii).eps = 0.4;      % Void fraction
    pbH(ii).AR  = 1.0;      % Tank aspect ratio
    
    % Solid properties
    % Replace these later with data from an input file
    pbH(ii).kS   = 41.7 ;    % Thermal conductivity, W/mK
    pbH(ii).rhoS = 5000 ;   % Density, kg/m3
    pbH(ii).cS   = 800 ;    % Specific heat capacity, J/kgK
    
    % Temperatures
    pbH(ii).TC = 550 + 273.15 ;     % Charged temperature, K
    pbH(ii).TD = 30 + 273.15 ;      % Discharged temperature, K
    
    % Fluid properties
    % Eventually this will be an input from external routines
    pbH(ii).kF   = 0.035 ; % Thermal conductivity, W/mK
    pbH(ii).rhoF = 1.003 ; % Density, kg/m3
    pbH(ii).cF   = 1000 ; % Specific heat capcaity, J/kgK
    pbH(ii).Pr   = 0.7 ;  % Prandtl number
    pbH(ii).mu   = 1.4e-5 ; % Viscosity, Pa.s
    pbH(ii).mdot = 10.0 ; % Fluid mass flow rate, kg/s
    pbH(ii).Pin  = 1e5 ;
    
    % Nominal charging time
    pbH(ii).tN = 8.0 * 3600.0 ; % Nominal charging time (to fully charge storage)
    pbH(ii).Ncyc = 50 ;
    
    % Cycle control
    pbH(ii).Ltime = false ; % End cycle after a certain time has elapsed
    pbH(ii).Ltext = true ; % End cycle after the exit temperature reaches a certain threshold
    
    pbH(ii).timeC = 8.0 * 3600. ; % End charge cycle after this time, s
    pbH(ii).timeD = 8.0 * 3600. ; % End discharge cycle after this time, s
    
    pbH(ii).textC = 0.25 ; % End charge cycle after this temperature threshold is exceeded
    pbH(ii).textD = 0.25 ; % End discharge cycle after this temperature threshold is exceeded
    
    % Grid and timestep controls
    pbH(ii).DELX = 0.01 ; % Size of grid-step relative to length of storage, dx / L
    pbH(ii).CFL  = 1000  ; % Courant-Friedrichs-Levy number, used to set time step size
    pbH(ii).TMAX = 2.0 ;  % Max time to run calculations for - this is a multiple of tN
    
    % Turn conduction terms on/off (can speed up and stabalize routine
    pbH(ii).Cond = 0 ; % Conduction is on if 1, off if 0
    
    % Plotting controls
    pbH(ii).Nprof = 5 ; % Number of temperature profiles to save
    
    % Set up packed bed
    pbH(ii) = PB_INITIALISE( pbH(ii) ) ;
    
end


% *** COLD PACKED BEDS *** %

