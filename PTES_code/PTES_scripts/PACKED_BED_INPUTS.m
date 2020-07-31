%%% SET PACKED BED INPUT VARIABLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lcyclic = true ; % Run the system until cyclic convergence reached
Ncyc    = 20 ; % Maximum number of cycles to run before giving up
Icyc    = 1 ; % Number of cycles that have elapsed so far

% Set up the classes
for ii = 1 : Nhot
    pbH(ii) = packbed_class('hot') ;
end

for ii = 1 : Ncld
    pbC(ii) = packbed_class('cold') ;
end

Ltime = false ; % End cycle after a certain time has elapsed
Ltext = true ; % End cycle after the exit temperature reaches a certain threshold

% *** HOT PACKED BEDS *** %
for ii = 1 : Nhot
   
    % Constant properties with temperature?
    pbH(ii).Lconst = false ;
    
    % Temperatures
    pbH(ii).TC = TH_chg0 ;     % Charged temperature, K
    pbH(ii).TD = TH_dis0 ;      % Discharged temperature, K
    
    % GEOMETRY
    pbH(ii).dp  = 0.02 ;    % Particle diameter
    pbH(ii).eps = 0.4;      % Void fraction
    pbH(ii).AR  = 1.0;      % Tank aspect ratio
    
    % Solid properties
    % Replace these later with data from an input file
    pbH(ii).Sname = 'Magnetite' ;
    
    % Nominal charging time
    pbH(ii).tN = 8.0 * 3600.0 ; % Nominal charging time (to fully charge storage)
    pbH(ii).Ncyc = 1 ;
    
    % Cycle control
    pbH(ii).Ltime = Ltime ; % End cycle after a certain time has elapsed
    pbH(ii).Ltext = Ltext ; % End cycle after the exit temperature reaches a certain threshold
    
    pbH(ii).timeC = 6.0 * 3600. ; % End charge cycle after this time, s
    pbH(ii).timeD = 6.0 * 3600. ; % End discharge cycle after this time, s
    
    pbH(ii).textC = 0.25 ; % End charge cycle after this temperature threshold is exceeded
    pbH(ii).textD = 0.25 ; % End discharge cycle after this temperature threshold is exceeded
    
    % Grid and timestep controls
    pbH(ii).DELX = (1/100) ; % Size of grid-step relative to length of storage, dx / L
    pbH(ii).CFL  = 200; % Courant-Friedrichs-Levy number, used to set time step size. (Choose 5 for oil with ideal gas routine. Choose 100 for ideal gas with ideal gas routine or 10 if using liquid routine. Choose 0.1 for liquids with liquid routine.
    pbH(ii).TMAX = 4.0 ;  % Max time to run calculations for - this is a multiple of tN
    
    % Turn conduction terms on/off 
    pbH(ii).Cond = 0 ; % Conduction is on if 1, off if 0
    
    % Plotting controls
    pbH(ii).Nprof = 5 ; % Number of temperature profiles to save
    
    pbH(ii).Lideal = true ; % Model gas in packed bed as an ideal gas?
        
end


% *** COLD PACKED BEDS *** %
for ii = 1 : Ncld
   
    % Constant properties with temperature?
    pbC(ii).Lconst = false ;
    
    % Temperatures
    pbC(ii).TC = TC_chg0 ;     % Charged temperature, K
    pbC(ii).TD = TC_dis0 ;      % Discharged temperature, K
    
    % GEOMETRY
    pbC(ii).dp  = 0.02 ;    % Particle diameter
    pbC(ii).eps = 0.4;      % Void fraction
    pbC(ii).AR  = 1.0;      % Tank aspect ratio
    
    % Solid properties
    % Replace these later with data from an input file
    pbC(ii).Sname = 'Magnetite' ;
    
    % Nominal charging time
    pbC(ii).tN = 8.0 * 3600.0 ; % Nominal charging time (to fully charge storage)
    pbC(ii).Ncyc = 1 ;
    
    % Cycle control
    pbC(ii).Ltime = Ltime ; % End cycle after a certain time has elapsed
    pbC(ii).Ltext = Ltext ; % End cycle after the exit temperature reaches a certain threshold
    
    pbC(ii).timeC = 6.0 * 3600. ; % End charge cycle after this time, s
    pbC(ii).timeD = 6.0 * 3600. ; % End discharge cycle after this time, s
    
    pbC(ii).textC = 0.25 ; % End charge cycle after this temperature threshold is exceeded
    pbC(ii).textD = 0.25 ; % End discharge cycle after this temperature threshold is exceeded
    
    % Grid and timestep controls
    pbC(ii).DELX = (1/100) ; % Size of grid-step relative to length of storage, dx / L
    pbC(ii).CFL  = 500; % Courant-Friedrichs-Levy number, used to set time step size. (Choose 5 for oil with ideal gas routine. Choose 100 for ideal gas with ideal gas routine or 10 if using liquid routine. Choose 0.1 for liquids with liquid routine.
    pbC(ii).TMAX = 4.0 ;  % Max time to run calculations for - this is a multiple of tN
    
    % Turn conduction terms on/off 
    pbC(ii).Cond = 0 ; % Conduction is on if 1, off if 0
    
    % Plotting controls
    pbC(ii).Nprof = 5 ; % Number of temperature profiles to save

    pbC(ii).Lideal = true ; % Model gas in packed bed as an ideal gas?
        
end
