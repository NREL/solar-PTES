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
   
    % Constant properties with temperature?
    pbH(ii).Lconst = false ;
    
    % Temperatures
    pbH(ii).TC = 250 + 273.15 ;     % Charged temperature, K
    pbH(ii).TD = 30 + 273.15 ;      % Discharged temperature, K
    
    % GEOMETRY
    pbH(ii).dp  = 0.02 ;    % Particle diameter
    pbH(ii).eps = 0.4;      % Void fraction
    pbH(ii).AR  = 1.0;      % Tank aspect ratio
    
    % Solid properties
    % Replace these later with data from an input file
    pbH(ii).Sname = 'Magnetite' ;
    pbH(ii).sld   = packbed_class.create_solid_table(pbH(ii).Sname) ;
    pbH(ii).kS    = pbH(ii).sld(1,6) ;    % Thermal conductivity, W/mK
    pbH(ii).rhoS  = 1./pbH(ii).sld(1,3) ;   % Density, kg/m3
    
    x  = pbH(ii).sld(:,1); %temperatures
    h1 = interp1(x,pbH(ii).sld(:,2),pbH(ii).TC);
    h2 = interp1(x,pbH(ii).sld(:,2),pbH(ii).TD);
        
    pbH(ii).cS   = (h1 - h2) / (pbH(ii).TC - pbH(ii).TD) ;    % Specific heat capacity, J/kgK
    
    % Fluid properties
    % Eventually this will be an input from external routines
    pbH(ii).kF   = 0.035 ; % Thermal conductivity, W/mK
    pbH(ii).rhoF = 857  ; % Density, kg/m3 - Need to calculate these properly!
    pbH(ii).cF   = 2300 ; % Specific heat capacity, J/kgK
    pbH(ii).Pr   = 0.7 ;  % Prandtl number
    pbH(ii).mu   = 5e-4 ; % Viscosity, Pa.s
    pbH(ii).mdot = 10.0 ; % Fluid mass flow rate, kg/s
    pbH(ii).Pin  = 1e5 ;
    
    % Nominal charging time
    pbH(ii).tN = 8.0 * 3600.0 ; % Nominal charging time (to fully charge storage)
    pbH(ii).Ncyc = 1 ;
    
    % Cycle control
    pbH(ii).Ltime = false ; % End cycle after a certain time has elapsed
    pbH(ii).Ltext = true ; % End cycle after the exit temperature reaches a certain threshold
    
    pbH(ii).timeC = 8.0 * 3600. ; % End charge cycle after this time, s
    pbH(ii).timeD = 8.0 * 3600. ; % End discharge cycle after this time, s
    
    pbH(ii).textC = 0.95 ; % End charge cycle after this temperature threshold is exceeded
    pbH(ii).textD = 0.95 ; % End discharge cycle after this temperature threshold is exceeded
    
    % Grid and timestep controls
    pbH(ii).DELX = (1/100) ; % Size of grid-step relative to length of storage, dx / L
    pbH(ii).CFL  = 0.1  ; % Courant-Friedrichs-Levy number, used to set time step size
    pbH(ii).TMAX = 4.0 ;  % Max time to run calculations for - this is a multiple of tN
    
    % Turn conduction terms on/off 
    pbH(ii).Cond = 1 ; % Conduction is on if 1, off if 0
    
    % Plotting controls
    pbH(ii).Nprof = 5 ; % Number of temperature profiles to save
    
    % Fluid
    FSname  = 'MineralOil';  % fluid name
    fluidS = fluid_class(FSname,'SF','TAB',NaN,1,30); % Storage fluid
    
    % Set up packed bed
    pbH(ii) = PB_INITIALISE( pbH(ii), fluidS ) ;
    
end


% *** COLD PACKED BEDS *** %

