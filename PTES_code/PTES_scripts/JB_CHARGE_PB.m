% Set stage indices
iG  = 1;  % keeps track of the gas stage number
iE1 = 1;  % keeps track of the heat rejection stream number 1
iE2 = 1;  % keeps track of the heat rejection stream number 2
iA  = 1;  % keeps track of the Air (heat rejection) stream number

% Expander inlet indices
iExp = 1 + Nc_ch*2;         % index expander inlet (after compression + cooling)

str1 = '        PACKED BED: TIMESTEP # %8i\n' ;

if design_mode == 1
    % Initial guess of charge conditions
    % Compressor inlet (regenerator hot outlet)
    gas.state(iL,1).p    = pbot; gas.state(iL,1).T = T1;
    gas.state(iL,1).mdot = Load.mdot(iL);
    [gas] = update(gas,[iL,1],1);
    
    gas.state(iL,iExp).T = T0;
    gas.state(iL,iExp).p = pbot;
    gas.state(iL,iExp).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iExp],1);
else
    % Set up matrix that guesses the converged solution
    for ii = 1 : numel(C(C~=0))/2
        gas.state(iL,ii).T    = C(1,ii) ;
        gas.state(iL,ii).p    = C(2,ii) ;
        gas.state(iL,ii).mdot = Load.mdot(iL) ;
        
        % For inventory control, assume that the pressure scales with the off-design mass flow rate
        gas.state(iL,ii).p = CCMP.Pin * Load.mdot(iL) / CCMP.mdot0 ;
        [gas] = update(gas,[iL,ii],1);
    end  
end

% Set initial values of packed beds 
pbH0 = pbH ;
pbC0 = pbC ;

% Set matrix of temperature and pressure points to test convergence
C_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
counter = 1;
while 1
    
    fprintf(1,['Charging JB PTES (PACKED BED). Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
        
    Lend = false ; % This is the criterion that decides when the packed beds are charged
    istep  = 1 ; % This counts timesteps
    iplotH = 1 ; % This counts out when hot thermal profiles are plotted
    iplotC = 1 ; % This counts out when cold thermal profiles are plotted
    time = 0 ;    
    
    % Set packed beds to original values at start of iteration
    pbH = pbH0 ;
    pbC = pbC0 ;
    
    % Calculate the energy in each packed bed at the start
    pbH = PB_ENERGY(pbH,gas,iL,'start') ;
    pbC = PB_ENERGY(pbC,gas,iL,'start') ;
    
    % Initialize stage counters that sum energy over timesteps (sets to zero)
    CCMPtot  = stage_class ;
    CEXPtot  = stage_class ;
    stagetot(1:numel(gas.stage(iL,:))) = stage_class ;
    
    while ~Lend
        
        if mod(istep,50) == 0
            fprintf(1,str1,istep) ;
        end
    
        for iN = 1:Nc_ch
            % COMPRESS
            if setTmax
                T_aim = Tmax;
                [CCMP(iN),gas,iG] = compexp_func (CCMP(iN),iL,gas,iG,'Taim',T_aim, design_mode) ; %#ok<*SAGROW>
                PRch = gas.state(iL,iG).p / gas.state(iL,iG-1).p ;
            else
                if design_mode == 1
                    PRc = (pmax/gas.state(iL,iG).p)^(1/(Nc_ch+1-iN)); % stage compression pressure ratio
                else
                    PRc = CCMP.pr0 ;
                end
                p_aim = gas.state(iL,iG).p*PRc;
                [CCMP(iN),gas,iG] = compexp_func (CCMP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ;
            end
            ptop  = gas.state(iL,iG).p;
            
            % Put hot gas through hot packed bed
            %pbH = PB_TIMESTEP(pbH, gas, 'chg') ;
            pbH = PB_TIMESTEP_IDEAL(pbH, gas, iL, iG, 'chg') ;
            
            pbH = PB_FLUX(pbH, T0, p0, gas, iL) ;
        
            % Fluid outlet
            iG = iG + 1 ;
            gas.state(iL,iG).T = pbH.TF(pbH.NX,1) ;
            gas.state(iL,iG).p = pbH.P(pbH.NX,1) ;
        
            % Update fluid
            [gas] = update(gas,[iL,iG],1);
        
            % Now put through a heat exchanger - reject heat
            T_aim = T0;
            [gas,environ,iG,iE1] = hex_set(gas,[iL,iG],environ,[iL,iE1],T_aim,1.0,0.01);
        
        end
        
        for iN = 1:Ne_ch
            % EXPAND
            PRe = (gas.state(iL,iG).p/pbot)^(1/(Ne_ch+1-iN)); % stage expansion pressure ratio
            p_aim = gas.state(iL,iG).p/PRe;
            [CEXP(iN),gas,iG] = compexp_func (CEXP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ;
            
            % Pass cold gas through cold packed bed
            %pbC = PB_TIMESTEP(pbC, gas, 'chg') ;
            pbC = PB_TIMESTEP_IDEAL(pbC, gas, iL, iG, 'chg') ;
            
            pbC = PB_FLUX(pbC, T0, p0, gas, iL) ;
        
            % Fluid outlet
            iG = iG + 1 ;
            gas.state(iL,iG).T = pbC.TF(pbC.NX,1) ;
            gas.state(iL,iG).p = pbC.P(pbC.NX,1) ;
        
            % Update fluid
            [gas] = update(gas,[iL,iG],1);
        
            % Now put through a heat exchanger - add heat
            T_aim = T0;
            [gas,environ,iG,iE2] = hex_set(gas,[iL,iG],environ,[iL,iE2],T_aim,1.0,0.01);
        end
        
        Lend = PB_STOP_PHASE(pbH,pbC,time,Load.type(iL)) ; % Stop charge or discharge?
                
        % Plot out results if appropriate
        if (istep == pbH.Tprof(iplotH)) || Lend
            pbH.TSprof(:,iplotH,iL) = pbH.TS(:,1) ;
            pbH.TFprof(:,iplotH,iL) = pbH.TF(:,1) ;
            iplotH        = iplotH + 1 ;
        end
        
        if (istep == pbC.Tprof(iplotC)) || Lend
            pbC.TSprof(:,iplotC,iL) = pbC.TS(:,1) ;
            pbC.TFprof(:,iplotC,iL) = pbC.TF(:,1) ;
            iplotC        = iplotC + 1 ;
        end
                
        % Now add the incremental energy terms to the total terms
        [~,CCMPtot]  = SUM_ENERGY(CCMP,CCMPtot,1,iL) ;
        [~,CEXPtot]  = SUM_ENERGY(CEXP,CEXPtot,1,iL) ;
        [~,stagetot] = SUM_ENERGY(gas.stage(iL,:),stagetot,2,numel(gas.stage(iL,:))) ;
                
        istep = istep + 1 ;
        time  = time + pbH.dt ;
        iG=1; iE1=1; iE2=1;
        
    end
    
    % Calculate the energy in each packed bed at the end
    pbH = PB_ENERGY(pbH,gas,iL,'end') ;
    pbC = PB_ENERGY(pbC,gas,iL,'end') ;
    
    % Assess losses
    pbH = PB_LOSS(pbH,iL) ;
    pbC = PB_LOSS(pbC,iL) ;
    
    % Assign total losses to actual components
    [CCMP,~]  = SUM_ENERGY(CCMP,CCMPtot,3,iL) ;
    [CEXP,~]  = SUM_ENERGY(CEXP,CEXPtot,3,iL) ;
    [gas.stage(iL,:),~] = SUM_ENERGY(gas.stage(iL,:),stagetot,4,numel(gas.stage(iL,:))) ;
        
    % Determine convergence and proceed
    C = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    
    if (all(abs((C(C~=0) - C_0(C~=0))./C(C~=0))*100 < 1e-3) || counter > 100) % is charge cycle converged?

        % Close working fluid streams
        gas.stage(iL,iG).type = 'end';
        gas = count_Nstg(gas);
        
        % Close air (heat rejection) streams
        iA_out = 1:2:(iA-1); iA_in  = iA_out + 1;
        for i=iA_in, air.stage(iL,i).type = 'end'; end
        air = count_Nstg(air);
        
       
        % Uncomment these lines to print states
        print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        
        % All heat exchangers now have their geometry set up properly
        for ihx = 1 : length(HX)
            HX(ihx).Lgeom_set = true ;
        end
        
        [pbH,pbC] = PB_REVERSE(pbH,pbC,Load,iL) ;
        
        % Exit loop
        break
    else
        
        if ~design_mode
            
            %print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
            
            % Adjust inlet pressure to try to reach convergence. The
            % 'smoothing' factor has to be quite small (<0.1, say) for this to be stable
%             fprintf('p1:   %13.8f \n',gas.state(iL,1).p/1e5)
%             fprintf('pEND: %13.8f \n\n',gas.state(iL,iG).p/1e5)
%             fprintf('T1:   %13.8f \n',gas.state(iL,1).T)
%             fprintf('TEND: %13.8f \n\n',gas.state(iL,iG).T)
            gas.state(iL,1).p = gas.state(iL,1).p - 0.10 * (gas.state(iL,iG).p - gas.state(iL,1).p) ;
            gas.state(iL,1).T = gas.state(iL,1).T + 0.10 * (gas.state(iL,iG).T - gas.state(iL,1).T) ;
                      
            gas.state(iL,1).mdot = Load.mdot(iL);
            [gas] = update(gas,[iL,1],1);
            
        else
            gas.state(iL,1) = gas.state(iL,iG);
        end
        
        C_0 = C;
        counter = counter + 1 ;
        
    end
end

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Atmospheric tanks
iA_out = 0; iA_in = 0;
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);