% Set stage indices
iG  = 1;  % keeps track of the gas stage number
iE1 = 1;  % keeps track of the Environment (heat rejection) stream number
iE2 = 1;  % keeps track of the Environment (heat rejection) stream number
iA  = 1;  % keeps track of the Air (heat rejection) stream number

str1 = '        PACKED BED: TIMESTEP # %8i\n' ;

% Expander inlet indices
iCmp = 2 + Nc_dis*2;         % index compressor outlet inlet (after cooling, compression, and heat rejection)

% Compute PR_dis based on charge pressure ratio and PRr
PRdis = PRr*PRch;

if design_mode == 1
    % Initial guess of discharge conditions
    % Expander outlet (regenerator hot inlet)
    gas.state(iL,1).p    = pbot; gas.state(iL,1).T = T1;
    gas.state(iL,1).mdot = Load.mdot(iL);
    [gas] = update(gas,[iL,1],1);
    gas.state(iL,iCmp).T = T0;
    gas.state(iL,iCmp).p = gas.state(iL,1).p*PRdis;
    gas.state(iL,iCmp).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iCmp],1);
else
    for ii = 1 : numel(D(D~=0))/2
        gas.state(iL,ii).T    = D(1,ii) ;
        gas.state(iL,ii).p    = D(2,ii) ;
        gas.state(iL,ii).mdot = Load.mdot(iL) ;
        
        % For inventory control, assume that the pressure scales with the off-design mass flow rate
        gas.state(iL,ii).p = (DEXP.Pin/DEXP.pr0) * Load.mdot(iL) / DEXP.mdot0 ;
        
        [gas] = update(gas,[iL,ii],1);
    end    
end

% Set initial values of packed beds 
pbH0 = pbH ;
pbC0 = pbC ;

% Set matrix of temperature and pressure points to test convergence
D_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
counter = 1 ;
while 1
    fprintf(1,['Discharging JB PTES (PACKED BED). Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    
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
    DCMPtot  = stage_class ;
    DEXPtot  = stage_class ;
    stagetot(1:numel(gas.stage(iL,:))) = stage_class ;
    
    while ~Lend
        
        if mod(istep,50) == 0
            fprintf(1,str1,istep) ;
        end

       
        for iN = 1:Nc_dis
            % REJECT HEAT (external HEX)
            % REPLACE THIS WITH hex_func call?
            T_aim = environ.T0;
            [gas,environ,iG,iE1] = hex_set(gas,[iL,iG],environ,[iL,iE1],T_aim,1.0,ploss);
            
            switch Load.mode
                case 0 % PTES
                    % COOL in the cold packed bed
                    %pbC = PB_TIMESTEP(pbC, gas, 'dis') ;
                    pbC = PB_TIMESTEP_IDEAL(pbC, gas, iL, iG, 'dis') ;
                    
                    pbC = PB_FLUX(pbC, T0, p0, gas, iL) ;
                    
                    % Fluid outlet
                    iG = iG + 1 ;
                    gas.state(iL,iG).T = pbC.TF(pbC.NX,1) ;
                    gas.state(iL,iG).p = pbC.P(pbC.NX,1) ;
                    
                    % Update fluid
                    [gas] = update(gas,[iL,iG],1);
                    
                case 1 % Heat engine only
            end
            
            % COMPRESS
            PRc_dis = (PRdis)^(1/Nc_dis)/(1-ploss)^2; % stage compression pressure ratio
            p_aim = gas.state(iL,iG).p*PRc_dis;
            [DCMP(iN),gas,iG] = compexp_func (DCMP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ;
            
            % Now put through a heat exchanger - reject heat
            T_aim = T0;
            [gas,environ,iG,iE2] = hex_set(gas,[iL,iG],environ,[iL,iE2],T_aim,1.0,0.01);
        end
        
        
        for iN = 1:Ne_dis
            % HEAT in the hot packed bed
            %pbH = PB_TIMESTEP(pbH, gas, 'dis') ;
            pbH = PB_TIMESTEP_IDEAL(pbH, gas, iL, iG, 'dis') ;
            
            pbH = PB_FLUX(pbH, T0, p0, gas, iL) ;
            
            % Fluid outlet
            iG = iG + 1 ;
            gas.state(iL,iG).T = pbH.TF(pbH.NX,1) ;
            gas.state(iL,iG).p = pbH.P(pbH.NX,1) ;
            
            % Update fluid
            [gas] = update(gas,[iL,iG],1);
            
            % EXPAND
            PRe_dis = (gas.state(iL,iG).p/pbot)^(1/(Ne_dis+1-iN));  % stage expansion pressure ratio
            p_aim = gas.state(iL,iG).p/PRe_dis;
            [DEXP(iN),gas,iG] = compexp_func (DEXP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ;
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
        [~,DCMPtot]  = SUM_ENERGY(DCMP,DCMPtot,1,iL) ;
        [~,DEXPtot]  = SUM_ENERGY(DEXP,DEXPtot,1,iL) ;
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
    [DCMP,~]  = SUM_ENERGY(DCMP,DCMPtot,3,iL) ;
    [DEXP,~]  = SUM_ENERGY(DEXP,DEXPtot,3,iL) ;
    [gas.stage(iL,:),~] = SUM_ENERGY(gas.stage(iL,:),stagetot,4,numel(gas.stage(iL,:))) ;
    
    
    % Determine convergence and proceed
    D = [[gas.state(iL,:).T];[gas.state(iL,:).p]];

    if all(abs((D(D~=0) - D_0(D~=0))./D(D~=0))*100 < 1e-3) || counter > 100 % is discharge cycle converged?
        % Close working fluid streams
        gas.stage(iL,iG).type = 'end';
        gas = count_Nstg(gas);
        
        % Close air (heat rejection) streams
        iA_out = 1:2:(iA-1); iA_in  = iA_out + 1;
        for i=iA_in, air.stage(iL,i).type = 'end'; end
        air = count_Nstg(air);
                
        % Uncomment these lines to print states
        print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        
        [pbH,pbC] = PB_REVERSE(pbH,pbC,Load,iL) ;
        
        % Exit loop
        break
    else
        % Set new initial conditions
        if ~design_mode
            
            %print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
%              fprintf('p1:   %13.8f \n',gas.state(iL,1).p/1e5)
%              fprintf('pEND: %13.8f \n\n',gas.state(iL,iG).p/1e5)
%              fprintf('T1:   %13.8f \n',gas.state(iL,1).T)
%              fprintf('TEND: %13.8f \n\n',gas.state(iL,iG).T)
            % Adjust inlet pressure to try to reach convergence. The
            % 'smoothing' factor has to be quite small (<0.1, say) for this to be stable
            gas.state(iL,1).p = gas.state(iL,1).p - 0.10 * (gas.state(iL,iG).p - gas.state(iL,1).p) ;
            gas.state(iL,1).T = gas.state(iL,1).T + 0.10 * (gas.state(iL,iG).T - gas.state(iL,1).T) ;
            gas.state(iL,1).mdot = Load.mdot(iL);
            [gas] = update(gas,[iL,1],1);
            
        else
            gas.state(iL,1) = gas.state(iL,iG);
        end
        
        D_0 = D;
        counter = counter + 1 ;
    end
end

% Find t_dis
t_dis  = time ;

Load.time(iL) = min([Load.time(iL),t_dis])*(1-1e-6);

[pbH,pbC,iL,Icyc] = PB_STOP_CYC(pbH,pbC,Load,iL,Icyc,Ncyc,Lcyclic) ; 

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Atmospheric tanks
iA_out = 0; iA_in = 0;
%AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);