% Set stage indices
iG = 1;  % keeps track of the gas stage number
iH = 1;  % keeps track of the Hot fluid stream number
iC = 1;  % keeps track of the Cold fluid stream number
iE = 1;  % keeps track of the heat rejection stream number
iA = 1;  % keeps track of the Air (heat rejection) stream number
iPMP = 1 ; % Keeps track of which pump is being used

% Regenerator inlet indices
iReg1 = 1 + Nc_ch*2;         % index hot regenerator inlet (after compression + cooling)
iReg2 = iReg1 + 2 + Ne_ch*2; % index cold regenerator inlet (after regeneration + heat_reject + expansion + heating)

if design_mode == 1
    % Initial guess of charge conditions
    % Compressor inlet (regenerator hot outlet)
    gas.state(iL,1).p    = pbot; gas.state(iL,1).T = T1;
    gas.state(iL,1).mdot = Load.mdot(iL);
    [gas] = update(gas,[iL,1],1);
    
    gas.state(iL,iReg2).T = T0;
    gas.state(iL,iReg2).p = pbot;
    gas.state(iL,iReg2).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,iReg2],1);
else
    % Set up matrix that guesses the converged solution
    for ii = 1 : numel(C(C~=0))/2
        gas.state(iL,ii).T    = C(1,ii) ;
        gas.state(iL,ii).p    = C(2,ii) ;
        gas.state(iL,ii).mdot = Load.mdot(iL) ;
        
        % For inventory control, assume that the pressure scales with the off-design mass flow rate
        %gas.state(iL,ii).p = CCMP.Pin * Load.mdot(iL) / CCMP.mdot0 ;
        gas.state(iL,ii).p = gas.state(iL,ii).p * Load.mdot(iL) / CCMP.mdot0 ;
        [gas] = update(gas,[iL,ii],1);
    end  
end

% Set matrix of temperature and pressure points to test convergence
C_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
max_iter=100;
for counter=1:max_iter
    
    fprintf(1,['Charging JB PTES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
        
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
        
        % COOL (gas-liquid)
        fluidH.state(iL,iH).T = HT.A(iL).T; fluidH.state(iL,iH).p = HT.A(iL).p;
        [fluidH] = update(fluidH,[iL,iH],1);

        [HX(ihx_hot(iN)),gas,iG,fluidH,iH] = hex_func(HX(ihx_hot(iN)),iL,gas,iG,fluidH,iH,1,1.0);
        % Now calculate pump requirements for moving fluidH
        [CPMP(iPMP),fluidH,iH] = compexp_func (CPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
        iH=iH+1;iPMP=iPMP+1;

    end    
        
    % REGENERATE (gas-gas)
    [HX(ihx_reg),gas,iG,~,~] = hex_func(HX(ihx_reg),iL,gas,iReg1,gas,iReg2,0,0);
        
    % REJECT HEAT (external HEX)
    T_aim = environ.T0 + T0_inc;    
    air.state(iL,iA).T = T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
    [HX(ihx_rej), gas, iG, air, iA] = hex_func(HX(ihx_rej),iL,gas,iG,air,iA,5,T_aim);
    [CFAN(1),air,iA] = compexp_func (CFAN(1),iL,air,iA,'Paim',p0,1);
        
    for iN = 1:Ne_ch
        % EXPAND
        PRe = (gas.state(iL,iG).p/pbot)^(1/(Ne_ch+1-iN)); % stage expansion pressure ratio
        p_aim = gas.state(iL,iG).p/PRe;
        [CEXP(iN),gas,iG] = compexp_func (CEXP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ; 
        
        % HEAT (gas-liquid)
        fluidC.state(iL,iC).T = CT.A(iL).T; fluidC.state(iL,iC).p = CT.A(iL).p;
        [fluidC] = update(fluidC,[iL,iC],1);

        [HX(ihx_cld(iN)),fluidC,iC,gas,iG] = hex_func(HX(ihx_cld(iN)),iL,fluidC,iC,gas,iG,2,1.0);
        [CPMP(iPMP),fluidC,iC] = compexp_func (CPMP(iPMP),iL,fluidC,iC,'Paim',fluidC.state(iL,1).p,1);
        iC=iC+1; iPMP=iPMP+1;

    end
    
    % REGENERATE (gas-gas)
    [HX(ihx_reg),~,~,gas,iG] = hex_func(HX(ihx_reg),iL,gas,iReg1,gas,iReg2,0,0);
    
    % Determine convergence and proceed
    C = [[gas.state(iL,:).T];[gas.state(iL,:).p]];

    if (all(abs((C(C~=0) - C_0(C~=0))./C(C~=0))*100 < 1e-3) || counter==max_iter) % is charge cycle converged?

        % Close working fluid streams
        gas.stage(iL,iG).type = 'end';
        gas = count_Nstg(gas);
        
        % Close air (heat rejection) streams
        iA_out = 1:3:(iA-1); iA_in  = iA_out + 2;
        for i=iA_in, air.stage(iL,i).type = 'end'; end
        air = count_Nstg(air);
        
        % Close storage fluid streams
        iH_out = 1:3:(iH-1); iH_in  = iH_out + 2;
        iC_out = 1:3:(iC-1); iC_in  = iC_out + 2;
        for i=iH_in, fluidH.stage(iL,i).type = 'end'; end
        for i=iC_in, fluidC.stage(iL,i).type = 'end'; end
        fluidH = count_Nstg(fluidH);
        fluidC = count_Nstg(fluidC);
        
        % Uncomment these lines to print states
        %%{
        print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
        print_states(fluidC,iL,1:fluidC.Nstg(iL)+1,Load);
        print_states(air,iL,1:air.Nstg(iL)+1,Load);
        %keyboard
        %}
        
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
            % Reduce smoothing factor with number of iterations
            smooth = 0.10;% / double(counter)^0.2 ; 
            gas.state(iL,1).p = gas.state(iL,1).p - smooth * (gas.state(iL,iG).p - gas.state(iL,1).p) ;
            gas.state(iL,1).T = gas.state(iL,1).T + smooth * (gas.state(iL,iG).T - gas.state(iL,1).T) ;
                      
            gas.state(iL,1).mdot = Load.mdot(iL);
            [gas] = update(gas,[iL,1],1);
            
        else
            gas.state(iL,1) = gas.state(iL,iG); 
        end
        
        C_0 = C;
        iG=1; iH=1; iC=1; iE=1; iA=1; iPMP=1;
        
    end
end
if counter==max_iter
    warning('Exiting JB_CHARGE cycle without having reached convergence');
end

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
HT = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
% Cold tanks
CT = run_tanks(CT,iL,fluidC,iC_out,iC_in,Load,T0);
% Atmospheric tanks
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);