% Regenerator cold inlet
ind.iReg1 = 1; % index regenerator hot inlet
switch Load.mode
    case 0 % PTES
        ind.iReg2 = ind.iReg1 + 1 + 3*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + cooling + compression)    
    case 2 % Heat engine only
        ind.iReg2 = ind.iReg1 + 1 + 2*Nc_dis; % index regenerator cold inlet (after regeneration + heat rejection + compression)    
end

% Other indices
ind.ihx_rejd = ihx_rejd;
ind.Nc_dis   = Nc_dis ;
ind.Ne_dis   = Ne_dis ;

% Compute PR_dis based on charge pressure ratio and PRr
TP.PRdis = PRr*PRch;
TP.ploss = ploss ;

if design_mode == 1
    % Initial guess of discharge conditions
    % Expander outlet (regenerator hot inlet)
    gas.state(iL,1).p    = pbot; gas.state(iL,1).T = T1;
    gas.state(iL,1).mdot = Load.mdot(iL);
    [gas] = update(gas,[iL,1],1);
    gas.state(iL,ind.iReg2).T = T0;
    gas.state(iL,ind.iReg2).p = gas.state(iL,1).p*TP.PRdis;
    gas.state(iL,ind.iReg2).mdot = gas.state(iL,1).mdot;
    [gas] = update(gas,[iL,ind.iReg2],1);
    TOLconv = 1e-3 ;
    environ.T0 = T0 ;

    switch HX_model
        case 'eff'
            HX_model_temp = 'eff';
        case 'geom'
            % Set the heat exchanger models to 'eff' temporarily, to obtain
            % approximated cycle points before using the 'geom' model
            HX_model_temp = 'eff';
            for ihx=1:length(HX)
                HX(ihx).model = HX_model_temp;
            end
        otherwise
            error('not implemented')
    end

    mprev  = 0 ; erprevM = 0 ; gradMM = 0;

else
    for ii = 1 : numel(D(D~=0))/2
        gas.state(iL,ii).T    = D(1,ii) ;
        gas.state(iL,ii).p    = D(2,ii) ;
        gas.state(iL,ii).mdot = Load.mdot(iL) ;
        
        % For inventory control, assume that the pressure scales with the off-design mass flow rate
        %gas.state(iL,ii).p = (DEXP.Pin/DEXP.pr0) * Load.mdot(iL) / DEXP.mdot0 ;
        i_dis = Design_Load.ind(any(Design_Load.type == {'dis'},2));
        gas.state(iL,ii).p = gas0.state(i_dis,ii).p * (Load.mdot(iL) / DEXP.mdot0) *  sqrt(Load.T0_off(iL) / T0) ; % Second ever run is discharging
        [gas] = update(gas,[iL,ii],1);
        
    end   
    
    p1prev = 0 ; erprevP = 0 ; gradPP = 0 ;
    T1prev = 0 ; erprevT = 0 ; gradTT = 0;

    environ.T0 = Load.T0_off(iL) ;
    TOLconv = 1e-4 ;
end


% Set matrix of temperature and pressure points to test convergence
D_0 = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
max_iter = 150;
for counter = 1:max_iter
    fprintf(1,['Discharging JB PTES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])

    [gas,fluidH,fluidC,HT,CT,air,DCMP,DEXP,DPMP,DFAN,HX,iG,iH,iC,iA] = ...
        run_JB_discharge_alt_Qrej(ind,gas,gas0,fluidH,fluidC,HT,CT,air,environ,DCMP,DEXP,DPMP,DFAN,HX,HX0,TP,Load,design_mode,iL);
    
    % Determine convergence and proceed
    D = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    convergence = all(abs((D(D~=0) - D_0(D~=0))./D(D~=0))*100 < TOLconv);
    
    if (convergence && strcmp(HX_model_temp,HX_model)) || counter >= max_iter % is discharge cycle converged?
        
        % Exit loop
        break
        
    elseif convergence
        % If convergence has been reach but HX_model_temp~=HX_model, set
        % the heat exchanger models back to the original using the new
        % converged state, and resume iteration
        HX_model_temp = HX_model;
        for ihx=1:length(HX)
            HX(ihx).model = HX_model_temp;
        end
        gas.state(iL,1) = gas.state(iL,iG);        
        D_0 = D;
        
    else
        % Set new initial conditions
        if ~design_mode
            
            %print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
%              fprintf('p1:   %13.8f \n',gas.state(iL,1).p/1e5)
%              fprintf('pEND: %13.8f \n\n',gas.state(iL,iG).p/1e5)
%              fprintf('T1:   %13.8f \n',gas.state(iL,1).T)
%              fprintf('TEND: %13.8f \n\n',gas.state(iL,iG).T)


            ernewP = gas.state(iL,1).p  - gas.state(iL,iG).p ;
            ernewT = gas.state(iL,1).T  - gas.state(iL,iG).T ;
                
            if counter == 1
                
                p1prev = gas.state(iL,1).p ;
                T1prev = gas.state(iL,1).T ;
                smooth = 0.025 ;
                
                gas.state(iL,1).p = gas.state(iL,1).p - smooth * (gas.state(iL,iG).p - gas.state(iL,1).p) ;
                gas.state(iL,1).T = gas.state(iL,1).T + smooth * (gas.state(iL,iG).T - gas.state(iL,1).T) ;
            else
                gradPP  = (ernewP - erprevP) / (gas.state(iL,1).p - p1prev) ;
                gradTT  = (ernewT - erprevT) / (gas.state(iL,1).T - T1prev) ;
                
                p1prev = gas.state(iL,1).p ;
                T1prev = gas.state(iL,1).T ;

                gas.state(iL,1).p = gas.state(iL,1).p - 0.2 * ernewP / gradPP;
                gas.state(iL,1).T = gas.state(iL,1).T - 0.2 * ernewT / gradTT;
            end
            
            erprevP = ernewP ;
            erprevT = ernewT ;
            [gas] = update(gas,[iL,1],1);
            
        else
            gas.state(iL,1) = gas.state(iL,iG);

            ernewM = fluidH.state(1,1).mdot * Design_Load.time(1) - fluidH.state(iL,1).mdot * Design_Load.time(iL);
                
            if counter == 1
                mprev = fluidH.state(iL,1).mdot ;
                smooth = 0.025 ;                
            else
                gradMM  = (ernewM - erprevM) / ((fluidH.state(iL,1).mdot - mprev));
                mprev = fluidH.state(iL,1).mdot ;
                TP.PRdis = TP.PRdis - 0.05 * ernewM / gradMM ;
            end
            
            erprevM = ernewM ;

        end
        
        D_0 = D;
        
    end
end
if counter==max_iter
    warning('Exiting JB_DISCHARGE cycle without having reached convergence');
end

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

print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
print_states(fluidC,iL,1:fluidC.Nstg(iL)+1,Load);
print_states(air,iL,1:air.Nstg(iL)+1,Load);


% Compute temperatures and mass of the sink tanks. Assume complete
% depletion of (at least) one of the source tanks to stablish
% discharge time.

% Find t_dis
MdotH = total_mdot(fluidH,iL,iH_out);
t_dis  = HT.B(iL).M/MdotH;
if Load.mode == 0
    MdotC = total_mdot(fluidC,iL,iC_out);
    tC_dis  = CT.B(iL).M/MdotC;
    t_dis   = min([t_dis,tC_dis]);
end
Load.time(iL) = min([Load.time(iL),t_dis])*(1-1e-6);

% Compute effect of fluid streams entering/leaving the sink/source tanks
% Hot tanks
HT = run_tanks(HT,iL,fluidH,iH_out,iH_in,Load,T0);
if Load.mode == 0
    % Cold tanks
    CT = run_tanks(CT,iL,fluidC,iC_out,iC_in,Load,T0);
end
% Atmospheric tanks
AT = run_tanks(AT,iL,air,iA_out,iA_in,Load,T0);


function [gas,fluidH,fluidC,HT,CT,air,DCMP,DEXP,DPMP,DFAN,HX,iG,iH,iC,iA] = run_JB_discharge(ind,gas,gas0,fluidH,fluidC,HT,CT,air,environ,DCMP,DEXP,DPMP,DFAN,HX,HX0,TP,Load,design_mode,iL)

    % Set stage indices
    iG = 1;  % keeps track of the gas stage number
    iH = 1;  % keeps track of the Hot fluid stream number
    iC = 1;  % keeps track of the Cold fluid stream number
    iE = 1;  % keeps track of the Environment (heat rejection) stream number
    iA = 1;  % keeps track of the Air (heat rejection) stream number
    iPMP = 1 ; % Keeps track of which pump is being used

    CSmode = 1 ; % Determines how cold store is discharged

    % REGENERATE (gas-gas)
    [HX(ind.ihx_reg),gas,iG,~,~] = hex_func(HX(ind.ihx_reg),iL,gas,ind.iReg1,gas,ind.iReg2,0,0);
    
    for iN = 1:ind.Nc_dis
        % REJECT HEAT (external HEX)
        %T_aim = environ.T0 + T0_inc;
        air.state(iL,iA).T = environ.T0; air.state(iL,iA).p = TP.p0; air = update(air,[iL,iA],1);
        
        if design_mode == 0 && Load.T0_off(iL) < TP.T0
            % If off-design and T0 is colder than the design value
            % Reduce air mass flow rate so that expander inlet temperature
            % remains at the design point.
            [HX(ind.ihx_rejd(iN)), gas, iG, air, iA] = hex_func(HX(ind.ihx_rejd(iN)),iL,gas,iG,air,iA,5,gas0.state(iL,iG+1).T);
        else
            [HX(ind.ihx_rejd(iN)), gas, iG, air, iA] = hex_func(HX(ind.ihx_rejd(iN)),iL,gas,iG,air,iA,1,0.5);
        end
        [DFAN(1),air,iA] = compexp_func (DFAN(1),iL,air,iA,'Paim',TP.p0,1);
        
        switch Load.mode
            case 0 % PTES
                % COOL (gas-liquid)
                fluidC.state(iL,iC).T = CT.B(iL).T; fluidC.state(iL,iC).p = CT.B(iL).p; %#ok<*SAGROW>
                [fluidC] = update(fluidC,[iL,iC],1);

                if CSmode == 0
                    % Equal mdot*cp on each side of heat exchanger
                    [HX(ind.ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ind.ihx_cld(iN)),iL,gas,iG,fluidC,iC,1,1.0);
                elseif CSmode == 1
                    % Force CS fluid to return to intial temperature
                    [HX(ind.ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ind.ihx_cld(iN)),iL,gas,iG,fluidC,iC,3,fluidC.state(1,1).T);
                elseif CSmode == 2 && design_mode == 0
                    % Discharge the CS at the same rate as the hot store so
                    % they discharge in the same amount of time
                    if isempty(HX(ind.ihx_hot(iN)).H(iL).mdot)
                        HSrat = 1;
                    else
                        HSrat = HX(ind.ihx_hot(iN)).H(iL).mdot / HX0(ind.ihx_hot(iN)).H(2).mdot ;
                    end
                    if HSrat <= 0 || isnan(HSrat)
                        HSrat = 1;
                    end
                    fluidC.state(iL,iC).mdot = HSrat * HX0(ind.ihx_cld(iN)).C(2).mdot ;
                    [HX(ind.ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ind.ihx_cld(iN)),iL,gas,iG,fluidC,iC,0,-1);

                elseif CSmode == 2 && design_mode == 1
                    % Discharge the CS at the same rate as the hot store so
                    % they discharge in the same amount of time
                    if isempty(HX(ind.ihx_hot(iN)).H(iL).mdot)
                        HSrat = 1;
                    else
                        HSrat = HX(ind.ihx_hot(iN)).H(iL).mdot / HX(ind.ihx_hot(iN)).C(1).mdot ;
                    end
                    if HSrat <= 0 || isnan(HSrat)
                        HSrat = 1;
                    end
                    fluidC.state(iL,iC).mdot = HSrat * HX(ind.ihx_cld(iN)).H(1).mdot ;
                    [HX(ind.ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ind.ihx_cld(iN)),iL,gas,iG,fluidC,iC,0,-1);
                
                end

                [DPMP(iPMP),fluidC,iC] = compexp_func (DPMP(iPMP),iL,fluidC,iC,'Paim',fluidC.state(iL,1).p,1);
                iC=iC+1; iPMP=iPMP+1;
            case 1 % Heat engine only
        end
        
        % COMPRESS
        PRc_dis = (TP.PRdis)^(1/ind.Nc_dis)/(1-TP.ploss)^2; % stage compression pressure ratio
        p_aim = gas.state(iL,iG).p*PRc_dis;
        [DCMP(iN),gas,iG] = compexp_func (DCMP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ; 
    end
    
    % REGENERATE (gas-gas)
    [HX(ind.ihx_reg),~,~,gas,iG] = hex_func(HX(ind.ihx_reg),iL,gas,ind.iReg1,gas,ind.iReg2,0,0);
    
    for iN = 1:ind.Ne_dis
        % HEAT (gas-fluid)
        fluidH.state(iL,iH).T = HT.B(iL).T; fluidH.state(iL,iH).p = HT.B(iL).p; THmin = HT.A(1).T;
        [fluidH] = update(fluidH,[iL,iH],1);
        Taim = THmin;

        [HX(ind.ihx_hot(iN)),fluidH,iH,gas,iG] = hex_func(HX(ind.ihx_hot(iN)),iL,fluidH,iH,gas,iG,4,THmin);
        [DPMP(iPMP),fluidH,iH] = compexp_func (DPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
        iH=iH+1; iPMP=iPMP+1;
        
        % EXPAND
        PRe_dis = (gas.state(iL,iG).p/TP.pbot)^(1/(ind.Ne_dis+1-iN));  % stage expansion pressure ratio
        p_aim = gas.state(iL,iG).p/PRe_dis;
        [DEXP(iN),gas,iG] = compexp_func (DEXP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ; 
    end
    

end


function [gas,fluidH,fluidC,HT,CT,air,DCMP,DEXP,DPMP,DFAN,HX,iG,iH,iC,iA] = run_JB_discharge_alt_Qrej(ind,gas,gas0,fluidH,fluidC,HT,CT,air,environ,DCMP,DEXP,DPMP,DFAN,HX,HX0,TP,Load,design_mode,iL)

    % Set stage indices
    iG = 1;  % keeps track of the gas stage number
    iH = 1;  % keeps track of the Hot fluid stream number
    iC = 1;  % keeps track of the Cold fluid stream number
    iE = 1;  % keeps track of the Environment (heat rejection) stream number
    iA = 1;  % keeps track of the Air (heat rejection) stream number
    iPMP = 1 ; % Keeps track of which pump is being used

    CSmode = 1 ; % Determines how cold store is discharged

    % REGENERATE (gas-gas)
    [HX(ind.ihx_reg),gas,iG,~,~] = hex_func(HX(ind.ihx_reg),iL,gas,ind.iReg1,gas,ind.iReg2,0,0);
    
    for iN = 1:ind.Nc_dis

        if ind.Nc_dis>1
            warning('Have not tested Nc_dis > 1.');
        end
        switch Load.mode
            case 0 % PTES
                % COOL (gas-liquid)
                fluidC.state(iL,iC).T = CT.B(iL).T; fluidC.state(iL,iC).p = CT.B(iL).p; %#ok<*SAGROW>
                [fluidC] = update(fluidC,[iL,iC],1);

                %if CSmode == 2 && design_mode == 1
                %    warning('CSmode == 2 only works in off-design operation. Changing to CSmode == 1')
                %    CSmode = 1 ;
                %end

                if CSmode == 0
                    % Equal mdot*cp on each side of heat exchanger
                    [HX(ind.ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ind.ihx_cld(iN)),iL,gas,iG,fluidC,iC,1,1.0);
                elseif CSmode == 1
                    % Force CS fluid to return to intial temperature
                    [HX(ind.ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ind.ihx_cld(iN)),iL,gas,iG,fluidC,iC,3,fluidC.state(1,1).T);
                elseif CSmode == 2 && design_mode == 0
                    % Discharge the CS at the same rate as the hot store so
                    % they discharge in the same amount of time
                    if isempty(HX(ind.ihx_hot(iN)).H(iL).mdot) || isempty(HX0(ind.ihx_hot(iN)).H(2).mdot)
                        HSrat = 1;
                    else
                        HSrat = HX(ind.ihx_hot(iN)).H(iL).mdot / HX0(ind.ihx_hot(iN)).H(2).mdot ;
                    end
                    if HSrat <= 0 || isnan(HSrat)
                        HSrat = 1;
                    end
                    fluidC.state(iL,iC).mdot = HSrat * HX0(ind.ihx_cld(iN)).C(2).mdot ;
                    [HX(ind.ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ind.ihx_cld(iN)),iL,gas,iG,fluidC,iC,0,-1);
                elseif CSmode == 2 && design_mode == 1
                    % Discharge the CS at the same rate as the hot store so
                    % they discharge in the same amount of time
                    if isempty(HX(ind.ihx_hot(iN)).H(iL).mdot)
                        HSrat = 1;
                    else
                        HSrat = HX(ind.ihx_hot(iN)).H(iL).mdot / HX(ind.ihx_hot(iN)).C(1).mdot ;
                    end
                    if HSrat <= 0 || isnan(HSrat)
                        HSrat = 1;
                    end
                    fluidC.state(iL,iC).mdot = HSrat * HX(ind.ihx_cld(iN)).H(1).mdot ;
                    [HX(ind.ihx_cld(iN)),gas,iG,fluidC,iC] = hex_func(HX(ind.ihx_cld(iN)),iL,gas,iG,fluidC,iC,0,-1);
                
                end

                [DPMP(iPMP),fluidC,iC] = compexp_func (DPMP(iPMP),iL,fluidC,iC,'Paim',fluidC.state(iL,1).p,1);
                iC=iC+1; iPMP=iPMP+1;
            case 1 % Heat engine only
        end
        
        % COMPRESS
        PRc_dis = (TP.PRdis)^(1/ind.Nc_dis)/(1-TP.ploss)^2; % stage compression pressure ratio
        p_aim = gas.state(iL,iG).p*PRc_dis;
        [DCMP(iN),gas,iG] = compexp_func (DCMP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ; 

        % REJECT HEAT (external HEX)
        air.state(iL,iA).T = environ.T0; air.state(iL,iA).p = TP.p0; air = update(air,[iL,iA],1);
        
        if design_mode == 0 && Load.T0_off(iL) < TP.T0
            % If off-design and T0 is colder than the design value
            % Reduce air mass flow rate so that expander inlet temperature
            % remains at the design point.
            [HX(ind.ihx_rejd(iN)), gas, iG, air, iA] = hex_func(HX(ind.ihx_rejd(iN)),iL,gas,iG,air,iA,5,gas0.state(iL,iG+1).T);
        else
            [HX(ind.ihx_rejd(iN)), gas, iG, air, iA] = hex_func(HX(ind.ihx_rejd(iN)),iL,gas,iG,air,iA,1,0.5);
        end
        [DFAN(1),air,iA] = compexp_func (DFAN(1),iL,air,iA,'Paim',TP.p0,1);
    end
    
    % REGENERATE (gas-gas)
    [HX(ind.ihx_reg),~,~,gas,iG] = hex_func(HX(ind.ihx_reg),iL,gas,ind.iReg1,gas,ind.iReg2,0,0);
    
    for iN = 1:ind.Ne_dis
        % HEAT (gas-fluid)
        fluidH.state(iL,iH).T = HT.B(iL).T; fluidH.state(iL,iH).p = HT.B(iL).p; THmin = HT.A(1).T;
        [fluidH] = update(fluidH,[iL,iH],1);
        
        [HX(ind.ihx_hot(iN)),fluidH,iH,gas,iG] = hex_func(HX(ind.ihx_hot(iN)),iL,fluidH,iH,gas,iG,4,THmin);
        [DPMP(iPMP),fluidH,iH] = compexp_func (DPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
        iH=iH+1; iPMP=iPMP+1;
        
        % EXPAND
        PRe_dis = (gas.state(iL,iG).p/TP.pbot)^(1/(ind.Ne_dis+1-iN));  % stage expansion pressure ratio
        p_aim = gas.state(iL,iG).p/PRe_dis;
        [DEXP(iN),gas,iG] = compexp_func (DEXP(iN),iL,gas,iG,'Paim',p_aim, design_mode) ; 
    end
    
end
