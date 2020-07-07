%%% CCES PLANT LAYOUT %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAES SIDE
%  ---------
%  1: Compressor inlet (atmosphere) - CCMP(1)
%  2: Hot HEX inlet                 - HX(1)
%  3: Medium HEX inlet              - HX(2)
%  4: Second Compressor inlet       - CCMP(2)
%  5: Second Medium HEX inlet       - HX(3)
%  6: Cold HEX inlet                - HX(4)
%  7: Coupler inlet (air side)      - HX(5)
%  8: Expander inlet
%  9: Liquid air tank inlet
%
%  PTES SIDE
%  ---------
%  1: Compressor inlet                  - CCMP(3)
%  2: Hot HEX inlet                     - HX(6)
%  3: Medium HEX inlet                  - HX(7)
%  4: Auxiliary compressor inlet        - CCMP(4)
%  5: Rejection unit inlet              - HX(8)
%  6: Regenerator inlet (high p side)   -
%  7: Expander inlet                    -
%  8: Coupler inlet (neon side)         -
%  9: Regenerator inlet (low p side)    -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set stage indices
iG1 = 1;  % keeps track of the gas stage number (LAES)
iG2 = 1;  % keeps track of the gas stage number (PTES)
iH  = 1;  % keeps track of the Hot fluid stream number (molten salt)
iM  = 1;  % keeps track of the Medium fluid stream number (mineral oil)
iC  = 1;  % keeps track of the Cold fluid stream number
iA  = 1;  % keeps track of the Air (heat rejection) stream number
iPMP = 1 ; % Keeps track of which pump is being used
iFAN = 1 ; % Keeps track of which fan is being used

% Initial guess of charge conditions (LAES side)
% Compressor inlet
gas1.state(iL,1).p = p0;
gas1.state(iL,1).T = T0;
gas1.state(iL,1).mdot = Load.mdot(iL);
[gas1] = update(gas1,[iL,1],1);

% Initial guess of charge conditions (PTES side). The mass flow rate is not
% know yet, as it is determined by the Coupler HEX.
% Compressor inlet (regenerator hot outlet)
gas2.state(iL,1).p = p0;
gas2.state(iL,1).T = T0;
[gas2] = update(gas2,[iL,1],1);
iCoup = 8;
gas2.state(iL,iCoup).p = p0;
gas2.state(iL,iCoup).T = 74;
[gas2] = update(gas2,[iL,iCoup],1);
% Set regenerator indices
iReg1 = iCoup-2; %high pressure side
iReg2 = iCoup+1; %low pressure side

% Set matrix of temperature and pressure points to test convergence
C_0 = [[gas1.state(iL,:).T];[gas1.state(iL,:).p]];
max_iter=100;
for counter=1:max_iter
    
    fprintf(1,['Charging CCES. Load period #',int2str(iL),'. Iteration #',int2str(counter),' \n'])
    
    %%%%%%%%%%%%%%%%
    %%%%% LAES %%%%%
    %%%%%%%%%%%%%%%%    
    % FIRST COMPRESSION
    T_aim = Tmax;
    [CCMP(1),gas1,iG1] = compexp_func (CCMP(1),iL,gas1,iG1,'Taim',T_aim, 1) ;
    
    % COOL (molten salt)
    % Set fluidH
    fluidH.state(iL,iH).T = HT.A(iL).T; fluidH.state(iL,iH).p = HT.A(iL).p;
    [fluidH] = update(fluidH,[iL,iH],1);
    % Run HEX
    [HX(1),gas1,iG1,fluidH,iH] = hex_func(HX(1),iL,gas1,iG1,fluidH,iH,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidH,iH] = compexp_func (CPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1); %#ok<*SAGROW>
    iH=iH+1;iPMP=iPMP+1;
    
    % COOL (mineral oil)
    % Set fluidM
    fluidM.state(iL,iM).T = MT.A(iL).T; fluidM.state(iL,iM).p = MT.A(iL).p;
    [fluidM] = update(fluidM,[iL,iM],1);
    % Run HEX
    [HX(2),gas1,iG1,fluidM,iM] = hex_func(HX(2),iL,gas1,iG1,fluidM,iM,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidM,iM] = compexp_func (CPMP(iPMP),iL,fluidM,iM,'Paim',fluidM.state(iL,1).p,1);
    iM=iM+1;iPMP=iPMP+1;
    
    % SECOND COMPRESSION
    p_aim = pmax_LA;
    [CCMP(2),gas1,iG1] = compexp_func (CCMP(2),iL,gas1,iG1,'Paim',p_aim, 1) ;
    
    % COOL (mineral oil)
    % Set fluidM
    fluidM.state(iL,iM).T = MT.A(iL).T; fluidM.state(iL,iM).p = MT.A(iL).p;
    [fluidM] = update(fluidM,[iL,iM],1);
    % Run HEX
    [HX(3),gas1,iG1,fluidM,iM] = hex_func(HX(3),iL,gas1,iG1,fluidM,iM,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidM,iM] = compexp_func (CPMP(iPMP),iL,fluidM,iM,'Paim',fluidM.state(iL,1).p,1);
    iM=iM+1;iPMP=iPMP+1;
    
    % COOL (cold fluid)
    % Set fluidC
    fluidC.state(iL,iC).T = CT.A(iL).T; fluidC.state(iL,iC).p = CT.A(iL).p;
    [fluidC] = update(fluidC,[iL,iC],1);
    % Run HEX
    [HX(4),gas1,iG1,fluidC,iC] = hex_func(HX(4),iL,gas1,iG1,fluidC,iC,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidC,iC] = compexp_func (CPMP(iPMP),iL,fluidC,iC,'Paim',fluidC.state(iL,1).p,1);
    iC=iC+1;iPMP=iPMP+1;
    
    % COUPLER HEX (air-neon)
    [HX(5),gas1,iG1,gas2,~] = hex_func(HX(5),iL,gas1,iG1,gas2,iCoup,1,1.0);
    % Update mass flow rate at the PTES side
    gas2.state(iL,1).mdot = gas2.state(iL,iCoup).mdot;
    
    % EXPANSION
    p_aim = max([RPN('QT_INPUTS',0.0,gas1.state(iL,iG1).T,'P',gas1),p0]);
    [CEXP(1),gas1,iG1] = compexp_func (CEXP(1),iL,gas1,iG1,'Paim',p_aim, 1) ;
    
    %%%%%%%%%%%%%%%%
    %%%%% PTES %%%%%
    %%%%%%%%%%%%%%%%
    % FIRST COMPRESSION
    T_aim = Tmax;
    [CCMP(3),gas2,iG2] = compexp_func (CCMP(3),iL,gas2,iG2,'Taim',T_aim, 1) ;
    
    % COOL (molten salt)
    % Set fluidH
    fluidH.state(iL,iH).T = HT.A(iL).T; fluidH.state(iL,iH).p = HT.A(iL).p;
    [fluidH] = update(fluidH,[iL,iH],1);
    % Run HEX
    [HX(6),gas2,iG2,fluidH,iH] = hex_func(HX(6),iL,gas2,iG2,fluidH,iH,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidH,iH] = compexp_func (CPMP(iPMP),iL,fluidH,iH,'Paim',fluidH.state(iL,1).p,1);
    iH=iH+1;iPMP=iPMP+1;
    
    % COOL (mineral oil)
    % Set fluidM
    fluidM.state(iL,iM).T = MT.A(iL).T; fluidM.state(iL,iM).p = MT.A(iL).p;
    [fluidM] = update(fluidM,[iL,iM],1);
    % Run HEX
    [HX(7),gas2,iG2,fluidM,iM] = hex_func(HX(7),iL,gas2,iG2,fluidM,iM,1,1.0);
    % Run Pump
    [CPMP(iPMP),fluidM,iM] = compexp_func (CPMP(iPMP),iL,fluidM,iM,'Paim',fluidM.state(iL,1).p,1);
    iM=iM+1;iPMP=iPMP+1;
    
    % SECOND COMPRESSION
    p_aim = gas2.state(iL,iG2).p*1.5;
    [CCMP(2),gas2,iG2] = compexp_func (CCMP(2),iL,gas2,iG2,'Paim',p_aim, 1) ;
    
    % REJECT HEAT (external HEX)
    T_aim = environ.T0 + T0_inc;    
    air.state(iL,iA).T = T0; air.state(iL,iA).p = p0; air = update(air,[iL,iA],1);
    [HX(8), gas2, iG2, air, iA] = hex_func(HX(8),iL,gas2,iG2,air,iA,5,T_aim);
    [CFAN(iFAN),air,iA] = compexp_func (CFAN(1),iL,air,iA,'Paim',p0,1);
    iFAN=iFAN+1;
    
    % REGENERATE (gas-gas)
    [HX(9),gas2,iG2,~,~] = hex_func(HX(9),iL,gas2,iG2,gas2,iCoup+1,0,0);
    
    % EXPANSION
    p_aim = p0;
    [CEXP(2),gas2,iG2] = compexp_func (CEXP(2),iL,gas2,iG2,'Paim',p_aim, 1) ;
    
    % COUPLER (this has already been done)
    iG2 = iG2+1;
    
    % REGENERATE (gas-gas)
    [HX(9),~,~,gas2,iG2] = hex_func(HX(9),iL,gas2,iReg1,gas2,iG2,0,0);
    
    print_states(gas1,iL,1:12,Load);
    print_states(gas2,iL,1:12,Load);
    keyboard
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
    %{
    
    % Determine convergence and proceed
    C = [[gas.state(iL,:).T];[gas.state(iL,:).p]];
    convergence = all(abs((C(C~=0) - C_0(C~=0))./C(C~=0))*100 < 1e-3);
    
    if (convergence && strcmp(HX_model_temp,HX_model)) || counter==max_iter % is charge cycle converged?

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
        
        if Load.mode==3
            % Close HTF streams
            iHTF_out = 1:4:(iHTF-1); iHTF_in  = iHTF_out + 3;
            for i=iHTF_in, HTF.stage(iL,i).type = 'end'; end
            HTF = count_Nstg(HTF);
        end
        
        % Uncomment these lines to print states
        %{
        print_states(gas,iL,1:gas.Nstg(iL)+1,Load);
        print_states(fluidH,iL,1:fluidH.Nstg(iL)+1,Load);
        print_states(fluidC,iL,1:fluidC.Nstg(iL)+1,Load);
        print_states(air,iL,1:air.Nstg(iL)+1,Load);
        print_states(HTF,iL,1:HTF.Nstg(iL)+1,Load);
        keyboard
        %}
        
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
        C_0 = C;
        iG=1; iH=1; iC=1; iA=1; iPMP=1; iHTF=1;
        
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
        iG=1; iH=1; iC=1; iA=1; iPMP=1; iHTF=1;
        
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
%}