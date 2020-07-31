% Charging function
function [obj, TS, TF, iCYC] = PB_RUN (obj, Npb, fld, iCYC, mode)

Lend = false;
i    = 1 ;
ii   = 1 ;
time = 0 ;

switch mode
    case 'chg'
        if obj(Npb).Ltext
            Tx   = obj(Npb).TD + obj(Npb).textC * (obj(Npb).TC - obj(Npb).TD) ;
        end
        str1 = 'CHARGING: TIMESTEP # %8i\n' ;
        str2 = 'CHARGING COMPLETE\n\n' ;
        Tend = obj(Npb).timeC ;
        ff    = 1.0 ;
    case 'dis'
        if obj(Npb).Ltext
            Tx = obj(Npb).TC + obj(Npb).textD * (obj(Npb).TD - obj(Npb).TC) ;
        end
        str1 = 'DISCHARGING: TIMESTEP # %8i\n' ;
        str2 = 'DISCHARGING COMPLETE\n\n' ;
        Tend = obj(Npb).timeD ;
        ff    = -1.0 ;
end

if obj.TC > obj.TD
    T0 = obj.TD ;
else
    T0 = obj.TC ;
end
% Reference points
HF0 = RPN('PT_INPUTS',1e5,T0,'H',fld) ;
SF0 = RPN('PT_INPUTS',1e5,T0,'S',fld) ;

while ~Lend
    
    if mod(i,100) == 0
        fprintf(1,str1,i) ;
    end
    
    % Step forward one time step
    for j = 1 : Npb
        obj(Npb) = PB_TIMESTEP(obj(Npb), fld, mode) ;
        %obj(Npb) = PB_TIMESTEP_IDEAL(obj(Npb), fld, mode) ;
        
        % >> Move this stuff into PB_TIMESTEP?
        % Calculate the mass, enthalpy, and entropy flux INTO the storage in that timestep
        Min = obj(Npb).u(1,1) * obj(Npb).rho(1,1) * obj(Npb).A * obj(Npb).dt ;
        obj(Npb).Mflux(iCYC, 1) = obj(Npb).Mflux(iCYC, 1) + Min ;
        obj(Npb).Hflux(iCYC, 1) = obj(Npb).Hflux(iCYC, 1) + Min * (RPN('PT_INPUTS',obj.P(1),obj.TF(1,1),'H',fld) - HF0);
        %obj(Npb).Hflux(iCYC, 1) = obj(Npb).Hflux(iCYC, 1) + Min * (obj(Npb).cF *(obj.TF(1,1)-T0));
        obj(Npb).Sflux(iCYC, 1) = obj(Npb).Sflux(iCYC, 1) + Min * (RPN('PT_INPUTS',obj.P(1),obj.TF(1,1),'S',fld) - SF0);
    
        % Calculate the mass, enthalpy, and entropy flux OUT OF the storage in that timestep
        Mout = obj(Npb).u(end,1) * obj(Npb).rho(end,1) * obj(Npb).A * obj(Npb).dt ;
        obj(Npb).Mflux(iCYC, 2) = obj(Npb).Mflux(iCYC, 2) + Mout ;
        obj(Npb).Hflux(iCYC, 2) = obj(Npb).Hflux(iCYC, 2) + Mout * (RPN('PT_INPUTS',obj.P(end),obj.TF(end,1),'H',fld) - HF0) ;
        %obj(Npb).Hflux(iCYC, 2) = obj(Npb).Hflux(iCYC, 2) + Mout * (obj(Npb).cF *(obj.TF(end,1)-T0));
        obj(Npb).Sflux(iCYC, 2) = obj(Npb).Sflux(iCYC, 2) + Mout * (RPN('PT_INPUTS',obj.P(end),obj.TF(end,1),'S',fld) - SF0) ;
            
    end
    
    
    % Check whether to end charging
    if obj(Npb).Ltime
        if time > obj(Npb).TMAX || time > Tend
            Lend = true ;
        end
    elseif obj(Npb).Ltext
        if time > obj(Npb).TMAX || ff*obj(Npb).TS(end,1) > ff*Tx
            Lend = true ;
        end
    end
    
    % Plot out results if appropriate
    if (i == obj.Tprof(ii)) || Lend
        TS(:,ii) = obj.TS(:,1) ;
        TF(:,ii) = obj.TF(:,1) ;
        ii = ii + 1 ;
    end
    
    i = i + 1 ;
    time = time + obj.dt ;
end

for j = 1 : Npb
    [obj(Npb), obj(Npb).H(iCYC), obj(Npb).S(iCYC)] = PB_ENERGY(obj(Npb), fld) ; % Evaluate energy in packed bed at end of charge
    
    % Calculate lossy stuff
    obj(Npb).W(iCYC)    = 0.0 ;
    obj(Npb).Q(iCYC)    = obj(Npb).Hflux(iCYC,1) - obj(Npb).Hflux(iCYC,2) ; % Heat transferred into/out of storage
    obj(Npb).DH(iCYC)   = obj(Npb).H(iCYC) - obj(Npb).H(iCYC-1) ; % Change in enthalpy of packed bed. Should equal Q (small errors occur though)
    obj(Npb).Sirr(iCYC) = obj(Npb).S(iCYC) - obj(Npb).S(iCYC-1) ; % Change in entropy of packed bed
    obj(Npb).Sirr(iCYC) = obj(Npb).Sirr(iCYC) - (obj(Npb).Sflux(iCYC,1) - obj(Npb).Sflux(iCYC,2)) ; % Net entropy increase of packed bed
    obj(Npb).Sirr(iCYC) = abs(obj(Npb).Sirr(iCYC)) ; % Has to be positive
end

switch mode
    case 'chg'
        obj.time(iCYC,1) = time ;
        fprintf(1,'Charging duration, h %8.3f\n',time/3600) ;
    case 'dis'
        obj.time(iCYC,2) = time ;
        fprintf(1,'Discharging duration,h  %8.3f\n',time/3600) ;
end

fprintf(1,str2) ;

% Flip all matrices ready for discharge cycle
for j = 1 : Npb
    obj(Npb).TS  = flip(obj(Npb).TS,1) ;
    obj(Npb).TF  = flip(obj(Npb).TF,1) ;
    obj(Npb).P   = flip(obj(Npb).P,1) ;
    obj(Npb).u   = flip(obj(Npb).u,1) ;
    obj(Npb).rho = flip(obj(Npb).rho,1) ;
end

% Probably also want to flip the order of the modules

iCYC = iCYC + 1 ;
end