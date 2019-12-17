% Charging function
function [obj, TS, TF, en, retu] = PB_RUN (obj, Npb, mode)

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

while ~Lend
    
    if mod(i,100) == 0
        fprintf(1,str1,i) ;
    end
    
    % Step forward one time step
    for j = 1 : Npb
        obj(Npb) = PB_TIMESTEP(obj(Npb), mode) ;
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
    
    retu(i,1) = obj(Npb).u(1,1) ;
    retu(i,2) = obj(Npb).u(end,1) ;
    i = i + 1 ;
    time = time + obj.dt ;
end

[obj, en] = PB_ENERGY(obj) ; % Evaluate energy in packed bed at end of charge

fprintf(1,str2) ;

% Flip all matrices ready for discharge cycle
for j = 1 : Npb
    obj(Npb).TS = flip(obj(Npb).TS,1) ;
    obj(Npb).TF = flip(obj(Npb).TF,1) ;
    obj(Npb).P  = flip(obj(Npb).P,1) ;
    obj(Npb).u  = flip(obj(Npb).u,1) ;
    obj(Npb).rho = flip(obj(Npb).rho,1) ;
end

% Probably also want to flip the order of the modules

end