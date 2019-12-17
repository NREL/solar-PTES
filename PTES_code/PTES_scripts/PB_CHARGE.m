% Charging function
function [obj, TS, TF, enC] = PB_CHARGE (obj, Npb)
Lend = false;
i    = 1 ;
ii   = 1 ;
time = 0 ;

if obj(Npb).Ltext
    Tx = obj(Npb).TD + obj(Npb).textC * (obj(Npb).TC - obj(Npb).TD) ;
end

while ~Lend
    
    if mod(i,100) == 0
        fprintf(1,'CHARGING: TIMESTEP # %8i\n',i) ;
    end
    
    % Step forward one time step
    for j = 1 : Npb
        obj(Npb) = PB_TIMESTEP(obj(Npb), 'chg') ;
    end
    
    % Check whether to end charging
    if obj(Npb).Ltime
        if time > obj(Npb).TMAX || time > obj(Npb).timeC
            Lend = true ;
        end
    elseif obj(Npb).Ltext
        if time > obj(Npb).TMAX || obj(Npb).TS(end,1) > Tx
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

[obj, enC] = PB_ENERGY(obj) ; % Evaluate energy in packed bed at end of charge

fprintf(1,'CHARGING COMPLETE\n\n') ;

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