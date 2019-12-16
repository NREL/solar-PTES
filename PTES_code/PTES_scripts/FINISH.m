% Close files and save copy of multi_run files
fprintf(1,'\n...END...\n');
fclose('all');
if ~multi_run
    diary off
    % Save plots
    if save_figs == 1
        save_fig(1,'./Outputs/T-s',0,0,0)
        save_fig(8,'./Outputs/Losses',0,0,0)
        %save_fig(10,'./Outputs/TQ_hex',0,0,0)
        if Load.mode==3
            save_fig(2,'./Outputs/T-s_Rankine',0,0,0)
        end
        if all([optimise,any(Load.mode==[0,2])])
            save_fig(3,'./Outputs/Golden_search',0,0,0)            
        end        
    end
end

% Releasing CoolProp AbstractState
ierr = 0; buffer_size = 10;
herr= char((1:1:buffer_size));
for i0=1:1000
    calllib('coolprop','AbstractState_free',i0, ierr,herr,buffer_size)
end