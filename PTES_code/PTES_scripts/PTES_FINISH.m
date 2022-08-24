% Close files and save copy of multi_run files
fprintf(1,'\n...END...\n');
fclose('all');
if ~multi_run
    diary off
    % Save plots
    if save_figs == 1
        formats = {'epsc','fig','svg'};
        save_fig(1, './Outputs/T-s',formats)
        save_fig(8, './Outputs/Losses',formats)
        save_fig(77,'./Outputs/Costs',formats)
        if Load.mode==3
            save_fig(2,'./Outputs/T-s_Rankine',formats)
        end
        % Also save workspace
        save('./Outputs/vars.mat')
    end
end

% Releasing CoolProp AbstractState
ierr = 0; buffer_size = 10;
herr= char((1:1:buffer_size));
for i0=1:1000
    calllib('coolprop','AbstractState_free',i0, ierr,herr,buffer_size)
end