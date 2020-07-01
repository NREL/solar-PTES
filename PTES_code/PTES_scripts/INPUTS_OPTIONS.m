% Code options
multi_run  = 0; % run cycle several times with different parameters?
optimise   = 0; % optimise cycle?
make_plots = 1; % make plots?
save_figs  = 0; % save figures at the end?
make_hex_plots = 0; % make plots of heat exchangers?


if multi_run==1 
    % Set variable along curves
    Vpnt = 'PRr';  % variable along curve
    Npnt = 10;            % points on curve
    pnt1 = 0.8;    % min value
    pnt2 = 1.6;    % max value
    Apnt = linspace(pnt1,pnt2,Npnt); % array
    
    % Set variable between curves
    Vcrv = 'eff';
    Acrv = [0.91, 0.95, 0.97];
    Ncrv = numel(Acrv);
       
    % Delete previous files
    %delete('./Outputs/Multi_run/*.csv')
    
    % Store information on the variables being changed along the multi-run
    % calls
    save('./Outputs/Multi_run/Multi_run_var.mat',...
        'Vpnt','Npnt','Apnt','Vcrv','Ncrv','Acrv');
    %csvwrite('./Outputs/Multi_run/Multi_run_var.csv',Apnt);

    
else 
    Npnt=1; Ncrv=1; Vpnt='Nil';
    % Start new logfile
    delete ./Outputs/log.txt
    diary  ./Outputs/log.txt
    
end


% Save initial values of Load structure (to be reset during the INITIALISE
% subroutine, in case any values changed during running time, e.g.
% discharge time)
