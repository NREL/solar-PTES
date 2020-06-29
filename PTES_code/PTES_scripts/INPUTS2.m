% Code options
multi_run  = 0; % run cycle several times with different parameters?
opti_run   = 0; % optimise cycle?
make_plots = 1; % make plots?
save_figs  = 0; % save figures at the end?
make_hex_plots = 0; % make plots of heat exchangers?


if multi_run==1 
    % Set variable along curves
    Vpnt = 'T0';  % variable along curve
    Npnt = 5;            % points on curve
    pnt1 = 15 + 273.15;    % min value
    pnt2 = 30 + 273.15;    % max value
    Apnt = linspace(pnt1,pnt2,Npnt); % array
    
    % Set variable between curves
    Vcrv = 'eta';
    Acrv = [0.91];
    Ncrv = numel(Acrv);
       
    % Delete previous files
    %delete('./Outputs/Multi_run/*.csv')
    
    % Store information on the variables being changed along the multi-run
    % calls
    save('./Outputs/Multi_run/Multi_run_var.mat',...
        'Vpnt','Npnt','Apnt','Vcrv','Ncrv','Acrv');
    %csvwrite('./Outputs/Multi_run/Multi_run_var.csv',Apnt);
 end
    
if multi_run ~=1 
    Npnt=1; Ncrv=1; Vpnt='Nil';
    % Start new logfile
    delete ./Outputs/log.txt
    diary  ./Outputs/log.txt
    
end

if opti_run ==1 
    Npnt=1; Ncrv=1; Vpnt='Nil'
end



% Save initial values of Load structure (to be reset during the INITIALISE
% subroutine, in case any values changed during running time, e.g.
% discharge time)
