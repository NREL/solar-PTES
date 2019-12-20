%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clear;

% Enter debugging mode if an error occurs
%dbstop if error
%dbclear all

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./PTES_scripts/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\PTES_scripts\')
end

% Set properties for plots
set_graphics

% Load CoolProp library (low-level interface)
load_coolprop

tic
PACKED_BED_INPUTS

i = 1 ;
Lcyc = false ;
den = 0.0 ;
en_prev = 0.0 ;

while ~Lcyc 

    [pbH, TsC, TfC, enC, retu] = PB_RUN(pbH, Nhot, 'chg');
    %[pbH, TsD, TfD, enD, retu] = PB_RUN(pbH, Nhot, 'dis');
    fprintf(1,'COMPLETED CYCLE %5i\n\n',i) ;
    
    den     = 100.0 * abs(enC - en_prev) / en_prev ;
    en_prev = enC ;
    
    if i == pbH.Ncyc || den < 0.1 
        Lcyc = true ;
    end
    
    i = i + 1 ;

end

toc