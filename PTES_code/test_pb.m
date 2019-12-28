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

i    = 1 ;
iCYC = 1 ;
Lcyc = false ;
den  = 0.0 ;
en_prev = 0.0 ;

for ii = 1 : Nhot
    [pbH(ii), pbH(ii).H(iCYC), pbH(ii).S(iCYC)] = PB_ENERGY(pbH(ii), fluidS) ; % Evaluate energy at the start of time
end
iCYC = iCYC + 1;

while ~Lcyc 
    
    [pbH, TsC, TfC, iCYC] = PB_RUN(pbH, Nhot, fluidS, iCYC, 'chg');
    %[pbH, TsD, TfD, iCYC] = PB_RUN(pbH, Nhot, fluidS, iCYC, 'dis');
    fprintf(1,'COMPLETED CYCLE %5i\n\n',i) ;
    
    den     = 100.0 * abs(pbH(1).H(1) - en_prev) / en_prev ;
    en_prev = pbH(1).H(1) ;
    
    if i == pbH.Ncyc || den < 0.1 
        Lcyc = true ;
    end
    
    i = i + 1 ;
    iCYC = 2 ;

end

toc