% Generation information
Np = 50;        % Population size
Nr = 25;        % Repository size
maxgen = 25;    % Maximum number of generations

%Want to use previous data? change previousdata to 1
previousdata=1;
Allfigures=0;

% Objective names
Obj(1).type = 'RTeff'; 
Obj(2).type = 'LCOS';
%Obj(3).type = 'CAPC';

% Decision variable names
Par(1).type = 'TH_dis0';
Par(2).type = 'PRr';
Par(3).type = 'eff';
Par(4).type = 'ploss';

nVar  = numel(Par);
Nobjs = numel(Obj); %Number of objectives = 2 or 3.

for i = 1 : Nobjs
   switch Obj(i).type
       
       case 'RTeff'
           Obj(i).title = '1 - Roundtrip efficiency';
           
       case 'LCOS'
           Obj(i).title = 'LCOS [USD/kWh]';
           
       case 'CAPC'
           Obj(i).title = 'Capital Cost [USD]';
   end
end

for i = 1 : nVar
   switch Par(i).type
       
       case 'TH_dis0'
           Par(i).title = 'Compressor Inlet Temperature';
           Par(i).var_min = 723 ;
           Par(i).var_max = 923 ;
           
       case 'PRr'
           Par(i).title = 'Discharge Pressure Ratio';
           Par(i).var_min = 0.8 ;
           Par(i).var_max = 1.5 ;
           
       case 'eff'
           Par(i).title = 'Effectiveness';
           Par(i).var_min = 0.8 ;
           Par(i).var_max = 0.99 ;
           
       case 'ploss'
           Par(i).title = 'Pressure loss';
           Par(i).var_min = 0.005 ;
           Par(i).var_max = 0.025 ;
           
       case 'NEXP'
           Par(i).title = 'Number of expansions';
           Par(i).var_min = 1 ;
           Par(i).var_max = 5 ;
           
   end
end
