%objectives = csvread('./Data/objectives');
if fluidH.state(2,3).T <TH_dis0
   f1=0.8+(0.9-0.8)*rand(1);
   f2=0.8+(0.9-0.8)*rand(1);
   extra=1000000000;
   f3=1000000000;
   err= zeros(1,1);
    
else
    f1=1-chi_PTES_para;
    f2=Cdata.lcosM;
    extra=Cdata.cap_costM;
    f3=Cdata.cap_costM;
    err= zeros(1,1);        
end
if Nobjs ==2
    fit = [f1 f2];  %For two objectives
end
if Nobjs==3
    fit=[f1 f2 f3]; %For three objectives
end


%To change the number of objectives
%For NSGA II, go to Optimize_NSGA2 and change the
%Value of M to 2.

