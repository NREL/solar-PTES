% If the hot store or cold store don't return to their original values (to within a small margin) then
% a constraint is exceeded
if fluidC.state(end,2).T < TC_dis0 - 5 || fluidH.state(end,2).T < TH_dis0 - 5 || Cdata.lcosM < 0
    
    f1 = 1e12 ;
    f2 = 1e12 ;
    f3 = 1e12 ;
    extra = 1e12 ;
    
else
    
    f1=1-chi_PTES_para;
    f2=Cdata.lcosM;
    extra=Cdata.cap_costM;
    f3=Cdata.cap_costM;
    
end
err= zeros(1,1);

if Nobjs ==2
    fit = [f1 f2];  %For two objectives
end
if Nobjs==3
    fit=[f1 f2 f3]; %For three objectives
end


