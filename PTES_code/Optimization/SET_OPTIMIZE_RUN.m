for i = 1 : numel(x)
   switch opt_par(i).type
       
       case 'TH_dis0'
           TH_dis0 = x(i); 
           
       case 'PRr'
           PRr = x(i) ;
           
       case 'eff'
           eff = x(i) ;
           
       case 'ploss'
           ploss = x(i) ;
           
       case 'NEXP'
          Ne_ch   = round (x(i)); 
           
   end
end
