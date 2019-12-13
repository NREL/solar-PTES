clear

PACKED_BED_INPUTS

i = 1 ;
Lcyc = false ;
den = 0.0 ;
en_prev = 0.0 ;

while ~Lcyc 

    [pbH, TsC, TfC, enC] = PB_CHARGE(pbH, Nhot);
    [pbH, TsD, TfD, enD] = PB_DISCHARGE(pbH, Nhot);
    fprintf(1,'COMPLETED CYCLE %5i\n\n',i) ;
    
    den     = 100.0 * abs(enC - en_prev) / en_prev ;
    en_prev = enC ;
    
    if i == pbH.Ncyc || den < 0.01 
        Lcyc = true ;
    end
    
    i = i + 1 ;

end
