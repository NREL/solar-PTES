switch Vpnt
    case 'PRch'
        PRch  = Apnt(ipnt);
    case 'pmax'
        pmax  = Apnt(ipnt);
    case 'PRr'
        PRr   = Apnt(ipnt);
    case 'TC_dis0'
        TC_dis0 = Apnt(ipnt);
    case 'TH_dis0'
        TH_dis0 = Apnt(ipnt);
    case 'Ran_TbotC'
        Ran_TbotC = Apnt(ipnt);
    case 'eff'
        eff = Apnt(ipnt);
    case 'ploss'
        ploss = Apnt(ipnt);
    case 'eta'
        eta = Apnt(ipnt);
    case 'Nc_ch'
        Nc_ch = Apnt(ipnt);
    case 'Ne_ch'
        Ne_ch = Apnt(ipnt);
    case 'mdot_off'
        Load0.mdot = Design_Load.mdot .* Apnt(ipnt) ;
    case 'T0_off'
        T0_off    = (T0 + Acrv(icrv)) * ones(Load.num,1) ;
    otherwise
        error('not implemented')
end

switch Vcrv
    case 'PRch'
        PRch = Acrv(icrv);
    case 'pmax'
        pmax = Acrv(icrv);
    case 'PRr'
        PRr   = Acrv(icrv);
    case 'TC_dis0'
        TC_dis0 = Acrv(icrv);
    case 'TH_dis0'
        TH_dis0 = Acrv(icrv);
    case 'Ran_TbotC'
        Ran_TbotC = Acrv(icrv);
    case 'eff'
        eff = Acrv(icrv);
    case 'ploss'
        ploss = Acrv(icrv);
    case 'eta'
        eta = Acrv(icrv);
    case 'Nc_ch'
        Nc_ch = Acrv(icrv);
    case 'Ne_ch'
        Ne_ch = Acrv(icrv);
    case 'mdot_off'
        Load0.mdot = Design_Load.mdot .* Acrv(icrv) ;
    case 'T0_off'
        T0_off    = (T0 + Acrv(icrv)) * ones(Load.num,1) ;
    otherwise
        error('not implemented')
end

if Lmulti_mdot
   Load0.mdot = multi_mdot(ipnt,icrv) * ones(Load0.num,1) ; 
end

fprintf(1,'\nMULTI RUN STEP: Point #%d of %d, Curve #%d of %d',ipnt,length(Apnt),icrv,length(Acrv))