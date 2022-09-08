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
        %Load0.mdot = Design_Load.mdot .* Apnt(ipnt) ;
        if ~Loffdesign
            Load0.mdot = ones(Design_Load.num,1) * Apnt(ipnt) ;
        else
            Load0.mdot = ones(Load0.num,1) * Apnt(ipnt) * Design_Load.mdot(1) ;
        end
    case 'T0_off'
        Load0.T0_off    = (T0 + Apnt(ipnt)) * ones(Load.num,1) ;
    case 'HT_A_dT'
        Load0.HT_A = ones(Load.num,1) * Apnt(ipnt) ;
    case 'HT_B_dT'
        Load0.HT_B = ones(Load.num,1) * Apnt(ipnt) ;
    case 'CT_A_dT'
        Load0.CT_A = ones(Load.num,1) * Apnt(ipnt) ;
    case 'CT_B_dT'
        Load0.CT_B = ones(Load.num,1) * Apnt(ipnt) ;
    case 'CSmode'
        Load0.CSmode = int8(Apnt(ipnt) + zeros(Load.num,1));
        if Load0.CSmode < 0 || Load0.CSmode > 2
            error('Parametric studies.\nCSmode must take a value of 0 , 1, or 2.')
        end
    case 'Wdis_req'
        Wdis_req = Apnt(ipnt);
    case 'stH'
        Design_Load.time = Apnt(ipnt) * 3600 * ones(Load.num,1);
        Load.time = Design_Load.time ;
        Load0.time = Design_Load.time ;
    case 'unbalanced'
        Design_Load.time = [stH/Apnt(ipnt);stH;stH/Apnt(ipnt);stH;stH/Apnt(ipnt);stH]*3600.;  % time spent in each load period, s
        Design_Load.mdot = [fac*Apnt(ipnt);fac;fac*Apnt(ipnt);fac;fac*Apnt(ipnt);fac];  % working fluid mass flow rate, kg/s
        Load  = Design_Load ;
        Load0 = Design_Load ;
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
        if ~Loffdesign
            Load0.mdot = ones(Design_Load.num,1) * Acrv(icrv) ;
        else
            Load0.mdot = ones(Load0.num,1) * Acrv(icrv) * Design_Load.mdot(1);
        end
    case 'T0_off'
        Load0.T0_off    = (T0 + Acrv(icrv)) * ones(Load.num,1) ;
    case 'HT_A_dT'
        Load0.HT_A = ones(Load.num,1) * Acrv(icrv) ;
    case 'HT_B_dT'
        Load0.HT_B = ones(Load.num,1) * Acrv(icrv) ;
    case 'CT_A_dT'
        Load0.CT_A = ones(Load.num,1) * Acrv(icrv) ;
    case 'CT_B_dT'
        Load0.CT_B = ones(Load.num,1) * Acrv(icrv) ;
    case 'CSmode'
        Load0.CSmode = int8(Acrv(icrv) + zeros(Load.num,1));
        if any(Load0.CSmode < 0) || any(Load0.CSmode > 2)
            error('Parametric studies.\nCSmode must take a value of 0 , 1, or 2.')
        end
    case 'Wdis_req'
        Wdis_req = Acrv(icrv);
    case 'stH'
        Design_Load.time = Acrv(ipnt) * 3600 * ones(Load.num,1);
        Load.time = Design_Load.time ;
        Load0.time = Design_Load.time ;
    otherwise
        error('not implemented')
end

if Lmulti_mdot
   Load0.mdot = multi_mdot(ipnt,icrv) * ones(Load0.num,1) ; 
end

fprintf(1,'\nMULTI RUN STEP: Point #%d of %d, Curve #%d of %d',ipnt,length(Apnt),icrv,length(Acrv))