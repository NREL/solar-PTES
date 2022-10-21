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
            Design_Load.mdot = ones(Design_Load0.num,1) * Apnt(ipnt) ;
        else
            OffD_Load.mdot = ones(OffD_Load.num,1) * Apnt(ipnt) * Design_Load0.mdot(1) ;
        end
    case 'T0_off'
        OffD_Load.T0_off    = (T0 + Apnt(ipnt)) * ones(OffD_Load.num,1) ;
    case 'HT_A_dT'
        OffD_Load.HT_A = ones(OffD_Load.num,1) * Apnt(ipnt) ;
    case 'HT_B_dT'
        OffD_Load.HT_B = ones(OffD_Load.num,1) * Apnt(ipnt) ;
    case 'CT_A_dT'
        OffD_Load.CT_A = ones(OffD_Load.num,1) * Apnt(ipnt) ;
    case 'CT_B_dT'
        OffD_Load.CT_B = ones(OffD_Load.num,1) * Apnt(ipnt) ;
    case 'CSmode'
        if any(int8(Apnt(ipnt)) < 0) || any(int8(Apnt(ipnt)) > 2)
            error('Parametric studies.\nCSmode must take a value of 0 , 1, or 2.')
        end
        Design_Load.CSmode = int8(Apnt(ipnt) + zeros(Design_Load.num,1));
        OffD_Load.CSmode = int8(Apnt(ipnt) + zeros(OffD_Load.num,1));
    case 'Wdis_req'
        Wdis_req = Apnt(ipnt);
    case 'stH'
        Design_Load.time = Apnt(ipnt) * 3600 * ones(Design_Load.num,1);
        OffD_Load.time = Apnt(ipnt) * 3600 * ones(OffD_Load.num,1);
    case 'HP_mult'
        %Design_Load.time = [stH/Apnt(ipnt);stH;stH/Apnt(ipnt);stH;stH/Apnt(ipnt);stH]*3600.;  % time spent in each load period, s
        %Design_Load.mdot = [fac*Apnt(ipnt);fac;fac*Apnt(ipnt);fac;fac*Apnt(ipnt);fac];  % working fluid mass flow rate, kg/s
        %OffD_Load  = Design_Load ;

        HP_mult = Apnt(ipnt) ;
        
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
            Design_Load.mdot = ones(Design_Load0.num,1) * Acrv(icrv) ;
        else
            OffD_Load.mdot = ones(OffD_Load.num,1) * Acrv(icrv) * Design_Load0.mdot(1) ;
        end
    case 'T0_off'
        OffD_Load.T0_off    = (T0 + Acrv(icrv)) * ones(OffD_Load.num,1) ;
    case 'HT_A_dT'
        OffD_Load.HT_A = ones(OffD_Load.num,1) * Acrv(icrv) ;
    case 'HT_B_dT'
        OffD_Load.HT_B = ones(OffD_Load.num,1) * Acrv(icrv) ;
    case 'CT_A_dT'
        OffD_Load.CT_A = ones(OffD_Load.num,1) * Acrv(icrv) ;
    case 'CT_B_dT'
        OffD_Load.CT_B = ones(OffD_Load.num,1) * Acrv(icrv) ;
    case 'CSmode'
        Design_Load.CSmode = int8(Acrv(icrv) + zeros(Design_Load.num,1));
        OffD__Load.CSmode = int8(Acrv(icrv) + zeros(OffD_Load.num,1));
        if any(int8(Acrv(icrv)) < 0) || any(int8(Acrv(icrv)) > 2)
            error('Parametric studies.\nCSmode must take a value of 0 , 1, or 2.')
        end
    case 'Wdis_req'
        Wdis_req = Acrv(icrv);
    case 'stH'
        Design_Load.time = Acrv(ipnt) * 3600 * ones(Design_Load.num,1);
        OffD_Load.time = Acrv(ipnt) * 3600 * ones(OffD_Load.num,1);
    case 'HP_mult'
        HP_mult = Acrv(ipnt) ;
    otherwise
        error('not implemented')
end

fprintf(1,'\nMULTI RUN STEP: Point #%d of %d, Curve #%d of %d',ipnt,length(Apnt),icrv,length(Acrv))