switch Vpnt
    case 'PRch'
        PRch  = Apnt(ipnt);
    case 'TC_0'
        TC_0 = Apnt(ipnt);
    case 'TH_0'
        TH_0 = Apnt(ipnt);
    case 'eff'
        eff = Apnt(ipnt);
    case 'eta'
        eta = Apnt(ipnt);
    case 'Nc_ch'
        Nc_ch = Apnt(ipnt);
    case 'Ne_ch'
        Ne_ch = Apnt(ipnt);
    otherwise
        error('not implemented')
end

switch Vcrv
    case 'PRch'
        PRch = Acrv(icrv);
    case 'TC_0'
        TC_0 = Acrv(icrv);
    case 'TH_0'
        TH_0 = Acrv(icrv);
    case 'eff'
        eff = Acrv(icrv);
    case 'eta'
        eta = Acrv(icrv);
    case 'Nc_ch'
        Nc_ch = Acrv(icrv);
    case 'Ne_ch'
        Ne_ch = Acrv(icrv);
    otherwise
        error('not implemented')
end