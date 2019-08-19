switch Vpnt
    case 'PR'
        PR  = Apnt(ipnt);
    case 'TC_0'
        TC_0 = Apnt(ipnt);
    case 'TH_0'
        TH_0 = Apnt(ipnt);
    case 'eff'
        eff = Apnt(ipnt);
    case 'eta'
        eta = Apnt(ipnt);
    otherwise
        error('not implemented')
end

switch Vcrv
    case 'PR'
        PR = Acrv(icrv);
    case 'TC_0'
        TC_0 = Acrv(icrv);
    case 'TH_0'
        TH_0 = Acrv(icrv);
    case 'eff'
        eff = Acrv(icrv);
    case 'eta'
        eta = Acrv(icrv);
    otherwise
        error('not implemented')
end