type = 'CMP';
aim_mode = 'Paim' ;
eff_mode = 'isen' ;
eta0 = 0.9 ;
numPeriods = 10 ;

CMP1 = compexp_class(type, aim_mode, eff_mode, eta0, numPeriods) ;
gas = fluid_class('Nitrogen','WF','CP','TTSE',Load.num,30);
gas.state(1,1).T = 300;
gas.state(1,1).p = 1e5;
gas = update(gas,[1,1],1);
[CMP1,gas,i] = compexp_func (CMP1,gas,[1,1],3e5);
