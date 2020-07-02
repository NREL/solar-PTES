% SIZING
tN   = tN * 3600 ;
E    = Pmax * tN ;
Vpcm = E / (Lpcm * rhopl) ;
Lp   = Vpcm ^ (1/3) ; % One way to estimate the pipe length for an approximately cubic TES
Lp   = 10. * Lp ; % Multiply by a factor to make more square
Lp   = 1.;
if PsatD > PsatC
    warning('Discharging steam pressure is larger than charging pressure. May cause problems');
end

steam = fluid_class('Water','WF','CP','HEOS',2,30);
TsatC = RP1('PQ_INPUTS',PsatC,1.0,'T',steam) ; % Saturation temperature
TsC   = TsatC ;
PsC   = PsatC ;

TsatD = RP1('PQ_INPUTS',PsatD,0.0,'T',steam) ; % Saturation temperature
TsD   = TsatD ;
PsD   = PsatD ;

Tpm   = TsC - 10 ; % Melting point of PCM
TpC   = TsC - 5 ; % Charged temperature of PCM
TpD   = TsC - 12 ; % Discharged temperature of PCM

cpD   = cps;
rhopD = rhops ;

hps   = cps * Tpm ;
hpl   = hps + Lpcm ;
if TpD < Tpm
    hpD   = cpD * TpD ;
else
    error('Not implemented')
end


% Read load data and find charge/discharge patterns
if Lreadload
    dat      = readmatrix(fload,'Range',[2 2]) ;   % Read load file
    Nload    = length(dat(:,1)) ;                  % Number of load periods
    dsg_pow  = dat(:,1) * 1e6 ; % Power coming from DSG plant, W
    dsg_q    = dat(:,2) ;       % Quality coming from DSG plant
    dsg_T    = dat(:,3) ;       % Temperature coming from DSG plant
    dsg_mdot = dat(:,4) ;       % Mass flow rate coming from DSG plant
    load_dur = 3600 ;           % Duration of each load cycle in seconds.
    
    Load     = strings(Nload,1) ;
    tes_pow  = dsg_pow - Pmax ;   % Power to/from TES. Positive for power in, negative for power out
        
    % Run through each load cycle and decide what load type it is
    for i = 1 : Nload
        if dsg_pow(i) < 0
            dsg_pow(i) = 0 ;
        end
        if dsg_pow(i) > Pmax
            Load(i) = "c" ;
        elseif dsg_pow(i) < Pmax
            Load(i) = "d" ;
        else
            Load(i) = "s" ;
        end
    end
    
    % What is the mass flow rate corresponding to Pmax
    diff      = abs(tes_pow) ; % When tes_pow = 0, this is when P = Pmax
    tes_pow_temp = tes_pow ;
    nbest     = 3 ;
    mdot_pmax = zeros(nbest,1) ;
    pow_pmax  = zeros(nbest,1) ;
    q_pmax    = zeros(nbest,1) ;
    mdot_tes  = zeros(nbest,1) ;
    pow_tes   = zeros(nbest,1) ;
    q_tes     = zeros(nbest,1) ;
    
    for i = 1 : nbest
        % Find the nbest lowest values
        [~,ind]      = min(diff) ;
        mdot_pmax(i) = dsg_mdot(ind) ;
        pow_pmax(i)  = dsg_pow(ind) ;
        q_pmax(i)    = dsg_q(ind) ;
        diff(ind)    = 1e11 ;
        
        % Also consider max power inputs to TES to find design power, mdot
        [~,ind]      = max(tes_pow_temp) ;
        mdot_tes(i) = dsg_mdot(ind) ;
        pow_tes(i)  = tes_pow(ind) ;
        q_tes(i)    = dsg_q(ind) ;
        tes_pow_temp(ind) = -1e11 ;
    end
    % Design values - these are the targets for PCM discharge
    dsg_mdot0 = mean(mdot_pmax) ;
    dsg_q0    = mean(q_pmax) ;
    dsg_P0    = Pmax ;
    
    % Design charging inputs
    tes_mdot0 = mean(mdot_tes) - dsg_mdot0 ;
    tes_pow0  = mean(pow_tes) ;
    tes_q0    = mean(q_tes) ;
    
    xsC       = tes_q0 ;
    
else
    load_dur = tN ; % The duration of each load cycle
end

Tload = load_dur * ones(Nload,1) ; % Array that contains the duration of each load
Iload       =       1 ;             % Load number counter


%**** SET UP SOME PARAMETERS ****%

% Saturation properties of steam 
hvC = RP1('PQ_INPUTS',PsatC,1.0,'H',steam) ; % Enthalpy of vapour
hlC = RP1('PQ_INPUTS',PsatC,0.0,'H',steam) ; % Enthalpy of liquid
vvC = 1./RP1('PQ_INPUTS',PsatC,1.0,'D',steam) ; % Specific volume of vapour
vlC = 1./RP1('PQ_INPUTS',PsatC,0.0,'D',steam) ; % Specific volume of liquid

hvD = RP1('PQ_INPUTS',PsatD,1.0,'H',steam) ; % Enthalpy of vapour
hlD = RP1('PQ_INPUTS',PsatD,0.0,'H',steam) ; % Enthalpy of liquid
vvD = 1./RP1('PQ_INPUTS',PsatD,1.0,'D',steam) ; % Specific volume of vapour
vlD = 1./RP1('PQ_INPUTS',PsatD,0.0,'D',steam) ; % Specific volume of liquid

csl   = RP1('PQ_INPUTS',PsC,0.0,'CPMASS',steam) ; 
rhosl = RP1('PQ_INPUTS',PsC,0.0,'D',steam) ; 
musl  = RP1('PQ_INPUTS',PsC,0.0,'VISCOSITY',steam) ; 
ktl   = RP1('PQ_INPUTS',PsC,0.0,'L',steam) ;  % Thermal conductivity of liquid

csv   = RP1('PQ_INPUTS',PsC,1.0,'CPMASS',steam) ; 
rhosv = RP1('PQ_INPUTS',PsC,1.0,'D',steam) ; 
musv  = RP1('PQ_INPUTS',PsC,1.0,'VISCOSITY',steam) ; 
ktv   = RP1('PQ_INPUTS',PsC,1.0,'L',steam) ;  % Thermal conductivity of liquid

Psc   = RP1(0,0,0,'PCRIT',steam);      % Critical pressure of steam


% Going to do something a bit crappy to avoid having to use Coolprops (very slow)
Ntab  = 50 ;
htabC = linspace(hvC,2*hvC,Ntab)' ;
htabD = linspace(hvD,2*hvD,Ntab)' ;
TtabC = zeros(Ntab,1) ;
VtabC = zeros(Ntab,1) ;
TtabD = zeros(Ntab,1) ;
VtabD = zeros(Ntab,1) ;

for i = 1 : Ntab
   TtabC(i) = RP1('HmassP_INPUTS',htabC(i),PsatC,'T',steam);
   VtabC(i) = 1./RP1('HmassP_INPUTS',htabC(i),PsatC,'D',steam);
   
   TtabD(i) = RP1('HmassP_INPUTS',htabD(i),PsatD,'T',steam);
   VtabD(i) = 1./RP1('HmassP_INPUTS',htabD(i),PsatD,'D',steam);
end

% Correlation between H and T, and H and V (specific volume)
HH = [ones(Ntab,1) htabC] ;
TcoefC = HH \ TtabC ;
VcoefC = HH \ VtabC ;

HH = [ones(Ntab,1) htabD] ;
TcoefD = HH \ TtabD ;
VcoefD = HH \ VtabD ;

% Design charge and discharge values of steam
hsC   = RP1('PQ_INPUTS',PsC,xsC,'H',steam) ;
ssC   = RP1('PQ_INPUTS',PsC,xsC,'S',steam) ; 
rhosC = RP1('PQ_INPUTS',PsC,xsC,'D',steam) ; 
musC  = RP1('PQ_INPUTS',PsC,xsC,'VISCOSITY',steam) ; 

hsD   = RP1('PQ_INPUTS',PsD,xsD,'H',steam) ;
ssD   = RP1('PQ_INPUTS',PsD,xsD,'S',steam) ; 
rhosD = RP1('PQ_INPUTS',PsD,xsD,'D',steam) ; 
musD  = RP1('PQ_INPUTS',PsD,xsD,'VISCOSITY',steam) ; 

