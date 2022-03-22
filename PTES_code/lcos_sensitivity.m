% Rubbish script to investigate the sensitivity of LCOS to parameters
% Have to run PTES first

Npt = 25 ;

CcapM  = Cdata.cap_costM ;
CcapV  = linspace(Cdata.cap_costM - Cdata.cap_costSD,Cdata.cap_costM + Cdata.cap_costSD, Npt) ;

PelM  = 0.5*(min(price) + max(price));
PelV  = linspace(min(price), max(price), Npt) ;

OM    = 0.5*(min(OnM) + max(OnM)) ;
OV    = linspace(min(OnM), max(OnM), Npt) ;

NcycM = 0.5*(250+365) ;
NcycV = linspace(250,365,Npt) ;

FCRM = 0.075 ;
FCRV = linspace(0.05,0.1,Npt) ;


Wout = (E_net_dis - heater_in)/(1e3*3600) ;
Win  = -E_net_chg/(1e3*3600) ;

lcos_var = zeros(Npt,5) ;


lcos_var(:,1) = (CcapV .* (FCRM + OM) + PelM.*Win.*NcycM) ./ (Wout.*NcycM) ;
lcos_var(:,2) = (CcapM .* (FCRV + OM) + PelM.*Win.*NcycM) ./ (Wout.*NcycM) ;
lcos_var(:,3) = (CcapM .* (FCRM + OV) + PelM.*Win.*NcycM) ./ (Wout.*NcycM) ;
lcos_var(:,4) = (CcapM .* (FCRM + OM) + PelV.*Win.*NcycM) ./ (Wout.*NcycM) ;
lcos_var(:,5) = (CcapM .* (FCRM + OM) + PelM.*Win.*NcycV) ./ (Wout.*NcycV) ;

x = linspace(0,1,Npt);
w = 2 ;

figure(1)
p(1) = plot(x,lcos_var(:,1)/lcos_var(round(Npt/2),1),'-','LineWidth',w) ; hold on;
p(2) = plot(x,lcos_var(:,2)/lcos_var(round(Npt/2),2),'-o','LineWidth',w) ;
p(3) = plot(x,lcos_var(:,3)/lcos_var(round(Npt/2),3),'-s','LineWidth',w) ;
p(4) = plot(x,lcos_var(:,4)/lcos_var(round(Npt/2),4),':','LineWidth',w) ;
p(5) = plot(x,lcos_var(:,5)/lcos_var(round(Npt/2),5),'--','LineWidth',w) ; hold off;

legend('Capital cost','FCR','O\&M','Electricity price','Number of cycles','Location','best')

xlabel('Normalized parameter, $(P - P_\mathrm{lo}) / (P_\mathrm{hi} - P_\mathrm{lo}$)')
ylabel('LCOS, $$\$/\mathrm{kWh_{e}}$$')