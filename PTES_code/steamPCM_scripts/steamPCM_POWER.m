% In this script add up the power inputs and outputs to the PCM charging
% cycle for each load period
HINave  = zeros(Nload,1) ;
HOUTave = zeros(Nload,1) ;

for i = 1 : Nload
   % Take the average of points that aren't zero 
   HINave(i) = 0.;
   HOUTave(i) = 0 ;
   for j = 1 : Ntank
       if numel(HINsave(HINsave(:,i,j)~=0,i,j)) > 0
           HINave(i)  = HINave(i) + mean(HINsave(HINsave(:,i,j)~=0,i,j)) ;
           HOUTave(i) = HOUTave(i) + mean(HOUTsave(HOUTsave(:,i,j)~=0,i,j)) ;
       end
   end
end

for i = 2 : Nload-1
    if isnan(HINave(i)) || isinf(HINave(i))
        HINave(i)  = 0.5 * (HINave(i+1) + HINave(i-1));
        HOUTave(i)  = 0.5 * (HOUTave(i+1) + HOUTave(i-1));
        warning('NaNs and Infs found when calculating power flows. Some dodgy interpolation has been done');
    end
end




figure(1)
plot(dsg_pow/1e6,'-'); hold on
plot(tes_pow/1e6,'-'); 
plot((HINave-HOUTave)/1e6,'o')
plot((dsg_pow - (HINave-HOUTave))/1e6,'.')
xlabel('Hour')
ylabel('Power, MW$$_\textrm{th}$$')
legend('Power from solar','Required power in/out of TES','Provided power in/out of TES','Power to load')
hold off

%{
hours = linspace(0,72,72);
figure(2)
plot(hours,dsg_pow(1:72)/1e6,'-','Color',c_pale_orange,'LineWidth',1.5); hold on
plot(hours,(HINave(1:72)-HOUTave(1:72))/1e6,'-','Color',c_dark_orange,'LineWidth',1)
plot(hours,(dsg_pow(1:72) - (HINave(1:72)-HOUTave(1:72)))/1e6,':','Color',c_pale_blue,'LineWidth',2)
legend('Power from solar','Power in/out of TES','Power to load')
xlabel('Hour')
ylabel('Power, MW$$_\textrm{th}$$')
ytickformat('%,.1f')
xlim([0 72])
ylim([-2 4.5])
hold off

figure(3)
plot(hours,dsg_pow(100:171)/1e6,'-','Color',c_pale_orange,'LineWidth',1.5); hold on
plot(hours,(HINave(100:171)-HOUTave(100:171))/1e6,'-','Color',c_dark_orange,'LineWidth',1)
plot(hours,(dsg_pow(100:171) - (HINave(100:171)-HOUTave(100:171)))/1e6,':','Color',c_pale_blue,'LineWidth',2)
legend('Power from solar','Power in/out of TES','Power to load')
xlabel('Hour')
ylabel('Power, MW$$_\textrm{th}$$')
ytickformat('%,.1f')
xlim([0 72])
ylim([-2 4.5])
hold off
%}