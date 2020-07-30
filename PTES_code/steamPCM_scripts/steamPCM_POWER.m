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


figure(1)
plot(dsg_pow,'-'); hold on
plot(tes_pow,'-'); 
plot(HINave-HOUTave,'o')
plot(dsg_pow - (HINave-HOUTave),'.')
legend('Power from solar','Required power in/out of TES','Provided power in/out of TES','Power to load')
hold off