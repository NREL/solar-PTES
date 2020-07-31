% In this script add up the power inputs and outputs to the PCM charging
% cycle for each load period
HINave  = zeros(Nload,1) ;
HOUTave = zeros(Nload,1) ;

for i = 1 : Nload
   % Take the average of points that aren't zero 
   if numel(HINsave(HINsave(:,i)~=0,i)) > 0
       HINave(i)  = mean(HINsave(HINsave(:,i)~=0,i)) ;
       HOUTave(i) = mean(HOUTsave(HOUTsave(:,i)~=0,i)) ;
   end
end