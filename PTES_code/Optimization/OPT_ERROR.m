f1= 1e12;%0.8+(0.9-0.8)*rand(1);
f2= 1e12;%0.8+(0.9-0.8)*rand(1);
extra=1e12;
f3=1e12;
err = 1;
    
if Nobjs ==2
    fit = [f1 f2];  %For two objectives
end
if Nobjs==3
    fit=[f1 f2 f3]; %For three objectives
end
   