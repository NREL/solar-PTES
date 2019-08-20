% Reset to 0 all the discharge states and discharge stages
[~,ind] = size(gas.state);
for i0=1:ind
    gas.state(2,i0).p=0; %#ok<*SAGROW>
    gas.state(2,i0).T=0;
    gas.state(2,i0).mdot=0;
    gas.state(2,i0).rho=0;
    gas.state(2,i0).h=0;
    gas.state(2,i0).s=0;
end

[~,ind] = size(gas.stage);
for i0=1:ind
    gas.stage(2,i0).w=0;
    gas.stage(2,i0).q=0;
    gas.stage(2,i0).sirr=0;
    gas.stage(2,i0).Dh=0;
    gas.stage(2,i0).type='0';
end

%Reset tanks
TANKS_STORAGE;

%Reset environment rejection streams
for i0=1:length(environ.sink(2,:)),environ.sink(2,i0).DHdot = 0;end