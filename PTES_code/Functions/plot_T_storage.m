function [] = plot_T_storage( tank, Load, load_types, fignum, line_specs, Celcius, varargin)
% PLOT THE STORAGE TEMPERATURES ON THE T-S DIAGRAM.
% Use the tank structure to obtain the storage temperatures at the desired
% points. The Load structure is used to obtain the index iL of the load
% period which is to be plot, as specified in load_type. Only the first
% load period that matches load_types is plotted (i.e. if there are several
% 'chg' cycles within Load.type and load_types = {chg}, only the first one
% is plot).
% Fignum is the figure number. line_specs define give the graphic
% specifications of the lines to be plot (one line for the source tank, one
% line for the sink tank).
% Celcius is a logical variable that determines whether to use degrees C
% (rather than K) in the plots.
% varargin is optinally specified and contains the arrays Namearray and
% ValueArray to specify further line properties, if desired.

% Usage example:
% plot_T_storage(HT,Load,'ran',2,{'r','r'},true,{'LineWidth','LineStyle'},{0.8,'--'});

% Open figure
figure(fignum)

% Set temperature in Celcius?
switch Celcius
    case false
        K_C = 0;
    case true
        K_C = -273.15;
end

% Obtain x axis range in plot
x = xlim();

% Screen the different Load periods
for iL=1:Load.num
    
    % Plot the first Load period that matches load_types, then break
    if any(strcmp(Load.type(iL),load_types))
        % Make plot of source tank temperature
        ref1 = plot(x,tank.A(iL).T*ones(size(x))  +K_C,line_specs{1});
        % Make plot of sink tank temperature
        ref2 = plot(x,tank.B(iL+1).T*ones(size(x))+K_C,line_specs{2});
        break
    end
    
end

if ~isempty(varargin)
    NameArray  = varargin{1};
    ValueArray = varargin{2};
    set(ref1,NameArray,ValueArray)
    set(ref2,NameArray,ValueArray)
end

end