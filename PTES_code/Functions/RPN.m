function [varargout] = RPN(input_pair,input1,input2,out,fluid)
%RPN Reads several thermophysical properties in a single time.
%
%   RPN supports three different modes, 'CP' (obtain properties from
%   CoolProp), 'TAB' (tabular interpolation) and 'IDL' (ideal gas).

% Set access to global structures TABS and LIB
global TABS
global LIB

% Load or create LIB file with tabular interpolation data (only first time)
if isempty(LIB)
    fileName = './LIB/LIB/LIB.mat';
    % Create Library if file does not exist or if the file is old
    renewalAge = 7; %days
    if isfile(fileName)
        % File exists. Check how old it is.
        fileInfo = dir(fileName);
        fileAge  = datenum(datetime) - fileInfo.datenum; %file age in days
        if fileAge > renewalAge
            create_LIB = 1;
        else
            create_LIB = 0;
        end
    else
        % File does not exist, needs to be created.
        create_LIB = 1;
    end
    
    % Create LIB file if necessary
    if create_LIB
        warning(['Creating LIB file with tabular interpolation data. ',...
            'This could take a couple of minutes but is only done once.']);
        tic
        LIB.Water.HmassP_INPUTS = LIB_THERM('Water','HmassP_INPUTS',100);
        toc
        tic
        LIB.Water.PQ_INPUTS     = LIB_THERM('Water','PQ_INPUTS',1e3);
        toc
        LIB.Water.name = 'Water';
        save('./LIB/LIB/LIB.mat','LIB')
    end
    
    S   = load(fileName,'LIB');
    LIB = S.LIB;
end

% Make sure that the 'out' array of query properties is a cell array and
% store its length into variable N.
if ~iscell(out)
    out = {out};
end
N = length(out);

%Preallocate the array of output arguments.
varargout = cell([1,N]);

switch fluid.read
    case 'CP' % Obtain properties from CoolProp
        
        % Select Coolprop handle. Control cases where some CoolProp backends
        % (other than HEOS) do not produce accurate results. This is seen to
        % happen at low pressures and close to the saturation curve
        handle = fluid.handle;
        if strcmp(fluid.name,'Water')
            switch input_pair
                case {'HmassP_INPUTS'}
                    if any(input2<0.8*1e5)
                        handle = fluid.HEOS;
                    end
                case {'PSmass_INPUTS','PT_INPUTS'}
                    if any(input1<0.8*1e5)
                        handle = fluid.HEOS;
                    end
                case {'PQ_INPUTS'}
                    if any(input1<0.8*1e5)
                        handle = fluid.HEOS;
                    end
                case {'QT_INPUTS'}
                    handle = fluid.HEOS;
                case {0}
                otherwise
                    error('not implemented')
            end
            
            %%{
            switch input_pair
                case {'HmassP_INPUTS','PQ_INPUTS'}
                    [outputs] = RLIB(LIB.Water.(input_pair), input_pair, input1, input2, out);
                    for io=1:N
                        varargout{io} = outputs(:,:,io);
                    end
                    %keyboard
                otherwise
                    % Call CoolProp and extract output
                    for io=1:N
                        varargout{io} = CP1(input_pair,input1,input2,out{io},handle);
                    end
                    
            end
            %}
            
            %{
            % Call CoolProp and extract output
            for io=1:N
                varargout{io} = CP1(input_pair,input1,input2,out{io},handle);
            end
            %}
            
        else
            % Call CoolProp and extract output
            for io=1:N
                varargout{io} = CP1(input_pair,input1,input2,out{io},handle);
            end
        end
        
        
    case 'TAB' % Use tabular interpolation
        
        TAB = TABS.(valid_name(fluid.name));
        
        % Obtain x array and query points (xq)
        switch input_pair
            case 'HmassP_INPUTS'
                % Tabulated data is pressure independent at the moment. Only
                % input1 (enthalpy) is used.
                x  = TAB(:,2); %enthalpy
                xq = input1;
                
            case 'PT_INPUTS'
                % Tabulated data is pressure independent at the moment. Only
                % input2 (temperature) is used.
                x  = TAB(:,1); %temperature
                xq = input2;
                
            case 'PSmass_INPUTS'
                % Tabulated data is pressure independent at the moment. Only
                % input2 (entropy) is used.
                x  = fluid.TAB(:,4); %entropy
                xq = input2;
                
            otherwise
                error('not implemented')
        end
        
        for io=1:N
            % Obtain y array
            switch out{io}
                case {'T'}                  % temperature
                    y = TAB(:,1);
                case {'H','HMASS','Hmass'}  % mass specific enthalpy
                    y = TAB(:,2);
                case {'D','DMASS','Dmass'}  % mass density
                    y = TAB(:,3);
                case {'S','SMASS','Smass'}  % mass specific entropy
                    y = TAB(:,4);
                case {'C','CPMASS','Cpmass'}             % isobaric specific heat capacity
                    y = TAB(:,5);
                case {'CONDUCTIVITY','L','conductivity'} % thermal conductivity
                    y = TAB(:,6);
                case {'VISCOSITY','V','viscosity'}       % dynamic viscosity
                    y = TAB(:,7);
                case {'PRANDTL','Prandtl'}  % Prandtl number
                    y = TAB(:,8);
                case {'Q'}                  % vapour quality
                    y = TAB(:,9);
                otherwise
                    error('not implemented')
            end
            
            % Interpolate
            varargout{io} = interp1(x,y,xq);
        end
        
    case 'IDL' % Calculate properties of an ideal gas
        
        % First obtain 'base' properties (h,p,T,s)
        switch input_pair
            case 'HmassP_INPUTS'
                h   = input1;
                p   = input2;
                
                T   = h / fluid.IDL.cp + fluid.IDL.T0 ;
                s   = fluid.IDL.cp * log(T) - fluid.IDL.R * log(p) - fluid.IDL.s0 ;
                
            case 'PT_INPUTS'
                p   = input1;
                T   = input2;
                
                h   = fluid.IDL.cp  * T - fluid.IDL.h0 ;
                s   = fluid.IDL.cp * log(T) - fluid.IDL.R * log(p) - fluid.IDL.s0 ;
                
            case 'PSmass_INPUTS'
                p   = input1;
                s   = input2;
                
                T   = fluid.IDL.T0 * exp( (s + fluid.IDL.R * log(p/fluid.IDL.P0))/fluid.IDL.cp) ;
                h   = fluid.IDL.cp  * T - fluid.IDL.h0 ;
        end
        
        % Then compute derived properties
        rho = p ./ (fluid.IDL.R * T);
        k   = fluid.IDL.k * ones(size(T));
        mu  = fvisc(fluid,T);
        cp  = fluid.IDL.cp * ones(size(T));
        cv  = fluid.IDL.cv * ones(size(T));
        Pr  = fluid.IDL.cp .* mu ./ k ;
        Q   = -1.0 * ones(size(T));
        
        % Store requested properties in outputs
        for io=1:N            
            switch out{io}
                case {'T'}
                    output1 = T;
                case {'H','HMASS','Hmass'}
                    output1 = h;
                case {'D','DMASS','Dmass'}
                    output1 = rho ;
                case {'S','SMASS','Smass'}
                    output1 = s;
                case {'C','CPMASS','Cpmass'}
                    output1 = cp ;
                case {'CVMASS'}
                    output1 = cv ;
                case {'CONDUCTIVITY','L','conductivity'}
                    output1 = k;
                case {'VISCOSITY','V','viscosity'}
                    output1 = mu;
                case {'PRANDTL','Prandtl'}
                    output1 = Pr;
                case {'Q'}
                    output1 = Q;
                otherwise
                    error('not implemented')
            end
            varargout{io} = output1;
        end
end

end