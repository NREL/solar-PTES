function [HX] = set_hex_geom(HX, varargin)
%SET_HEX Determine the geometry of a heat exchanger
%
%   There are two scenarios when this is required:
%
%   1) The heat exchangers are defined in 'geom' mode so geometry is
%   required. The geometry is obtained to satisfy the performance
%   objectives set by DT and ploss.
%
%   2) The heat exchanger was already solved using the 'eff' or 'DT' modes
%   but the geometry should now be estimated for economic calculations.
%
%   Note that SET_HEX has the same inputs as HEX_FUNC, which is used in the
%   first scenario.
%
%   The sizing procedure follows the steps below:
%
%   (a) To simplify things, assume that both sides have the same hydraulic
%   diameter and cross-sectional area (this leads to slightly sub-optimal
%   designs and conservative cost estimates, but seems necessary to speed
%   things up at run-time while obtaining geometries that accurately match
%   the performance objectives)
%
%   (b) For a given hydraulic diameter, iterate over different values of
%   Af.
%
%   (c) For each value of Af, obtain in the following order: G1 and G2, Re1
%   and Re2, St2 and St2, Cf1 and Cf2, and the overall heat transfer
%   coefficient at each heat exchanger section.
%
%   (d) Compute required heat transfer area and total pressures loss.

% Select inputs according to HX.model
switch HX.model
    case 'geom'
        if nargin~=8
            error('incorrect number of inputs')
        end
        iL     = varargin{1};
        fluidH = varargin{2};
        iH     = varargin{3};
        fluidC = varargin{4};
        iC     = varargin{5};
        mode   = varargin{6};
        par    = varargin{7};
        
    case {'eff','DT'}
        if nargin~=1
            error('incorrect number of inputs')
        end
        
    otherwise
        error('not implemented')
end

% This function should only be called if the geometry has not been set.
if HX.Lgeom_set
    error('The heat exchanger geometry has already been set!')
end

% Set geometry according to scenario
switch HX.model
    case 'geom'
        
        % The design process is done twice because the pressure loss of
        % only one of the channels (the channel with the largest pressure
        % loss) is known a priory.
        for i0=1:2
            % Use hex_func to obtain temperature profiles. Temporarily set
            % HX.model to 'eff' in order to achieve this.
            HX.model = 'eff';
            [HX,~] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par);
            HX.model = 'geom';
            
            % Employ the compute_pressure function to determine for which
            % value of Af the pressure loss is the same as the objective.
            % Check whether f1 changes sign over interval. If it doesn't,
            % choose Af that is closest to zero.
            f1 = @(Af) hex_compute_area(HX,iL,'Af',Af,0);
            Afmin = 1e-3;
            Afmax = 1e+3;
            f1min = f1(Afmin) ;
            f1max = f1(Afmax) ;
            if f1min*f1max >= 0
                if abs(f1min) < abs(f1max)
                    Af0 = Afmin ;
                else
                    Af0 = Afmax ;
                end
            else
                opt = optimset('TolX',Afmin/1e3,'Display','notify');
                Af0 = fzero(f1,[Afmin,Afmax],opt);
            end
            [~,HX] = f1(Af0);
            
            % Update the 'design' values of ploss (either plossH0 = ploss
            % or plossC0 = ploss, with the other one being lower) and
            % repeat design process once.
            HX.plossH0 = HX.DppH(iL);
            HX.plossC0 = HX.DppC(iL);
        end
        
    case {'eff','DT'}
        
        if strcmp(HX.name,'rej')
            return
        end
        
        iL = [];
        
        % Size the heat exchanger for the first load period during which it
        % was employed
        for i = 1:length(HX.H)
            if ~isempty(HX.H(i).T)
                iL = i;
                break
            end
        end
        if isempty(iL)
            return
        end
        
        % Employ the compute_pressure function to determine for which
        % value of Af the pressure loss is the same as the objective.
        % Check whether f1 changes sign over interval. If it doesn't,
        % choose Af that is closest to zero.
        f1 = @(Af) hex_compute_area(HX,iL,'Af',Af,0);
        Afmin = 1e-3;
        Afmax = 1e+3;
        f1min = f1(Afmin) ;
        f1max = f1(Afmax) ;
        if f1min*f1max >= 0
            if abs(f1min) < abs(f1max)
                Af0 = Afmin ;
            else
                Af0 = Afmax ;
            end
        else
            opt = optimset('TolX',Afmin/1e4,'Display','notify');
            Af0 = fzero(f1,[Afmin,Afmax],opt);
        end
        [~,HX] = f1(Af0);
        
    otherwise
        error('not implemented')
end

% Note that geometry is now set
HX.Lgeom_set = true;

end