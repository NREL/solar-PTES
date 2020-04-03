% This script determines whether the geometry of the heat exchanger should
% be evaluated. There are two scenarios when this is required:
%   1) The heat exchangers are defined in 'geom' mode so geometry is
%   required.
%   2) The heat exchanger is in 'eff' mode but the geometry should be
%   estimated for economic calculations

% It has the same inputs and outputs as hex_func that is subsequently called

% At the moment this is all a bit approximate, because the geometry gets
% set up on the very first call, before the cycle has converged.

function [HX, fluidH, iH, fluidC, iC] = set_hex(HX, iL, fluidH, iH, fluidC, iC, mode, par)


% If geometry is not set, then set it! Only do this once. 
if ~HX.Lgeom_set
      
    HX.D1    = 0.025 ; % Seems to be necessary to set this as a constant, otherwise the recuperators iterate to have almost zero size
    switch HX.model
        case 'eff' % If eff model, then run hex_func as normal
            [HX, fluidH, iH, fluidC, iC] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par) ;
            orig = 'eff' ;
        case 'geom' % If geom model, then change to eff model to obtain estimates for UA, NTU etc.
            HX.model = 'eff' ; % Temporary
            [HX, fluidH, ~, fluidC, ~] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par) ;
            orig = 'geom' ;
    end
    
    % Extract data 
    HX.UA0 = HX.UA(iL) ;
    HX.NTU0 = HX.NTU(iL) ;
    HX.LMTD0 = HX.LMTD(iL) ;
    
    % Set up geometry
    HX = set_hex_geom(HX, iL, fluidH, iH, fluidC, iC, mode, par, HX.NTU0, HX.ploss, HX.D1); % HXt is a temporary class
    [~,~,HX] = shell_and_tube_geom(HX.C(iL), HX.H(iL), HX) ;
    HX.Lgeom_set = true;
    
    % Set the HX mode back to its original type
    HX.model = orig ;
    
    % Run geom model with proper parameters
    switch HX.model
        case 'geom'
            [HX, fluidH, iH, fluidC, iC] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par) ;
    end
            
else
    
    % If geometry is set, then call hex_func as normal.
    [HX, fluidH, iH, fluidC, iC] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par) ;

end

% 
% 
% switch HX.model
%     case 'eff'
%     
%         % Only want to estimate the area once - the first time it's called
%         if ~HX.Lgeom_set
%             
%             [HX, fluidH, iH, fluidC, iC] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par) ;
%             
%             HX.UA0 = HX.UA(iL) ;
%             HX.NTU0 = HX.NTU(iL) ;
%             HX.LMTD0 = HX.LMTD(iL) ;
%             
%             HX = set_hex_geom(HX, iL, fluidH, iH, fluidC, iC, mode, par, HX.NTU0, HX.ploss, HX.D1); % HXt is a temporary class
%             [~,~,HX] = shell_and_tube_geom(HX.C(iL), HX.H(iL), HX) ;
%             
%             HX.Lgeom_set = true ;
%         
%         else
%             [HX, fluidH, iH, fluidC, iC] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par) ;
%         end
%         
%     case 'geom'
%         
%         % The first time around want to find the geometry for the
%         % approximate value of effectiveness
%         if ~HX.Lgeom_set
%             HX.model = 'eff' ; % Temporary
%             [HX, fluidH, ~, fluidC, ~] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par) ;
%             
%             % Extract key data
%             HX.UA0 = HX.UA(iL) ;
%             HX.NTU0 = HX.NTU(iL) ;
%             HX.LMTD0 = HX.LMTD(iL) ;
%             
%             HX = set_hex_geom(HX, iL, fluidH, iH, fluidC, iC, mode, par, HX.NTU0, HX.ploss, HX.D1); % HXt is a temporary class
%             [~,~,HX] = shell_and_tube_geom(HX.C(iL), HX.H(iL), HX) ;
%             
%             HX.model = 'geom' ; % Reset
%             HX.Lgeom_set = true ;
%         end
%         
%         % Run the geometry for this case
%         [HX, fluidH, iH, fluidC, iC] = hex_func(HX, iL, fluidH, iH, fluidC, iC, mode, par) ;
%         
%         
% end
