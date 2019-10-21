function [ DH, GC, GH, AxC, AxH, Axm, AC, AH] = printed_circuit_geom( M, DC, pC, pH, mC, mH)
%function [M, C, H] = printed_circuit_geom(M, C, H)

%   Use base parameters to compute derived geometric parameters

%DC = C.D;
%pC = C.pin;
%pH = H.pin;
%mC = C.mdot;
%mH = H.mdot;
Vm = M.V;
L  = M.L;
sigma   = M.sigma;
t_D_min = M.t_D_min;
shape   = M.shape;
VR      = M.VR;

if strcmp(shape,'squared')
    % For squared pipes
    tA_DC = max([2*pC/sigma,t_D_min]);
    tB_DC = max([(pC + sqrt(4*pC^2+8*pC*sigma))/(4*sigma),t_D_min]);
    pxC   = tA_DC + 2*tB_DC;
    tA_DH = max([2*pH/sigma,t_D_min]);
    tB_DH = max([(pH + sqrt(4*pH^2+8*pH*sigma))/(4*sigma),t_D_min]);
    pxH   = tA_DH + 2*tB_DH;
    Kx  = (1 + tA_DH)*pxC/((1 + tA_DC)*pxH);
elseif strcmp(shape,'circular')
    % For circular pipes inside flat plates
    t_DC = max([pC/(2*sigma),t_D_min]);
    pxC  = 4/pi*(1 + t_DC)^2 - 1;
    t_DH = max([pH/(2*sigma),t_D_min]);
    pxH  = 4/pi*(1 + t_DH)^2 - 1;    
    Kx  = (1 + t_DH)*pxC/((1 + t_DC)*pxH);
else
    error('not implemented')
end

% Compute VmC, VmH and DH
VmC = Vm/(1 + VR);
VmH = Vm/(1 + 1/VR);
DH = Kx*VR*DC;

% Compute mass flux, cross-sectional areas and heat transfer areas
GC  = L*mC*pxC/VmC;
GH  = L*mH*pxH/VmH;
AxC = mC/GC;
AxH = mH/GH;
Axm = Vm/L;
AC = 4*L*AxC/DC;
AH = 4*L*AxH/DH;
% AxT = Ax1 + Ax2 + Axm;

% Export computed variables for output
%M.Ax = Axm;

%C.G  = GC;
%C.Ax = AxC;
%C.A  = AC;

%H.D  = DH;
%H.G  = GH;
%H.Ax = AxH;
%H.A  = AH;
end

