%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCM with steam
% This code models a Phase Change Material (PCM) which surrounds a pipe 
% through which evaporating or condensing steam passes. The code uses 
% routines developed by Pau Farres-Antunez (Cambridge University) and
% Josh McTigue (NREL). The method is based partly on the work of Sharan et
% al. (2018) and partly on McTigue's PhD thesis.
%
% In v2 we no longer assume that both steam and PCM are changing phase at
% the same time. They may have variable temperatures.
%
% In v3, unsteady steam terms seem to be small - try a new steam routine
%
% In v4, we make further simplifications. Assume a steady-state solution
% for each instance in time. Step forward only the PCM equation in time,
% and iterate the steam profile based on this. Also include proper
% calculation of the heat transfer coefficients.
%
% In v4.1, the PCM equation is written in terms of enthalpy which improves
% the stability of the solution. Have also added charge and discharge
% modes.
%
% In v5, the program is seperated into different files to make a better
% program structure.
%
% Author: Josh McTigue, Pau Farres-Antunez
% 13 May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% START PROGRAM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Workspace variables
clear;

% Determine Operating System
c = computer();

% Add paths
switch computer
    case 'GLNXA64' %Linux
        addpath('./Classes/','./Generic/','./Functions/','./steamPCM_scripts/','./Other/')
    case 'PCWIN64' %Windows
        addpath('.\Classes\','.\Generic\','.\Functions\','.\steamPCM_scripts\','.\Other\')
end
% Set properties for plots
set_graphics
% Load CoolProp library (low-level interface)
load_coolprop
tic

% Load input file
steamPCM_INPUTS
steamPCM_INITIALISE

% Set up the geometry of each tube
steamPCM_SETUP
mdot(2) = (1/3) * mdot(2) ;
% Run through each load cycle
while Iload <= Nload
    steamPCM_LOADCYCLE
    steamPCMv5
    Iload = Iload + 1 ; 
end

if Lreadload
    steamPCM_POWER
end

toc

