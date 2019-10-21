  
% LOAD COOLPROP LIBRARY

% Determine Operating System
c = computer();

% Load corresponding library
switch computer
    case 'GLNXA64'
        % Addpath to linux CoolProp library
        % From old SWIG version (allows high-level interface)
        %path_to_lib = '../CoolProp/linux/coolprop_swig/build'; %specify path to coolprop shared library
        %path_to_include = '../CoolProp/linux/coolprop_swig/include'; %specify path to coolprop's include folder
        % From newer DLL version (low-level interface only)
        path_to_lib = '../CoolProp/linux/coolprop_dll'; %specify path to coolprop shared library
        path_to_include = '../CoolProp/linux/coolprop_dll'; %specify path to coolprop's include folder
        % Loading shared library
        if ~libisloaded('coolprop') %checking whether library is already loaded
            addpath(path_to_lib)
            addpath(path_to_include)
            libname = 'libCoolProp'; %OSX and linux
            if ispc
                libname = 'CoolProp';
            end
            loadlibrary(libname,'CoolPropLib.h','includepath',path_to_include,'alias','coolprop'); % loading library with alias coolprop
            disp('loaded CoolProp shared library.')
            disp('loaded these functions: ')
            libfunctions coolprop
        end        
    case 'PCWIN64'
        % Addpath to windows CoolProp folder
        %addpath('.\CoolProps\');
        %addpath('..\CoolProp\windows\');
        
        path_to_lib = '../CoolProp/windows/libs/shared'; %specify path to coolprop shared library
        path_to_include= '../CoolProp/linux/coolprop/include'; %specify path to coolprop's include folder

        % Loading shared library
        if ~libisloaded('coolprop') %checking whether library is already loaded
            addpath(path_to_lib)
            addpath(path_to_include)
            libname = 'libCoolProp'; % OSX and linux
            if ispc
                libname = 'CoolProp';
            end
            loadlibrary(libname,'CoolPropLib.h','includepath',path_to_include,'alias','coolprop'); % loading library with alias coolprop
            disp('loaded CoolProp shared library.')
            disp('loaded these functions: ')
            libfunctions coolprop
        end
end