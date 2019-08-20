This version of CoolProp for Linux can be run in Matlab via the low-level interface.
The high-level interface needs to be run via the Python interface (to be installed separately).

It only requires the files CoolPropLib.h and libCoolProp.so

Add paths to folder containing the above two files, e.g.:
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
