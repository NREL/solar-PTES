classdef stage_class
    properties
        w    = 0 %specific work transfer (output)
        q    = 0 %specific heat transfer (input)
        Dh   = 0 %specific enthalpy change
        sirr = 0 %specific entropy generation
        type = '0' %either 'hex', 'regen', 'compexp', 'heat_reject' or '0'
    end
end