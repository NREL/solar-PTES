function [new_name] = valid_name(name)
%VALID_NAME Simplifies strings, such that the new string only contains
%valid characters.
%
%   Any characters that would be invalid as part of file or variable names 
%   (such as '.', ':', '[' or ']') are substituted by an underscore '_' or
%   a blank space.
%
%   USAGE:
%   [new_name] = valid_name(name)

% Copy name into new_name
new_name = name;

% Substitute any invalid characters
new_name(new_name=='.') = '';
new_name(new_name==':') = '_';
new_name(new_name=='[') = '_';
new_name(new_name==']') = '';

% If two underscores are found together, only leave one
ind  = new_name=='_';
ind2 = ind(1:end-1) & ind(2:end);
new_name(ind2) = '';

end

