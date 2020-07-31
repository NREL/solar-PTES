function [result] = log_mean(x)
% LOG_MEAN Returns a "logarithmic mean" of the array x.
%
% The function employs the logarithmic and the exponential functions, and
% is designed to return a weighted average such that the result stands at
% an intermediate order of magnitude (rather than at an intermediate
% absolute value).
%
% The function is only well defined for positive (real) numbers.
%
% If x is an array, log_mean(x) returns a scalar. If x is a 2D matrix,
% log_mean(x) returns a value for each column.
%
% Examples:
% log_mean([0.01, 100]) = 1.0
% log_mean([2e-6,8e4])  = 0.40
% log_mean([1,5,125,625])    = 25.0
% log_mean([100 2; 10 4; 1 8]) = [10.0 4.0]

% Check that x has no more than two dimensions
if ~ismatrix(x)
    error(strcat('The supplied array has more than 2 dimensions, but',...
        'the log_mean function is only designed for 1D and 2D matrices'));
end

if any(size(x) == 1)
    % 1D matrix
    div = length(x);
else
    % 2D matrix
    div = size(x,1);
end

result = exp(sum(log(x))/div);

end