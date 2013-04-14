% logsamples - Generate a list of logarithmically-spaced values
%
% Return a list of n sample values ranging from firstVal to lastVal. The values
% are evenly spaced on a logarithmic scale and may include fractional values.
%
% Usage:
%    logsamples(firstVal, lastVal, n)
%
% Inputs:
%    firstVal   = the smallest (first) sample value
%    lastVal    = the largest (last) sample value
%    n          = number of sample value to return
%
% Outputs:
%    none
%
%
% Created:   8/16/2011, Paul King
%--------------------------------------------------------------------------
function valueList = logsamples( firstVal, lastVal, n )

    % calculate power multiplier to span range from min to max
    z = (lastVal/firstVal) ^ (1/(n-1));
    
    % compute the value list
    valueList = firstVal * z.^(0:n-1);
end
