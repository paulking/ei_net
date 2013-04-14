% RescaleArray - Rescale an array to the [0 1] range
%
% Rescale an array so that it uses the full [0, 1] range.
% The rescaling is linear.
%
% Usage:
%    x = RescaleArray( x );
%
% Inputs:
%    x   = the array to rescale
%
% Output:
%    x   = the rescaled array

%--------------------------------------------------------------------------
function x = RescaleArray( x )

    minValue = min( x(:) );
    maxValue = max( x(:) );

    if maxValue > minValue
        coeff = 1 / (maxValue - minValue);   % multiply is faster than divide
        x = coeff * (x - minValue);
    else
        % special case: all elements have the same value!
        newValue = max(min(maxValue, 1), 0);  % constrain to 0-1 range
        x(:) = newValue;
    end
end
