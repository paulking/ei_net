% GenerateCSFilterArray - Generate an array of center-surround filters
%
% The filter generated is normalized to the unit variance downscaled
% by the number of pixels.
%
% Usage:
%    filterArray            = GenerateCSFilterArray(patchDims, numFilters, params)
%    [filterArray,location] = GenerateCSFilterArray(patchDims, numFilters, params)
%
%
% Inputs:
%    patchDims    = array(1,2) indicating row & column dimensions of input image patch
%    numFilters   = how many filters to generate
%    params       = configuration parameters
%        sigmaCenter         = filter radius in pixels; lower value is sharper
%                              (default = 1.2)
%        sigmaSurroundFactor = size of surround relative to center (default = 1.5)
%        amplifyCenter       = how much to amplify the center relative to the surround
%                              (default = 1)
%
% Outputs:
%    filterArray  = array(numPixels,numFilters) that is the array of filters
%    location     = array(numFilters,2) location of center of each filter (optional)
%
% Created:   3/31/2012, Paul King
%--------------------------------------------------------------------------
function [filterArray, location] = GenerateCSFilterArray( patchDims, numFilters, params )

    % initialize and apply default parameter values
    if nargin < 3
        params = struct();
    end
    defaultValues.sigmaCenter         = 1;
    defaultValues.sigmaSurroundFactor = 1.5;
    defaultValues.amplifyCenter       = 1;
    params = ApplyDefaultValues(params, defaultValues);


    % assign geometrically tiled center locations to each pixel
    % (this will be a grid if numFilters is square, and staggered otherwise)
    numCols  = round(sqrt(numFilters));
    idx      = (1:numFilters)';
    col      = floor((idx - .5) * numCols / numFilters);   % row # (zero-based)
    X        = (col + .5) / numCols;
    Y        = (idx - .5) * numCols / numFilters - col;
    Y        = Y + (1 - max(Y) - min(Y)) / 2;              % adjust to equalize margins

    % output: filterArray(numPixels,numFilters) - RF filter set
    sigmaCenter    = params.sigmaCenter;
    sigmaSurround  = params.sigmaSurroundFactor * sigmaCenter;
    numPixels      = prod(patchDims);
    filterArray    = zeros(numPixels, numFilters);
    [X_in, Y_in]   = meshgrid(1:patchDims(2), 1:patchDims(1));
    for i = 1:numFilters
        xCenter        = X(i) * patchDims(2) + .5;
        yCenter        = Y(i) * patchDims(1) + .5;
        filterCenter   = normpdf(X_in, xCenter, sigmaCenter)   .* normpdf(Y_in, yCenter, sigmaCenter);
        filterSurround = normpdf(X_in, xCenter, sigmaSurround) .* normpdf(Y_in, yCenter, sigmaSurround);
        filter         = filterCenter * params.amplifyCenter - filterSurround;
        filter         = filter - mean(filter(:));
        filter         = filter / std(filter(:)) / numPixels;
        filterArray(:,i) = filter(:);
    end

    % return optional location argument
    if nargout >= 2
        location = [Y, X];
    end
end
