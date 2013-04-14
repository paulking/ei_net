% AnalyzeVisualResponse - Analyze 2D receptive field properties for a population
%
% Analyze the receptive fields of visual cells against static sine grating and
% center-surround stimuli.
%
% Usage:
%    results = AnalyzeVisualResponse(model, params)
%
% Inputs:
%    model    = spiking network model
%    params   = configuration parameters:
%        cellGroupName        = name of cell group to analyze (default = 'cg_output')
%        spatialFreqSteps     = number of spatial frequencies to analyze (default = 4)
%        spatialFreqRange     = range of spatial frequencies to analyze
%                               in cycles-per-patchWidth (default = (1 4])
%        orientationSteps     = number of orientations to analyze (0 - 180 degrees)
%                               (default = 10)
%        phaseSteps           = number of phases to analyze (0 - 360 degrees)
%                               (default = 10)
%        sigmaCenter          = width of center of center-surround RF; lower value is sharper
%                               (default = 1.2)
%        sigmaSurroundFactor  = width of surround relative to center (default = 1.5)
%        csSizeSteps          = number of size step increments
%        csSizeRange          = array(1,2) of minimum and maximum center-surround size
%                               (units are fraction of image patch size)
%        csPositionSteps      = number of position steps in both X and Y direction
%        debugDisplay         = show debug images: test stimuli patches and ideal results
%
% Output:
%    results = computed results
%        cellStats            = struct array(numCells,1) containing per-cell stats
%            gaborOrientationML       = 
%            gaborOrientationMean     = 
%            gaborOrientationSigma    = 
%            gaborPhaseML             = 
%            gaborPhaseMean           = 
%            gaborPhaseSigma          = 
%            gaborSpatialFreqML       = 
%            gaborSpatialFreqMean     = 
%            gaborSpatialFreqSigmaLog = 
%            csResponsePolarityML     = is cell ON-center (1), OFF-center (-1)?
%        gaborOrientationSigma    =
%        gaborPhaseSigma          =
%        gaborSpatialFreqMean     =
%        gaborSpatialFreqSigmaLog =
%        
%
% Created:   3/31/12, Paul King
%--------------------------------------------------------------------------
function results = AnalyzeVisualResponse( model, params )

    % initialize and apply default parameter values
    if nargin < 2
        params = struct();
    end
    defaultValues.cellGroupName       = 'cg_output';
    defaultValues.spatialFreqSteps    = 4;
    defaultValues.spatialFreqRange    = [1 3];
    defaultValues.orientationSteps    = 10;
    defaultValues.phaseSteps          = 10;
    defaultValues.csSizeSteps         = 4;
    defaultValues.csSizeRange         = [.1 .5];
    defaultValues.csPositionSteps     = 5;
    defaultValues.sigmaCenter         = 1.2;
    defaultValues.sigmaSurroundFactor = 1.5;
    defaultValues.amplifyCenter       = 1;
    defaultValues.debugDisplay        = false;
    params = ApplyDefaultValues(params, defaultValues);


    % ==============================  INITIALIZE  ===================================

    i_cg_input  = NetModel_FindElement(model, 'cg_input');
    i_cg_target = NetModel_FindElement(model, params.cellGroupName);
    numCells    = model.cellGroup{i_cg_target}.numCells;


    % ==============================  ANALYZE GABOR RFs  =============================

    % generate a matrix of sine gratings at various resolutions and orientations
    theta        = (((1:params.orientationSteps) - 1) / params.orientationSteps) * pi;
    phase        = (((1:params.phaseSteps) - 1) / params.phaseSteps) * 2 * pi;
    sf_logsteps  = (((1:params.spatialFreqSteps) - 1) / (params.spatialFreqSteps-1)) ...
                       * log(params.spatialFreqRange(2) / params.spatialFreqRange(1));
    spatialFreq  = params.spatialFreqRange(1) * exp(sf_logsteps);
    numPixels    = prod(model.inputDims);
    gratings     = zeros(numPixels, numel(theta), numel(phase), numel(spatialFreq));
    [Xpos, Ypos] = meshgrid(1:model.inputDims(2), 1:model.inputDims(1));
    Xnorm        = (Xpos - model.inputDims(2)/2 - .5) / model.inputDims(2);
    Ynorm        = (Ypos - model.inputDims(1)/2 - .5) / model.inputDims(1);
    for i = 1:numel(theta)
        X = Xnorm * cos(theta(i)) - Ynorm * sin(theta(i));
        Y = Ynorm * cos(theta(i)) + Xnorm * sin(theta(i));
        for j = 1:numel(phase)
            for k = 1:numel(spatialFreq)
                Z = sin(X * 2*pi * spatialFreq(k) + phase(j));
                gratings(:,i,j,k) = Z(:);
            end
        end
    end
    inputData = reshape(gratings, numPixels, []);

    % normalize to zero mean and unit variance
    inputData = bsxfun(@minus, inputData, mean(inputData));
    inputData = bsxfun(@times, inputData, 1./std(inputData));

    % DEBUG: Show stimulus patches used for response measurement
    if params.debugDisplay
        figure(1); ShowImagePatchGrid(reshape(inputData, model.inputDims(1), model.inputDims(2), []));
    end

    % run input data through model
    testModel                        = model;
    testModel.learningRate           = 0;
    testModel.numSamplesPerBatch     = size(inputData,2);
    testModel.stats.keepSpikeHistory = true;
    testModel = RunVisionNetworkSimulation({'numInputSamples',testModel.numSamplesPerBatch, ....
            'visualInput',inputData}, testModel);

    % convert results back to response matrix
    numSpikes = sum(testModel.snapshot.spikeHistory{i_cg_target}, 3);
    response  = reshape(numSpikes, numCells, numel(theta), numel(phase), numel(spatialFreq));


    % ---------------------  analyze cell response results  -----------------------

    % find mean and std-dev orientation, phase, and log spatial freq for each cell

    % find maximum likelihood values
    [maxResponses, maxIndices] = max(reshape(response, numCells,[]), [], 2);
    siz = size(response);
    [maxThetaIdx, maxPhaseIdx, maxSpatialFreqIdx] = ind2sub(siz(2:end), maxIndices);
    
    response_theta       = reshape( mean(mean(response, 4), 3), numCells, []); 
    response_phase       = reshape( mean(mean(response, 4), 2), numCells, []); 
    response_spatialFreq = reshape( mean(mean(response, 3), 2), numCells, []); 
    
    % rescale all response matrices so that the total for each cell (each row) adds up to 1
    rescale = 1 ./ sum(response_theta, 2);
    rescale(~isfinite(rescale)) = 0;
    response_theta = bsxfun(@times, response_theta, rescale);

    rescale = 1 ./ sum(response_phase, 2);
    rescale(~isfinite(rescale)) = 0;
    response_phase = bsxfun(@times, response_phase, rescale);

    rescale = 1 ./ sum(response_spatialFreq, 2);
    rescale(~isfinite(rescale)) = 0;
    response_spatialFreq = bsxfun(@times, response_spatialFreq, rescale);

    % compute circular mean and circular sigma for theta
    % (periodic variable properties are determined from von Mises distribution)
    % (note: it is necessary to multiple and later divide theta by 2 to use 180 deg circle)
    % circular mean:  theta_mean  = atan( E[sin(theta)] / E[cos(theta)] )
    % circular sigma: theta_sigma = sqrt(1 - E[cos(theta - theta_mean)])
    theta2          = theta * 2;
    theta_ML        = theta(maxThetaIdx);
    theta_E_sin     = sum(bsxfun(@times, sin(theta2), response_theta), 2);
    theta_E_cos     = sum(bsxfun(@times, cos(theta2), response_theta), 2);
    theta_mean      = atan2(theta_E_sin, theta_E_cos) / 2;            % range = [-pi/2, pi/2]
    theta_E_cos_dev = sum( cos(bsxfun(@minus, theta,  theta_mean)) .* response_theta, 2);
    theta_sigma     = sqrt(1 - theta_E_cos_dev);

    % compute circular mean and circular sigma for phase
    phase_ML        = phase(maxPhaseIdx);
    phase_E_sin     = sum(bsxfun(@times, sin(phase), response_phase), 2);
    phase_E_cos     = sum(bsxfun(@times, cos(phase), response_phase), 2);
    phase_mean      = atan2(phase_E_sin, phase_E_cos);
    phase_E_cos_dev = sum( cos(bsxfun(@minus, phase,  phase_mean)) .* response_phase, 2);
    phase_sigma     = sqrt(1 - phase_E_cos_dev);

    % spatial frequency measures determined in log normal space
    % spatialFreq_mean     = exp( E[log(spatialFreq)] )
    % spatialFreq_sigmaLog = sqrt( E[log(spatialFreq/spatialFreq_mean)^2] )
    spatialFreq_ML       = spatialFreq(maxSpatialFreqIdx);
    spatialFreq_mean     = exp(sum(bsxfun(@times, log(spatialFreq), response_spatialFreq), 2));
    spatialFreq_logdev   = log( bsxfun(@times, spatialFreq, 1./spatialFreq_mean) );
    spatialFreq_sigmaLog = sqrt( sum(bsxfun(@times, spatialFreq_logdev.^2, response_spatialFreq), 2) );

    % DEBUG: show ideal patches using both mean and maximum likelihood response estimate vlaues
    if params.debugDisplay
        % generate "ideal" patches from estimated gabor parameters
        idealPatches = zeros(model.inputDims(1), model.inputDims(2), numCells);
        for i = 1:numCells
            X = Xnorm * cos(theta_mean(i)) - Ynorm * sin(theta_mean(i));
            Y = Ynorm * cos(theta_mean(i)) + Xnorm * sin(theta_mean(i));
            Z = sin(X * 2*pi * spatialFreq_mean(i) + phase_mean(i));
            idealPatches(:,:,i) = Z;
        end
        ShowFigure(3, 'ideal gabors averaged'); ShowImagePatchGrid(idealPatches);

        % generate "ideal" ML patches from estimated gabor parameters
        idealPatches = zeros(model.inputDims(1), model.inputDims(2), numCells);
        for i = 1:numCells
            X = Xnorm * cos(theta_ML(i)) - Ynorm * sin(theta_ML(i));
            Y = Ynorm * cos(theta_ML(i)) + Xnorm * sin(theta_ML(i));
            Z = sin(X * 2*pi * spatialFreq_ML(i) + phase_ML(i));
            idealPatches(:,:,i) = Z;
        end
        ShowFigure(4, 'ideal gabors maximum likelihood'); ShowImagePatchGrid(idealPatches);
    end


    % =======================  ANALYZE CENTER-SURROUND RFs  ==========================

    % generate an array of shifted center-surround filters
    size_logsteps  = (((1:params.csSizeSteps) - 1) / (params.csSizeSteps-1)) ...
                       * log(params.csSizeRange(2) / params.csSizeRange(1));
    sizeVals       = params.csSizeRange(1) * exp(size_logsteps);
    xCenter        = ((1:params.csPositionSteps) - .5) / params.csPositionSteps;
    yCenter        = ((1:params.csPositionSteps) - .5) / params.csPositionSteps;
    numPixels      = prod(model.inputDims);
    csRF           = zeros(numPixels, numel(yCenter), numel(xCenter), numel(sizeVals), 2);
    [Xpos, Ypos]   = meshgrid(1:model.inputDims(2), 1:model.inputDims(1));
    X              = (Xpos - .5) / model.inputDims(2);
    Y              = (Ypos - .5) / model.inputDims(1);
    for i = 1:numel(sizeVals)
        sigmaCenter   = params.sigmaCenter * sizeVals(i);
        sigmaSurround = params.sigmaSurroundFactor * sigmaCenter;
        for j = 1:numel(yCenter)
            for k = 1:numel(xCenter)
                zCenter   = normpdf(X, xCenter(k), sigmaCenter)   .* normpdf(Y, yCenter(j), sigmaCenter);
                zSurround = normpdf(X, xCenter(k), sigmaSurround) .* normpdf(Y, yCenter(j), sigmaSurround);
                Z         = zCenter * params.amplifyCenter - zSurround;
                csRF(:,j,k,i,1) = Z(:);
                csRF(:,j,k,i,2) = -Z(:);
            end
        end
    end
    inputData = reshape(csRF, numPixels, []);

    % normalize to zero mean and unit variance
    inputData = bsxfun(@minus, inputData, mean(inputData));
    inputData = bsxfun(@times, inputData, 1./std(inputData));

    % figure(2); ShowImagePatchGrid(reshape(inputData, model.inputDims(1), model.inputDims(2), []));

    % run input data through model
    testModel                        = model;
    testModel.learningRate           = 0;
    testModel.numSamplesPerBatch     = size(inputData,2);
    testModel.stats.keepSpikeHistory = true;
    testModel = RunVisionNetworkSimulation({'numInputSamples',testModel.numSamplesPerBatch, ....
            'visualInput',inputData}, testModel);

    % convert results back to response matrix
    numSpikes = sum(testModel.snapshot.spikeHistory{i_cg_target}, 3);
    response  = reshape(numSpikes, numCells, numel(sizeVals), numel(yCenter), numel(xCenter), 2);


    % ---------------------  analyze cell response results  -----------------------

    % find mean and std-dev orientation, phase, and log spatial freq for each cell

    % find maximum likelihood values
    [maxResponses, maxIndices] = max(reshape(response, numCells,[]), [], 2);
    siz = size(response);
    [maxSizeIdx, maxYIdx, maxXIdx, maxPolarityIdx] = ind2sub(siz(2:end), maxIndices);
    
    response_size     = reshape( mean(mean(mean(response, 5), 4), 3), numCells, []); 
    response_y        = reshape( mean(mean(mean(response, 5), 4), 2), numCells, []); 
    response_x        = reshape( mean(mean(mean(response, 5), 3), 2), numCells, []); 
    response_polarity = reshape( mean(mean(mean(response, 4), 3), 2), numCells, []); 
    
    % rescale all response matrices so that the total for each cell (each row) adds up to 1
    rescale = 1 ./ sum(response_size, 2);
    rescale(~isfinite(rescale)) = 0;
    response_size = bsxfun(@times, response_size, rescale);

    rescale = 1 ./ sum(response_y, 2);
    rescale(~isfinite(rescale)) = 0;
    response_y = bsxfun(@times, response_y, rescale);

    rescale = 1 ./ sum(response_x, 2);
    rescale(~isfinite(rescale)) = 0;
    response_x = bsxfun(@times, response_x, rescale);

    rescale = 1 ./ sum(response_polarity, 2);
    rescale(~isfinite(rescale)) = 0;
    response_polarity = bsxfun(@times, response_polarity, rescale);

    % size measures determined in log normal space
    % size_mean     = exp( E[log(sizeVal)] )
    % size_sigmaLog = sqrt( E[log(sizeVal/size_mean)^2] )
    size_ML       = sizeVals(maxSizeIdx);
    size_mean     = exp(sum(bsxfun(@times, log(sizeVals), response_size), 2));
    size_logdev   = log( bsxfun(@times, sizeVals, 1./size_mean) );
    size_sigmaLog = sqrt( sum(bsxfun(@times, size_logdev.^2, response_size), 2) );

    % polarity is the cell's bias towards ON-center (1) vs. OFF-center (-1)
    % polarity_mean  = E[polarity]
    % polarity_invar = E[polarity_ON] * E[polarity_OFF] / max
    polarity_ML    = 2 - 2*maxPolarityIdx;
    polarity_mean  = sum(bsxfun(@times, [1 -1], response_polarity), 2);
    polarity_invar = response_polarity(:,1) .* response_polarity(:,2) / .5^2;

    % position measures determined in normal space
    % TODO use sigma log for position variability?
    % pos_mean     = exp( E[pos] )
    % pos_sigma    = sqrt( E[(pos - posMean)^2] )
    y_ML       = yCenter(maxYIdx);
    y_mean     = sum(bsxfun(@times, yCenter, response_y), 2);
    y_dev      = bsxfun(@minus, yCenter, y_mean);
    y_sigma    = sqrt( sum(bsxfun(@times, y_dev.^2, response_y), 2) );

    x_ML       = yCenter(maxXIdx);
    x_mean     = sum(bsxfun(@times, yCenter, response_x), 2);
    x_dev      = bsxfun(@minus, yCenter, x_mean);
    x_sigma    = sqrt( sum(bsxfun(@times, x_dev.^2, response_x), 2) );

    % combine x and y into pos
    pos_ML    = [y_ML,    x_ML];
    pos_mean  = [y_mean,  x_mean];
    pos_sigma = sqrt(y_sigma.^2 + x_sigma.^2);


    % ============================  RETURN RESULTS  ==============================

    % return results for cells individually and population in aggregate
    results = struct();
    for i = 1:numCells
        results.cellStats(i).gaborOrientationML       = theta_ML(i);
        results.cellStats(i).gaborOrientationMean     = theta_mean(i);
        results.cellStats(i).gaborOrientationSigma    = theta_sigma(i);
        results.cellStats(i).gaborPhaseML             = phase_ML(i);
        results.cellStats(i).gaborPhaseMean           = phase_mean(i);
        results.cellStats(i).gaborPhaseSigma          = phase_sigma(i);
        results.cellStats(i).gaborSpatialFreqML       = spatialFreq_ML(i);
        results.cellStats(i).gaborSpatialFreqMean     = spatialFreq_mean(i);
        results.cellStats(i).gaborSpatialFreqSigmaLog = spatialFreq_sigmaLog(i);
        results.cellStats(i).csSizeML                 = size_ML(i);
        results.cellStats(i).csSizeMean               = size_mean(i);
        results.cellStats(i).csSizeSigmaLog           = size_sigmaLog(i);
        results.cellStats(i).csPolarityML             = polarity_ML(i);
        results.cellStats(i).csPolarityMean           = polarity_mean(i);
        results.cellStats(i).csPolarityInvariance     = polarity_invar(i);
        results.cellStats(i).csXML                    = x_ML(i);
        results.cellStats(i).csXMean                  = x_mean(i);
        results.cellStats(i).csXSigma                 = x_sigma(i);
        results.cellStats(i).csYML                    = y_ML(i);
        results.cellStats(i).csYMean                  = y_mean(i);
        results.cellStats(i).csYSigma                 = y_sigma(i);
        results.cellStats(i).csPosSigma               = pos_sigma(i);
    end
    results.gaborOrientationSigma    = mean(theta_sigma);
    results.gaborPhaseSigma          = mean(phase_sigma);
    results.gaborSpatialFreqMean     = mean(spatialFreq_mean);
    results.gaborSpatialFreqSigmaLog = mean(spatialFreq_sigmaLog);
    results.csSizeMean               = mean(size_mean);
    results.csSizeSigmaLog           = mean(size_sigmaLog);
    results.csPolarityMean           = mean(polarity_mean);
    results.csPolarityInvar          = mean(polarity_invar);
    results.csPosSigma               = mean(pos_sigma);
end

