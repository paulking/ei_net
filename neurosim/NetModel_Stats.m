% NetModel_Stats - Chart various model statistics
%
% Compute and display various model statistics in table or plot form.
% The figures to be drawn are described in the model data structure.
% Calling this function periodically will update the figures.
%
% Measure - model.stats.measure.<measureName> has the following subfields:
%    measureType     = type of measure
%        'resError'      = residual error from linear reconstruction
%        'STA'           = spike-triggered average
%        'timeSTA'       = spike-triggered average for movies
%        'correlation'   = RMS correlation of cells (or between two cell groups)
%        'spikeRate'     = moving average spike rate of a population
%        'sparseness'    = sparsity measures (population, lifetime, and activity sparseness)
%        'deltaWeight'   = compute moving average of RMS weight change to test for convergence
%        'timelineStats' = track metrics over extended time
%        'custom'        = user-defined callback function
%
% Measure - special fields for particular measureTypes:
%    'resError':
%        sourceName      = input block for reconstruction, e.g. 'cg_V1e.in_input'
%        numAvgSamples   = how many input samples to include in the moving average window
%        normalize       = normalize reconstruction magnitude to unit variance?
%                          (default = true)
%        timeSeries      = assume time-series input (default = false)
%        sigmaWeighting  = gaussian sigma for local (time-series) averaging. Only
%                          applies when timeSeries==true (default = 10)
%        rmsResErr       = (out) moving average RMS reconstruction error
%    'STA':
%        sourceName      = cell group to analyze, e.g. 'cg_V1e'
%        numAvgSamples   = number of spike samples to include in the moving average window
%                          (default = model.stats.numSTASamples)
%        STA             = (out) measured spike-triggered average (moving average)
%    'timeSTA':
%        timeInterval    = array(1,2), the time interval of the movie (-/+ # iterations).
%                          (default = [-numInterationsPerSample/2, 0])
%        numFrames       = number of movie frames to collect over timeInterval
%        STA             = (out) measured spike-triggered average movie (moving average)
%    'correlation':
%        sourceName      = a cell group ('cg_xxx') or a cell array of two
%                          cell groups (e.g. {'cg_xxx', 'cg_yyy'} indicating
%                          which cells to analyze for correlation.
%        rmsCorrelation  = (out) the measured RMS correlation
%    'spikeRate':
%        sourceName      = cell group to analyze, e.g. 'cg_V1e'
%        spikeRate       = (out) measured spike rate (moving average)
%    'sparseness':
%        sourceName      = cell group to analyze, e.g. 'cg_V1e'
%        numAvgSamples   = how many test samples to include in the moving average windows
%        timeSparseness  = (out) lifetime sparseness (moving average)
%        popSparseness   = (out) population sparseness (moving average)
%        activitySparseness = (out) activity sparseness (moving average)
%    'deltaWeight':
%        deltaInterval   = the interval in time units for calculating dW (default = 1000)
%        windowSize      = the moving average window for RMS calculation in time units
%                          (default = deltaInterval)
%        ib(i)           = information collected on input block #i (i is arbitrary)
%            name           = a human-readable name of this input block, e.g. 'V1e->V1i'
%            cgId           = cell group id of this input block
%            ibId           = cell input block id of this input block
%            srcId          = cell group id of the input source
%            dW          = the moving average RMS weight change (units = dW / timeUnit)
%        cg(i)           = information collected on cell group #i (i is arbitrary)
%            name           = a human-readable name of this cell group, e.g. 'V1e'
%            cgId           = cell group id of this input block
%            dThresh     = the moving average RMS threshold change (units = dThresh / timeUnit)
%        dW              = the moving average RMS weight change across all weight sets
%                          (units = dW / timeUnit)
%        ddW             = the derivative of dW
%        dThresh         = the moving average RMS spike threshold change across all cell groups
%                          (units = dThresh / timeUnit)
%        ddThresh        = the derivative of dThresh
%    'timelineStats':
%        metricExpr      = The metrics to track (model-specific MatLab expressions)
%        historySize  = number of historical samples to collect before recycling
%                       (default = 1000).
%    'custom':
%        measureFn       = function pointer of form m = fn(m, model, spikeHistory)
%                          of the measure to calculate.
%
% Print - model.stats.print.<printName> must have exactly one of the following subfields:
%    measureName         = the name of a measure to use (uses the main metric)
%    matlabExpr          = a string containing a MatLab expression to evaluate
%    builtin             = a string naming a built-in field to evaluate
%        'sampleCount'       = total number of samples that have been used for training
%
% Figure - model.stats.figure{} contains subfields that describe the figure
% be drawn. See NetModel_Plot for details.
%
% Usage:
%    model = NetModel_Stats(model, spikeHistory)
%
% Inputs:
%    model          = the simulation model state, which includes the full
%                     simulation state of the model, plus the following field:
%       stats           = struct with display-specific configuration parameters
%           measure            = struct of sub-structs identifying measurements
%                                to calculate. The struct field names are the measure
%                                identifiers. Each sub-struct contains fields unique
%                                to that measure.                                 
%                                the measurements to calculate. (See above for details.)
%           figure             = cell array (numFigures,1) of structs describing
%                                the figures.  (See NetModel_Plot for details.)
%           print              = struct of sub-structs identifying console
%                                print columns
%           numSTASamples      = number of spike samples to include in the moving
%                                average window for STA measure calculation.
%                                Can be overridden for each STA measure.
%           printStatsEvery    = number of displayTimeUnits between stats printing.
%                                (default = 0)
%           updateFiguresEvery = number of displayTimeUnits between figure updates.
%                                (default = 1000)
%    spikeHistory   = cell array (numCellGroups,1) of spike history matrices, each
%                     with dimensions matrix(numCells,numNetworks,numIterations),
%                     if present. This parameter and all spike histories are optional.
%
% Output:
%    model          = updated model (only stats-tracking fields may be updated)
%
% See also: NetModel_Plot
%
% Created:   10/21/11, Paul King
%--------------------------------------------------------------------------
function model = NetModel_Stats( model, spikeHistory )

    % perform stats module initialization the first time through
    if ~isfield(model.stats, 'initialized') || ~model.stats.initialized

        % initialize figures
        for i = 1:numel(model.stats.figure)
            fig = model.stats.figure{i};
            if ~isempty(fig) && (~isfield(fig, 'figureNum') || fig.figureNum > 0)
                [fig, model] = InitializeFigure(fig, model);
                model.stats.figure{i} = fig;
            end
        end

        model.stats.initialized = true;
    end


    % ======================  PREPROCESS SPIKE HISTORY  =====================

    [model, spikeHistory, numTotalIterations] = PreprocessSpikeHistory(model, spikeHistory);

    % figure out the number of complete samples represented by this call
    if numTotalIterations == (model.numIterationsPerSample * model.numSamplesPerBatch)
        numSamples = model.numSamplesPerBatch;
    elseif numTotalIterations <= model.numIterationsPerSample
        % 1 if we have a complete iteration set, otherwise 0
        numSamples = double(model.snapshot.currentIterationIdx == model.numIterationsPerSample);
    else
        error('partial simulations not yet supported');
    end

    % update running count of training samples
    if model.learningRate > 0
        model.stats.trainingSampleCount = model.stats.trainingSampleCount + numSamples;
    end

    % save spike history (for debugging, can take up a lot of memory)
    if model.stats.keepSpikeHistory && ~isempty(spikeHistory) 
        % store snapshot.spikeHistory as 'single' to save on memory
        % ('logical' would be even more efficient, but generates errors in some uses)
        sh = [];
        for i = 1:numel(spikeHistory)
            sh{i} = single(spikeHistory{i});    % single uses 1/2 memory, logical 1/8
        end
        model.snapshot.spikeHistory = sh;
    end


    % ======================  CALCULATE MEASUREMENTS  =====================

    % iterate through measurements identified by struct field name
    if ~isempty(spikeHistory)
        measureList = fieldnames(model.stats.measure);
        for i = 1:numel(measureList)
            m = model.stats.measure.(measureList{i});
            if isempty(m) || (isfield(m, 'enabled') && ~m.enabled)
                continue;
            end

            switch m.measureType

            case 'resError'
                % calculate moving average RMS residual error
                m = Measure_LinearResidualError(m, model, spikeHistory);

            case 'STA'
                % calculate spike-triggered average
                m = Measure_STA(m, model, spikeHistory);

            case 'timeSTA'
                % calculate spike-triggered average that uses time dimension
                m = Measure_TimeSTA(m, model, spikeHistory);

            case 'correlation'
                % calculate cell activity correlation
                m = Measure_Correlation(m, model, spikeHistory);

            case 'spikeRate'
                % calculate spikeRate for a cell group
                m = Measure_SpikeRate(m, model, spikeHistory);

            case 'sparseness'
                % calculate various sparseness measures for a cell group
                m = Measure_Sparseness(m, model, spikeHistory);

            case 'deltaWeight'
                % calculate the RMS change in weights over time to determine convergence
                m = Measure_DeltaWeight(m, model, spikeHistory);

            case 'timelineStats'
                % record network metrics and statistics over extended time
                m = Measure_TimelineStats(m, model, spikeHistory);

            case 'custom'
                % execute custom measurement function
                m = feval(m.measureFn, m, model, spikeHistory);

            otherwise
                error('Unknown measure type "%s"', m.measureType);
            end
            model.stats.measure.(measureList{i}) = m;
        end
    end


    % ===========================  PRINT STATISTICS  ==========================

    % print lines of statistics every so often (model.stats.printStatsEvery)
    model.stats.printStatsWait = model.stats.printStatsWait ...
            + numTotalIterations / model.displayTimeScale;
    if model.stats.printStatsEvery > 0 ...
            && model.stats.printStatsWait >= model.stats.printStatsEvery
        model.stats.printStatsWait = 0;                  % reset counter
        model = PrintStatsLine(model);
    end


    % =============================  DISPLAY FIGURES  ==========================

    % display figures every so often
    model.stats.updateFiguresWait = model.stats.updateFiguresWait ...
            + numTotalIterations / model.displayTimeScale;
    if model.stats.updateFiguresEvery > 0 ...
            && model.stats.updateFiguresWait >= model.stats.updateFiguresEvery
        model.stats.updateFiguresWait = 0;   % reset counter

        % draw or update each figure
        for i = 1:numel(model.stats.figure)
            fig = model.stats.figure{i};
            
            % determine figureNum and whether or not to display figure
            if isempty(fig)
                continue;               % disabled figure
            elseif ~isfield(fig, 'figureNum') || isempty(fig.figureNum)
                fig.figureNum = 100+i;  % assign a figure number
            elseif fig.figureNum < 1
                continue;               % disabled figure
            end

            % only show the figure windows the first time through (for speed)
            if ~isfield(fig, 'title')
                fig.title = [];
            end
            showFigures = ~model.stats.areFiguresShown || ~ishandle(fig.figureNum) ...
                    || ~strcmp(get(fig.figureNum, 'Name'), fig.title);
            if showFigures
                figure(fig.figureNum);
            elseif get(0,'CurrentFigure') ~= fig.figureNum
                set(0, 'CurrentFigure', fig.figureNum);
            end

            % draw the plot inside the figure window
            fig = NetModel_Plot(fig, model, spikeHistory);

            % add a title to the figure window
            if showFigures
                set(fig.figureNum, 'Name',fig.title, 'NumberTitle','off');
            end

            model.stats.figure{i} = fig;
        end
        if numel(model.stats.figure) > 0
            drawnow();
        end

        model.stats.areFiguresShown = true;
    end
end


% PreprocessSpikeHistory - Aggregate multiple iterations into a single spike history if necessary
%
% Inputs:
%    model              = the model
%    spikeHistory       = spikeHistory as provided to NetModel_Stats (may contain
%                         only one iteration)
%
% Outputs:
%    model              = the model (possibly updated)
%    spikeHistory       = spikeHistory, possibly updated to be either null or complete
%    numTotalIterations = total number of iterations represented by initial spikeHistory
%
function [model, spikeHistory, numTotalIterations] = PreprocessSpikeHistory( model, spikeHistory )

    % initialize model.snapshot if missing
    if ~isfield(model, 'snapshot')
        model.snapshot.inputData            = [];
        model.snapshot.spikeHistory         = [];
        model.snapshot.potentialHistory     = [];
    end

    % calculate the total number if iterations contained in spikeHistory
    numTotalIterations = 0;
    activeCellGroupIds = [];
    for i = 1:numel(spikeHistory)
        [numCells, numNetworks, numIterations] = size(spikeHistory{i});
        if numCells > 0
            numTotalIterations = max(numTotalIterations, numNetworks * numIterations);
            activeCellGroupIds(end+1) = i;
        end
    end

    % special case: spikes are delivered a few iterations (or one) at a time
    if numTotalIterations < (model.numIterationsPerSample * model.numSamplesPerBatch)
    
        % initialize spike history aggregation if necessary
        if ~isfield(model.snapshot, 'currentIterationIdx')
            model.snapshot.currentIterationIdx  = 0;
            model.snapshot.currentSampleIdx     = 1;
            model.snapshot.spikeHistory         = [];
            for k = 1:numel(model.cellGroup)
                if ~strcmp(model.cellGroup{k}.cellType, 'disabled')
                    model.snapshot.spikeHistory{k} = zeros(model.cellGroup{k}.numCells, ...
                            model.numSamplesPerBatch, model.numIterationsPerSample, 'single');
                end
            end
        end

        % increment iteration and sample counters
        if model.snapshot.currentIterationIdx >= model.numIterationsPerSample
            model.snapshot.currentIterationIdx = 0;
            model.snapshot.currentSampleIdx = model.snapshot.currentSampleIdx + 1;
            if model.snapshot.currentSampleIdx > model.numSamplesPerBatch
                model.snapshot.currentSampleIdx = 1;
            end
        end
        
        if model.snapshot.currentIterationIdx + numTotalIterations > model.numIterationsPerSample
            error('iteration sets must fall on even batch boundaries');
        end

        % copy these iteration(s) into the spike history snapshot
        for k = activeCellGroupIds
            model.snapshot.spikeHistory{k}(:, model.snapshot.currentSampleIdx, ...
                    model.snapshot.currentIterationIdx+(1:numTotalIterations)) = spikeHistory{k};
        end
        model.snapshot.currentIterationIdx = model.snapshot.currentIterationIdx + numTotalIterations;

        % convert spikeHistory into either full spike history snapshot (if complete) or null
        if model.snapshot.currentSampleIdx == model.numSamplesPerBatch ...
                && model.snapshot.currentIterationIdx == model.numIterationsPerSample
            spikeHistory = model.snapshot.spikeHistory;
        else
            spikeHistory = [];
        end
    end
end


% ======================================================================================
% ===========================   MEASUREMENT ROUTINES   =================================
% ======================================================================================


% Measure_LinearResidualError
%
% Calculate residual error from spiking activity by backpropagating
% through the input weights (e.g. assuming a linear model). Supports both
% time-series input and static input (mean spike rate across sample).
% Also supports spiking input.
%
% Reconstructed input is normalized to the magnitude of the input before comparing.
% Input samples are assumed to have unit standard deviation.
%
% Note: Requires all iterations for at least one image sample to be present,
%       and so cannot work on a single iteration, e.g. from CellGroup_Update.
%       FIX: accumulate numSpikes until last sample is collected
% Note: timeSeries mode produces slightly different results, even if the input
%       is unchanging. The reason is that timeSeries mode uses spikeRate local
%       to a moving gaussian region rather than a spikeRate global across the sample.
%
% Inputs:
%    m            = the measure descriptor (model.stats.measure.<measureName>)
%    model        = the model (read-only)
%--------------------------------------------------------------------------
function m = Measure_LinearResidualError( m, model, spikeHistory )

    % first-time validation and variable initialization
    if ~isfield(m, 'initialized') || ~m.initialized
        defaultValues.numAvgSamples   = 2000;       % default # samples to average
        defaultValues.normalize       = true;       % normalize reconstruction magnitude?
        defaultValues.timeSeries      = false;      % default is to average across iterations
        defaultValues.sigmaWeighting  = 10;         % gaussian sigma for local (time-series) averaging
        m = ApplyDefaultValues(m, defaultValues);

        [i_cg, i_ib]      = NetModel_FindElement(model, m.sourceName);
        m.src_i_cg        = i_cg;
        m.src_i_ib        = i_ib;
        m.biasedRmsResErr = 0;
        m.biasedStats_max = 0;
        m.rmsResErr       = Inf;
        m.minResErr       = Inf;
        m.initialized     = true;
    end

    % look up cellGroup and inputBlock referants
    [numCells, numNetworks, numIterations] = size(spikeHistory{m.src_i_cg});
    W         = model.cellGroup{m.src_i_cg}.inputBlock(m.src_i_ib).weight;
    spikeHist = spikeHistory{m.src_i_cg};
    if isempty(spikeHist)
        return;
    end

    % locate inputData to be reconstructed
    input_i_cg    = model.cellGroup{m.src_i_cg}.inputBlock(m.src_i_ib).srcId;
    useModelInput = strcmp(model.cellGroup{input_i_cg}.name, 'input');
    if useModelInput
        inputData = model.snapshot.inputData;
    elseif ~m.timeSeries
        % collapse spikes over time into a spikeRate average
        inputData = mean(spikeHistory{input_i_cg}, 3);
    else
        error('cellgroup input not yet supported for time-series resError');
    end

    % calculate reconstructed input (rInputData)
    if m.timeSeries
        % generate a time-series spike weighting kernel
        [X, Y] = meshgrid(1:numIterations, 1:numIterations);
        weightingKernel = exp( -(X - Y).^2 / m.sigmaWeighting^2 );
        %weightingKernel(Y > X) = 0;    % don't count future spikes (TODO could try moving average)
        %weightingKernel = bsxfun(@times, weightingKernel, 1./sum(weightingKernel,1)); % normalize

        % generate time-series mean spike rate history using weighting kernel
        spikeHist_weighted = reshape(spikeHist, numCells*numNetworks, numIterations) * weightingKernel;

        % calculate rInputData
        rInputData = W' * reshape(spikeHist_weighted, numCells, numNetworks*numIterations);
        rInputData = reshape(rInputData, [], numNetworks, numIterations);
        if model.inputPreprocess.splitInvert
            numVals = size(model.inputData, 1);
            assert( size(rInputData, 1) == 2 * numVals );
            rInputData = rInputData(1:numVals,:,:) - rInputData(numVals+1:end,:,:);
        end
    else
        numSpikes  = sum(spikeHist, 3);
        rInputData = W' * numSpikes;
        if model.inputPreprocess.splitInvert
            numVals = size(inputData,1);
            assert( useModelInput && (size(rInputData,1) == 2*numVals) );
            rInputData = rInputData(1:numVals,:) - rInputData(numVals+1:end,:);
        end
    end

    % rescale rInput to unit variance (TODO mean seems to be near zero naturally)
    % note: std(inputData,0,1) == 1 by prior normalization, so no need to normalize here
    if m.normalize
        rescale = 1 ./ std(rInputData,0,1);
        rescale(~isfinite(rescale)) = 1;
        rInputData = bsxfun(@times, rInputData, rescale);   % per input sample
        if ~useModelInput
            % if the input is spikes and we're normalizing the reconstruction
            % then we need to normalize the input also
            rescale_in = 1 ./ std(inputData,0,1);
            rescale_in(~isfinite(rescale_in)) = 1;
            inputData = bsxfun(@times, inputData, rescale_in);   % per input sample
        end
    end
    
    % accumulate the set of reconstructed sample standard deviations for analysis
    % (activated by using the 'rStdHistogram' plot)
    if isfield(m, 'stdPoints')
        r_std = std(rInputData,0,1);
        if isempty(m.stdPoints)
            m.stdPoints = zeros( numel(r_std) * model.stats.updateFiguresEvery / numNetworks, 1);
            m.stdCount  = 0;
        elseif m.stdCount + numel(r_std) >= size(m.stdPoints, 1)
            m.stdCount = 0;
        end
        m.stdPoints(m.stdCount+1:m.stdCount+numel(r_std),1) = r_std(:);
        m.stdCount = m.stdCount + numel(r_std);
    end

    % calculate RMS residual error (RMS reconstruction error)
    resError          = inputData - rInputData;
    thisRmsResErr     = sqrt(mean(resError(:).^2));

    % calculate biased and unbiased moving average residual error
    stats_eta         = 1 - exp(- numNetworks / m.numAvgSamples);
    m.biasedStats_max = (1-stats_eta) .* m.biasedStats_max + stats_eta;
    m.biasedRmsResErr = (1-stats_eta) .* m.biasedRmsResErr + stats_eta * thisRmsResErr;
    m.rmsResErr       = m.biasedRmsResErr / m.biasedStats_max;
    m.minResErr       = min(m.minResErr, m.rmsResErr);
end


% Measure_STA
%
% Calculate or update moving spike-triggered average (STA).
%
% Inputs:
%    m            = the measure descriptor (model.stats.measure.<measureName>)
%        numAvgSamples = number of samples to include in moving average
%                        (default = stats.numSTASamples)
%    model        = the model (read-only)
%--------------------------------------------------------------------------
function m = Measure_STA( m, model, spikeHistory )

    assert(~isempty(spikeHistory), 'spikeHistory is required for STA tracking');

    % first-time variable initialization
    if ~isfield(m, 'initialized') || ~m.initialized
        assert( ~isfield(m, 'referantName') || strcmp(m.referantName, 'input'), ...
                'STA measure only supports referantName = "input"');
        % m.referantCellGroupId = 0;      % could be used later for customizable referant
        if ~isfield(m, 'numAvgSamples')
            m.numAvgSamples = model.stats.numSTASamples;         % default
        end
        m.cellGroupId = NetModel_FindElement(model, m.sourceName);
        if strcmp(model.cellGroup{m.cellGroupId}.cellType, 'disabled')
            fprintf('WARNING: Can''t track STA on disabled cell group %s\n', m.sourceName);
        end
        if size(model.snapshot.inputData, 3) > 1
            error('can''t compute STA for movie input');
        end
        numCells        = model.cellGroup{m.cellGroupId}.numCells;
        numVals         = prod(model.inputDims);
        m.biasedSTA     = zeros(numVals, numCells);
        m.biasedSTA_max = zeros(1,numCells);
        m.STA           = zeros(numVals, numCells);
        m.initialized   = true;
    end
    if isempty(spikeHistory{m.cellGroupId})
        return;
    end

    % update biased STA for this data set
    numSpikes       = sum(spikeHistory{m.cellGroupId}, 3);    
    numCellSpikes   = sum(numSpikes, 2)';
    STA_eta         = 1 - exp(- numCellSpikes ./ m.numAvgSamples);
    spikeWeight     = 1 ./ numCellSpikes;
    spikeWeight(numCellSpikes == 0) = 0;                       % no spikes = no weight!
    m.biasedSTA     = bsxfun(@times, (1-STA_eta), m.biasedSTA) ...
                    + bsxfun(@times, STA_eta .* spikeWeight, (model.snapshot.inputData * numSpikes'));

    % calculate unbiased STA
    m.biasedSTA_max = (1-STA_eta) .* m.biasedSTA_max + STA_eta;
    scale           = 1 ./ m.biasedSTA_max;
    scale(m.biasedSTA_max == 0) = 0;           % biasedSTA_max == 0 until spikes occur
    m.STA           = bsxfun(@times, m.biasedSTA, scale);
end


% Measure_TimeSTA
%
% Calculate or update moving spike-triggered average (STA) for time-series input
% (spike-triggered movie).
%
% timeInterval is an array [t_minus, t_plus] indicating the time interval range
% for the constructed STA movie relative to the spike trigger. The STA movie
% will be generated over the period [t0-t_minus, t0+t_plus], where t0 is the time
% that the spike occurred.
%
% Note: The STA move could be subsamples to a smaller number of frames,
% but this would increase rather than decrease compute cost.
%
% Inputs:
%    m            = the measure descriptor (model.stats.measure.<measureName>)
%        timeInterval  = array(1,2), the time interval of the movie (-/+ # iterations).
%                        (default = [-numInterationsPerSample/2, 0])
%        numAvgSamples = number of samples to include in moving average
%                        (default = stats.numSTASamples)
%    model        = the model (read-only)
%--------------------------------------------------------------------------
function m = Measure_TimeSTA( m, model, spikeHistory )

    assert(~isempty(spikeHistory), 'spikeHistory is required for STA tracking');

    % first-time variable initialization
    if ~isfield(m, 'initialized') || ~m.initialized
        assert( ~isfield(m, 'referantName') || strcmp(m.referantName, 'input'), ...
                'STA measure only supports referantName = "input"');
        defaultValues.numAvgSamples  = model.stats.numSTASamples;
        defaultValues.timeInterval   = [-ceil(model.numIterationsPerSample/2), 0];
        defaultValues.numFrames      = 5;
        m = ApplyDefaultValues(m, defaultValues);

        m.timeOffsets = round( m.timeInterval(1) + ((0:m.numFrames-1) / (m.numFrames-1)) ...
                              * (m.timeInterval(2) - m.timeInterval(1)) );
        m.cellGroupId = NetModel_FindElement(model, m.sourceName);
        if strcmp(model.cellGroup{m.cellGroupId}.cellType, 'disabled')
            fprintf('WARNING: Can''t track STA on disabled cell group %s\n', m.sourceName);
        end
        numCells        = model.cellGroup{m.cellGroupId}.numCells;
        numVals         = size(model.snapshot.inputData, 1);
        m.biasedSTA     = zeros(numVals, numCells, m.numFrames);
        m.biasedSTA_max = zeros(1,numCells, m.numFrames);
        m.STA           = zeros(numVals, numCells, m.numFrames);
        m.initialized   = true;
    end
    if isempty(spikeHistory{m.cellGroupId})
        return;
    end

    spikeHist = spikeHistory{m.cellGroupId};
    inputData = model.snapshot.inputData;
    [numCells, numNetworks, numIterations] = size(spikeHist);
    numVals   = size(model.snapshot.inputData, 1);

    % compute spike-triggered sum for this data set
    % (buffer contains the time-series data for constructing the STA movie over m.timeInterval)
    bufferSum   = zeros(numVals, numCells, m.numFrames);
    bufferCount = zeros(numCells, m.numFrames);
    for k = 1:m.numFrames
        % generate a single STA frame from inputData on interval [i0,i1] and spikeHist [s0,s1]
        t  = m.timeOffsets(k);
        i0 = max(1+t,1);
        i1 = min(numIterations+t, numIterations);
        s0 = i0-t;
        s1 = i1-t;
        subInputData     = reshape(inputData(:,:,i0:i1), numVals,  []);
        subSpikeHist     = reshape(spikeHist(:,:,s0:s1), numCells, []);
        bufferSum(:,:,k) = subInputData * subSpikeHist';
        bufferCount(:,k) = sum(subSpikeHist, 2);
    end

    % update biased STA for this data set
    bufferCount     = reshape(bufferCount, 1, numCells, m.numFrames);
    STA_eta         = 1 - exp(- bufferCount ./ m.numAvgSamples);
    spikeWeight     = 1 ./ bufferCount;
    spikeWeight(bufferCount == 0) = 0;                       % no spikes = no weight!
    m.biasedSTA     = bsxfun(@times, (1-STA_eta), m.biasedSTA) ...
                    + bsxfun(@times, STA_eta.*spikeWeight, bufferSum);

    % calculate unbiased STA
    m.biasedSTA_max = (1-STA_eta) .* m.biasedSTA_max + STA_eta;
    scale           = 1 ./ m.biasedSTA_max;
    scale(m.biasedSTA_max == 0) = 0;           % biasedSTA_max == 0 until spikes occur
    m.STA           = bsxfun(@times, m.biasedSTA, scale);
end


% Measure_Correlation
%
% Calculate the correlation between two cell populations, or within a cell
% population.
%
% Inputs:
%    m            = the measure descriptor (model.stats.measure.<measureName>)
%    model        = the model (read-only)
%
% Outputs:
%    m.rmsCorrelation = the RMS correlation
%--------------------------------------------------------------------------
function m = Measure_Correlation( m, model, spikeHistory )

    assert(~isempty(spikeHistory), 'spikeHistory is required for correlation tracking');

    % first-time variable initialization
    if ~isfield(m, 'initialized') || ~m.initialized
        if ~isfield(m, 'numAvgSamples')
            m.numAvgSamples = 2000;           % TODO use model.stats.numSTASamples?
        end
        if ~isfield(m,'sourceName') || isempty(m.sourceName)
            error('sourceName must be provided for correlation measure');
        elseif ischar(m.sourceName)
            m.cellGroupIds(1:2) = NetModel_FindElement(model, m.sourceName);
        else
            m.cellGroupIds(1)   = NetModel_FindElement(model, m.sourceName{1});
            m.cellGroupIds(2)   = NetModel_FindElement(model, m.sourceName{2});
        end
        m.rmsCorrelation       = NaN;
        m.biasedRmsCorrelation = 0;
        m.biasedStats_max      = 0;
        m.initialized          = true;
    end

    % measure the cell correlation as a moving average
    [numCells1, numNetworks, numIterations] = size(spikeHistory{m.cellGroupIds(1)});
    stats_eta  = 1 - exp(- numNetworks/m.numAvgSamples);

    % TODO for time-series covariance, need to use gaussian time windows for numSpikes

    % calculate correlation coefficient matrix
    numSpikes1 = sum(spikeHistory{m.cellGroupIds(1)}, 3);
    if m.cellGroupIds(1) == m.cellGroupIds(2)
        % self-correlation - easy!
        cor = corrcoef(numSpikes1');
        cor(1 : numCells1+1 : end) = NaN;       % set diagonal (self-correlation) to NaN
    else
        % cross-correlation: find concatenated self-correlation and take only one quadrant
        numSpikes2 = sum(spikeHistory{m.cellGroupIds(2)}, 3);
        cor        = corrcoef( cat(2, numSpikes1', numSpikes2') );
        cor        = cor(1:numCells1, numCells1+1:end);   % select one cross-correlation quadrant
    end

    % calculate moving RMS average of correlation coefficients
    %thisRmsCorrelation = sqrt(mean(max(cor(isfinite(cor)),0).^2));   % treat neg correlation as zero
    thisRmsCorrelation = sqrt(mean(cor(isfinite(cor)).^2));
    if ~isnan(thisRmsCorrelation)
        m.biasedRmsCorrelation = (1-stats_eta) .* m.biasedRmsCorrelation ...
                + stats_eta * thisRmsCorrelation;
    end
    m.biasedStats_max = (1-stats_eta) .* m.biasedStats_max + stats_eta;
    m.rmsCorrelation  = m.biasedRmsCorrelation / m.biasedStats_max;
end


% Measure_SpikeRate - Track the average spike rate in spikes-per-timeUnit
%
% Inputs:
%    m            = the measure descriptor (model.stats.measure.<measureName>)
%        numAvgSamples  = number of samples to include in moving average
%    model        = the model (read-only)
%--------------------------------------------------------------------------
function m = Measure_SpikeRate( m, model, spikeHistory )

    assert(~isempty(spikeHistory), 'spikeHistory is required for spike rate tracking');

    % first-time variable initialization
    if ~isfield(m, 'initialized') || ~m.initialized
        if ~isfield(m,'sourceName') || isempty(m.sourceName)
            error('sourceName must be provided for spike-rate measure');
        else
            i_cg = NetModel_FindElement(model, m.sourceName);
            m.cellGroupId = i_cg;
        end
        if ~isfield(m, 'numAvgSamples')
            m.numAvgSamples = 200;
        end
        m.biasedSpikeRate = 0;
        m.biasedStats_max = 0;
        m.spikeRate       = 0;
        m.popSpikeRate    = 0;
        m.initialized     = true;
    end

    % calculate the mean spike rate for this sample set
    spikeHist = spikeHistory{m.cellGroupId};
    [numCells, numSamples, numIterations] = size(spikeHist);
    thisSpikeRate = mean(reshape(spikeHist, numCells, []), 2) / model.simTimeStep;
    if isempty(thisSpikeRate)
        return;    % skip disabled cell group
    end

    % update the moving average of the spike rate
    stats_eta         = 1 - exp(- numSamples/m.numAvgSamples);
    m.biasedSpikeRate = (1-stats_eta) .* m.biasedSpikeRate + stats_eta * thisSpikeRate;
    m.biasedStats_max = (1-stats_eta) .* m.biasedStats_max + stats_eta;
    m.spikeRate       = m.biasedSpikeRate ./ m.biasedStats_max;
    m.popSpikeRate    = mean(m.spikeRate);

    % accumulate the spikeRate per sample data for later analysis
    % This is activated by using the 'spikeRateHistogram' plot. An array containing
    % the total number of spikes generated by each training sample will be maintained
    % so that the distribution of total spike rates can be analyzed. The number stored
    % is based on the frequency of figure updating by default (?)
    if isfield(m, 'historySize') && m.historySize > 0
        if ~isfield(m, 'srPointsPop') || m.srPointsPop
            % default is to track sample activity across all cells (srPointsPop mode)
            % (count spikes-per-sample and convert to average spikes-per-timeUnit)
            srSamples = sum(sum(spikeHist,3),1);
            srSamples = srSamples / (numCells * numIterations * model.simTimeStep);
        else
            % otherwise collect a 2D array of spike rates for all cells across samples
            %srSamples = sum(spikeHist,3);           % count spikes rather than mean spike rate
            srSamples = mean(spikeHist,3) / model.simTimeStep;
        end
        if ~isfield(m, 'srPoints') || isempty(m.srPoints) || size(m.srPoints,2) ~= m.historySize
            m.srPoints     = zeros(size(srSamples,1), m.historySize);
            m.srCount      = 0;               
            m.srCurrentPtr = 0;
        elseif m.srCurrentPtr + numSamples > m.historySize
            m.srCurrentPtr = 0;
        end
        m.srPoints(:, m.srCurrentPtr+1:m.srCurrentPtr+numSamples) = srSamples;
        m.srCount      = max(m.srCount, m.srCurrentPtr + numSamples);
        m.srCurrentPtr = mod(m.srCurrentPtr + numSamples, m.historySize);
    end
end


% Measure_Sparseness - Calculate various sparseness measures
%
% Lifetime sparseness (Vinje & Gallant, Treves & Rolls):
%     sparseness = E_p[ (1 - E_t[r]^2 / E_t[r^2]) / (1 - 1/N_t) ]
%
% Population sparseness:
%     sparseness = E_t[ (1 - E_p[r]^2 / E_p[r^2]) / (1 - 1/N_p) ]
%
% Where r = response (spike rate) of the neuron; E_t[] = expected value
% over time, calculated as a moving average over simulation time;
% E_p[] = expected value over the cell population; N is the number of cells.
%
% Activity sparseness (Berkes P et al, NIPS, 2009):
%     sparseness = 1 - n_t / N
%
% Where n_t = number of neurons with activity larger than threshold at time t.
% threshold is set to one stddev = 68th percentile
%
% Inputs:
%    m            = the measure descriptor (model.stats.measure.<measureName>)
%        sourceName   = name of cell group, e.g. 'cg_V1e'
%        numSamples   = number of sample to include in moving averages
%    model        = the model (read-only)
%
% Outputs:
%    m.timeSparseness     = lifetime sparseness measure
%    m.popSparseness      = population sparseness measure
%    m.activitySparseness = activity sparseness measure
%
% Q. (TODO) Since spikeRate is regulated, is E_t[r] == p?  what about E_p[r] ?
%--------------------------------------------------------------------------
function m = Measure_Sparseness( m, model, spikeHistory )

    assert(~isempty(spikeHistory), 'spikeHistory is required for sparsity measurement');

    % first-time variable initialization
    if ~isfield(m, 'initialized') || ~m.initialized
        if ~isfield(m, 'numAvgSamples')
            m.numAvgSamples = 2000;         % default # input samples to average
        end
        m.cellGroupId               = NetModel_FindElement(model, m.sourceName);
        m.biasedStats_max           = 0;
        m.r_ave_biased              = 0;
        m.r2_ave_biased             = 0;
        m.popSparseness_biased      = 0;
        m.activitySparseness_biased = 0;
        m.initialized               = true;
    end

    % precalculate some values
    % r(numCells,numNetworks) = mean spike rate (spikes/timeUnit) for each neuron x sample
    [numCells, numNetworks, numIterations] = size(spikeHistory{m.cellGroupId});
    r                 = mean(spikeHistory{m.cellGroupId}, 3) / model.simTimeStep;
    stats_eta         = 1 - exp(- numNetworks/m.numAvgSamples);
    m.biasedStats_max = (1-stats_eta) .* m.biasedStats_max + stats_eta;

    % compute lifetime sparseness
    % r_ave  = moving average spikeRate (r) for each cell
    % r2_ave = moving average spikeRate^2 (r^2) for each cell
    m.r_ave_biased    = (1 - stats_eta) * m.r_ave_biased  + stats_eta * mean(r,    2);
    m.r2_ave_biased   = (1 - stats_eta) * m.r2_ave_biased + stats_eta * mean(r.^2, 2);
    m.r_ave           = m.r_ave_biased  / m.biasedStats_max;
    m.r2_ave          = m.r2_ave_biased / m.biasedStats_max;
    N_effective       = m.numAvgSamples * m.biasedStats_max;
    m.timeSparseness  = mean( (1 - m.r_ave.^2 ./ m.r2_ave) / (1 - 1/N_effective) );

    % compute population sparseness
    popSparseness_this     = (1 - mean(r,1).^2 ./ mean(r.^2,1)) / (1 - 1/numCells);
    popSparseness_this(~isfinite(popSparseness_this)) = 0;    % remove NaNs
    m.popSparseness_biased = (1 - stats_eta) * m.popSparseness_biased ...
                                + stats_eta * mean(popSparseness_this);
    m.popSparseness        = m.popSparseness_biased / m.biasedStats_max;

    % compute activity sparseness
    % numActive(numNetworks) = number of cells with spikeRate > 1 standard deviation
    numActive                   = sum(bsxfun(@gt, r, std(r,1)), 1);
    activitySparseness_this     = 1 - numActive / numCells;
    m.activitySparseness_biased = (1 - stats_eta) * m.activitySparseness_biased ...
                                     + stats_eta * mean(activitySparseness_this);
    m.activitySparseness        = m.activitySparseness_biased / m.biasedStats_max;
end


% Measure_DeltaWeight - Track the moving average of the RMS weight change
%
% Inputs:
%    m            = the measure descriptor (model.stats.measure.<measureName>)
%        deltaInterval = the interval in time units for calculating dW (default = 1000)
%        windowSize    = the moving average window for RMS calculation in time units
%                        (default = deltaInterval)
%        cg(i)         = information collected on cell group #i (i is arbitrary)
%            name         = a human-readable name of this cell group, e.g. 'V1e'
%            cgId         = cell group id of this input block
%            dThresh   = the moving average RMS threshold change (units = dThresh / timeUnit)
%        ib(i)         = information collected on input block #i
%            name         = a human-readable name of this input block, e.g. 'V1e->V1i'
%            cgId         = cell group id of this input block
%            ibId         = cell input block id of this input block
%            srcId        = cell group id of the input source
%            dW        = the moving average RMS weight change (units = dW / timeUnit)
%        dW         = the moving average RMS weight change across all weight sets
%                        (units = dW / timeUnit)
%        ddW        = the derivative of dW
%        dThresh    = the moving average RMS spike threshold change across all cell groups
%                     (units = dW / timeUnit)
%        ddThresh   = the derivative of dThresh
%    model        = the model (read-only)
%--------------------------------------------------------------------------
function m = Measure_DeltaWeight( m, model, spikeHistory )

    cg = model.cellGroup;

    % initialize: generate a list of all active/enabled input blocks
    if ~isfield(m, 'initialized') || ~m.initialized
        if ~isfield(m, 'deltaInterval')
            m.deltaInterval = 2000;            % dW interval in time units
        end
        if ~isfield(m, 'windowSize')
            m.windowSize = m.deltaInterval;    % moving average window in time units
        end
        m.intervalWait = 0;
        m.cg           = [];
        m.ib           = [];
        for i = 1:numel(cg)
            if cg{i}.isExternal || isempty(cg{i}.inputBlock) || strcmp(cg{i}.cellType, 'disabled')
                continue;
            end
            m.cg(end+1).cgId         = i;
            m.cg(end).name           = cg{i}.name;
            m.cg(end).N              = cg{i}.numCells;
            m.cg(end).thresh0        = cg{i}.spikeThresh;
            m.cg(end).dThresh_biased = 0;
            m.cg(end).dThresh        = 0;
            ibIds = find( ~strcmp({cg{i}.inputBlock.connectionType}, 'disabled') );
            for k = ibIds
                m.ib(end+1).cgId    = i;
                m.ib(end).ibId      = k;
                m.ib(end).srcId     = cg{i}.inputBlock(k).srcId;
                m.ib(end).name      = sprintf('%s->%s', cg{i}.inputBlock(k).name, cg{i}.name);
                m.ib(end).N         = numel(cg{i}.inputBlock(k).weight);
                m.ib(end).W0        = cg{i}.inputBlock(k).weight;
                m.ib(end).dW_biased = 0;
                m.ib(end).dW        = 0;
            end
        end
        m.biasedStats_max = 0;
        m.dW              = 0;
        m.ddW             = 0;
        m.dThresh         = 0;
        m.ddThresh        = 0;
        m.initialized     = true;
    end

    % see if it is time to update the RMS measurement
    [numCells, numNetworks, numIterations] = size(spikeHistory{model.outputCellGroupId});
    numTimeUnits   = numNetworks * numIterations * model.simTimeStep;
    m.intervalWait = m.intervalWait + numTimeUnits;

    % every so often calculate the weight change and update the moving average rms
    if m.deltaInterval > 0 && m.intervalWait >= m.deltaInterval

        % update the moving average of the dW measures
        stats_eta         = 1 - exp(- numTimeUnits/m.windowSize);
        m.biasedStats_max = (1-stats_eta) .* m.biasedStats_max + stats_eta;

        % calculate threshold convergence for each cell group
        for i = 1:numel(m.cg)
            thresh = cg{m.cg(i).cgId}.spikeThresh;
            
            % calculate dW for this interval and update the W0 reference
            dThresh = thresh - m.cg(i).thresh0;
            if m.intervalWait == m.deltaInterval
                % simulation interval exactly matches the target interval boundary
                m.cg(i).thresh0 = thresh;
            else
                % otherwise interpolate to get consistency across samples
                dThresh         = dThresh .* (m.deltaInterval / m.intervalWait);
                m.cg(i).thresh0 = thresh - dThresh;
            end
            
            % calculate the rms dThresh (normalized to thresh magnitude) and add to moving average
            threshRms              = max( sqrt(mean( thresh(:).^2 )), .01);
            dThresh_this           = sqrt(mean( dThresh(:).^2 )) / threshRms;
            m.cg(i).dThresh_biased = (1-stats_eta) .* m.cg(i).dThresh_biased + stats_eta * dThresh_this;
            m.cg(i).dThresh        = m.cg(i).dThresh_biased ./ m.biasedStats_max;
        end

        % update average (weighted mean) across all threshold sets
        dThresh_prev = m.dThresh;
        m.dThresh    = sum([m.cg.dThresh] .* [m.cg.N]) / sum([m.cg.N]);    % weighted mean
        m.ddThresh   = m.dThresh - dThresh_prev;

        % calculate the weight convergence for each input block
        for i = 1:numel(m.ib)
            W = cg{m.ib(i).cgId}.inputBlock(m.ib(i).ibId).weight;
            
            % calculate dW for this interval and update the W0 reference
            dW = W - m.ib(i).W0;
            if m.intervalWait == m.deltaInterval
                % simulation interval exactly matches the target interval boundary
                m.ib(i).W0 = W;
            else
                % otherwise interpolate to get consistency across samples
                dW         = dW .* (m.deltaInterval / m.intervalWait);
                m.ib(i).W0 = W - dW;
            end
            
            % calculate the rms dW (normalized to W magnitude) and add to moving average
            Wrms              = max( sqrt(mean( W(:).^2 )), .01);
            dW_this           = sqrt(mean( dW(:).^2 )) / Wrms;
            m.ib(i).dW_biased = (1-stats_eta) .* m.ib(i).dW_biased + stats_eta * dW_this;
            m.ib(i).dW        = m.ib(i).dW_biased ./ m.biasedStats_max;
        end

        % update average (weighted mean) across all weight sets
        dW_prev = m.dW;
        m.dW    = sum([m.ib.dW] .* [m.ib.N]) / sum([m.ib.N]);    % weighted mean
        m.ddW   = m.dW - dW_prev;

        m.intervalWait = 0;
    end
end


% Measure_TimelineStats - Record network metrics over extended time
%
% Collects measurements which are maintaind in a rotating buffer
% (to reduce memory shifting and reallocation).
%
% Inputs:
%    m            = the measure descriptor (model.stats.measure.<measureName>)
%        metricExpr   = a cell array of MatLab expressions for each metric
%                       to track.  See NetModel_EvalExpression for details
%                       on automatic variables in the evaluation scope.
%        historySize  = number of historical samples to collect before recycling
%                       (default = 1000).
%    model        = the model (read-only)
%
% Outputs:
%    m            = updated measure descriptor
%        history       = matrix(numStats, historySize) of measurement samples
%        sampleNum     = array(1, historySize) indicating the sample counter
%                        at the time of measurement
%        numHistory    = number of historical points available
%        currentPtr    = pointer to the current value
%        bufferWrapped = has history wrapped around to beginning?
%--------------------------------------------------------------------------
function m = Measure_TimelineStats( m, model, spikeHistory )

    % first-time variable initialization
    if ~isfield(m, 'initialized') || ~m.initialized
        if ~isfield(m, 'historySize')
            m.historySize = 1000;
        end
        if ischar(m.metricExpr)
            m.metricExpr = { m.metricExpr };     % convert string to cell array
        end
        m.history       = zeros(numel(m.metricExpr), m.historySize);
        m.sampleNum     = zeros(1, m.historySize);
        m.currentPtr    = 0;
        m.bufferWrapped = false;
        m.initialized   = true;
    end

    % update history with new metric samples
    m.historyWrapped = (m.currentPtr == m.historySize);
    m.currentPtr     = mod(m.currentPtr, m.historySize) + 1;
    vals = NetModel_EvalExpression(model, m.metricExpr);
    for i = 1:numel(m.metricExpr)
        m.history(i, m.currentPtr) = vals{i};
    end
    m.sampleNum(m.currentPtr) = model.stats.trainingSampleCount;
end



% ======================================================================================
% ================================   PRINT ROUTINES   ==================================
% ======================================================================================

% PrintStatsLine - Print a line of statistics
%
function model = PrintStatsLine( model )

    if isempty(model.stats.print)
        return;
    end

    if ~model.stats.areTitlesPrinted
        
        % determine column display order, store in model.stats.printNames
        printList   = fieldnames(model.stats.print);
        colPosition = zeros(1, numel(printList));
        for i = 1:numel(printList)
            p = model.stats.print.(printList{i});
            if isfield(p, 'position')
                colPosition(i) = p.position;
            end
        end
        [sortedColPositions, sortedColIndices] = sort(colPosition);
        model.stats.printNames = printList( sortedColIndices(sortedColPositions > 0) );
        
        % print column titles first time through
        for i = 1:numel(model.stats.printNames)
            p = model.stats.print.(model.stats.printNames{i});
            if ~isfield(p, 'title') || isempty(p.title)
                p.title = model.stats.printNames{i};     % default title is column field name
            end
            if ~isfield(p, 'formatTitle') || isempty(p.formatTitle)
                p.formatTitle = '%8s';      % default title format
            end
            if ~isfield(p, 'formatData') || isempty(p.formatData)
                p.formatData = '%8.4f';     % default data format
            end
            fprintf(p.formatTitle, p.title);
            model.stats.print.(model.stats.printNames{i}) = p;
        end
        if numel(model.stats.printNames) > 0
            fprintf('\n');
        end
        model.stats.areTitlesPrinted = true;
    end


    % print a single line of statistics data
    for i = 1:numel(model.stats.printNames)
        p = model.stats.print.(model.stats.printNames{i});
        if isfield(p, 'measureName') && ~isempty(p.measureName)
            if ~isfield(model.stats.measure, p.measureName)
                error('Error processing stats.print.%s: measureName "%s" not found', ...
                        model.stats.printNames{i}, p.measureName);
            end
            m = model.stats.measure.(p.measureName);
            switch m.measureType
            case 'resError'
                val = m.rmsResErr;
            case 'correlation'
                val = m.rmsCorrelation;
            case 'spikeRate'
                val = m.popSpikeRate;
            case 'deltaWeight'
                val = m.dW;
            otherwise
                error('Unknown measure type "%s"', m.measureType);
            end
        elseif isfield(p, 'builtin') && ~isempty(p.builtin)
            switch (p.builtin)
            case 'sampleCount'
                val = model.stats.trainingSampleCount;
            otherwise
                error('Unknown builtin field "%s"', p.builtin);
            end
        elseif isfield(p, 'matlabExpr') && ~isempty(p.matlabExpr)
            if isfield(p, 'sourceName') && ~isempty(p.sourceName)
                val = NetModel_EvalExpression(model, p.matlabExpr, p.sourceName);
            else
                val = NetModel_EvalExpression(model, p.matlabExpr);
            end
        else
            error('unrecognized print column type at "%s"', model.stats.printNames{i});
        end
        fprintf(p.formatData, val);
        model.stats.print.(model.stats.printNames{i}) = p;
    end
    if numel(model.stats.printNames) > 0
        fprintf('\n');
    end
end


% ======================================================================================
% ===============================   FIGURE ROUTINES   ==================================
% ======================================================================================


% InitializeFigure - Perform first-time initialization of a figure
%
%--------------------------------------------------------------------------
function [fig, model] = InitializeFigure( fig, model )

    switch fig.plotType

    case 'STA'
        % if measureName is missing, fill in from sourceName
        if ~isfield(fig, 'measureName') || isempty(fig.measureName)
            assert(isfield(fig, 'sourceName') && ~isempty(fig.sourceName), ...
                'Either sourceName or measureName must be provided when creating STA figure');
            [model, measureName] = FindOrCreateSTA(model, fig.sourceName, 'STA');
            fig.measureName = measureName;
        end

    case 'STA_movie'
        % if measureName is missing, fill in from sourceName
        if ~isfield(fig, 'measureName') || isempty(fig.measureName)
            assert(isfield(fig, 'sourceName') && ~isempty(fig.sourceName), ...
                'Either sourceName or measureName must be provided when creating STA figure');
            [model, measureName] = FindOrCreateSTA(model, fig.sourceName, 'timeSTA');
            fig.measureName = measureName;
            if isfield(fig, 'numFrames')
                model.stats.measure.(measureName).numFrames = fig.numFrames;
            end
            if isfield(fig, 'timeInterval')
                model.stats.measure.(measureName).timeInterval = fig.timeInterval;
            end
        end

    case 'STAWRatioHistogram'
        if ~isfield(fig, 'measureName') || isempty(fig.measureName)
            assert(isfield(fig, 'sourceName') && ~isempty(fig.sourceName), ...
                'Either sourceName or measureName must be provided when creating STAWRatioHistogram');
            [model, measureName] = FindOrCreateSTA(model, fig.sourceName, 'STA');
            fig.measureName = measureName;
        end

    case 'rStdHistogram'
        if ~isfield(fig, 'measureName') || isempty(fig.measureName)
            fig.measureName = 'resErr_out';
        end
        model.stats.measure.(fig.measureName).stdPoints = [];    % turn on data collection

    case 'spikeRateHistogram'
        if ~isfield(fig, 'measureName') || isempty(fig.measureName)
            fig.measureName = 'spikeRate_out';
        end
        if ~isfield(model.stats.measure, fig.measureName)
            m = struct();
            m.measureType = 'spikeRate';
            m.sourceName  = 'cg_output';
        else
            m = model.stats.measure.(fig.measureName);
        end
        if ~isfield(fig, 'dimension')
            fig.dimension = 'perSample';           % default is perSample spike rate data collection
        end
        if strcmp(fig.dimension, 'perSample')
            if isfield(fig, 'historySize')
                m.historySize = fig.historySize;
            elseif ~isfield(m, 'historySize')
                m.historySize = model.stats.updateFiguresEvery * model.displayTimeScale ...
                        / model.numIterationsPerSample;
            end
        end
        if isfield(fig, 'srPointsPop')
            m.srPointsPop = fig.srPointsPop;
        end
        m.initialized = false;
        model.stats.measure.(fig.measureName) = m;

    case 'timelineStats'
        if ~isfield(fig, 'measureName') || isempty(fig.measureName)
            fig.measureName = 'timelineStats';      % standard measure name
            if ~isfield(model.stats.measure, fig.measureName)
                model.stats.measure.(fig.measureName).measureType = 'timelineStats';
                model.stats.measure.(fig.measureName).metricExpr = fig.metricExpr;
                fig = rmfield(fig, 'metricExpr');
                if isfield(fig, 'historySize')
                    model.stats.measure.(fig.measureName).historySize = fig.historySize;
                    fig = rmfield(fig, 'historySize');
                end
            end
        end

    case 'composite'
        for i = 1:numel(fig.subFigure)
            subfig = fig.subFigure{i};
            if ~isempty(subfig)
                 [subfig, model]  = InitializeFigure(subfig, model);
                 fig.subFigure{i} = subfig;
            end
        end
    end
end


% FindOrCreateSTA - Find or create an STA measure
function [model, measureName] = FindOrCreateSTA( model, sourceName, measureType )

    measureName = [];
    
    % if sourceName has more than one part, keep only first part
    dots = strfind(sourceName, '.');
    if ~isempty(dots)
        sourceName = sourceName(1:dots(1)-1);
    end
        
    % search for a pre-existing STA measure on the right targets
    measureList = fieldnames(model.stats.measure);
    for j = 1:numel(measureList)
        m = model.stats.measure.(measureList{j});
        if strcmp(m.measureType, measureType) && strcmp(m.sourceName, sourceName) ...
                && (~isfield(m, 'referantName') || strcmp(m.referantName, 'input'))
            measureName = measureList{j};
            break;
        end
    end

    % if no STA measure found, create a default one
    if isempty(measureName)
        measureName = [measureType '_' sourceName(4:end)];
        if isfield(model.stats.measure, measureName)
            error('Can''t create STA measure; measure with name "%s" already exists', measureName);
        end
        model.stats.measure.(measureName).measureType = measureType;
        model.stats.measure.(measureName).sourceName  = sourceName;
    end
end

