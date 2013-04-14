% RunParameterizedExperiment - Run a set of experiments varying one or two parameters
%
% Inputs:
%    simParams    = the simulation (and model) parameters
%    params       = configuration parameters:
%        propNames         = cell array{numProps} of one or two parameter names to vary,
%                            e.g. { 'model.cg_V1e.numCells', 'model.cg_V1i.numCells' }.
%                            (required).
%        propValues        = cell array{numProps} of one or two value sets for each
%                            parameter. Each value set is an array(numValues). (required)
%        actionPre         = string or cell array of statement(s) to execute
%                            prior to simulation. (default = [])
%        actionPost        = string or cell array of statement(s) to execute
%                            after results measurement simulation. (default = [])
%        captureSparseness = cell array of cell group names for capturing
%                            sparseness information, if any. (default = [])
%
% Outputs:
%    results      = experiment results
%
% Example:
%    results = RunParameterizedExperiment(simParams, params);
%
% Created:   1/6/2013, Paul King
%--------------------------------------------------------------------------
function results = RunParameterizedExperiment( simParams, params )

    % initialize and apply default parameter values
    if nargin < 2
        params = struct();
    end
    defaultValues = struct();
    defaultValues.numTestIterations    = 10000;
    defaultValues.captureSparseness    = [];
    defaultValues.measureAvgSamples    = 5000;
    defaultValues.sparseAvgSamples     = 5000;
    defaultValues.propNames            = [];
    defaultValues.propValues           = [];
    defaultValues.actionPre            = [];
    defaultValues.actionPost           = [];

    params = ApplyDefaultValues(params, defaultValues);

    assert(~isempty(params.propNames),  'parameter propNames must be provided');
    assert(~isempty(params.propValues), 'parameter propValues must be provided');

    % convert optional strings to cell arrays for processing consistency
    if ischar(params.captureSparseness)
        params.captureSparseness = { params.captureSparseness };
    end
    if ischar(params.actionPre)
        params.actionPre = { params.actionPre };
    end
    if ischar(params.actionPost)
        params.actionPost = { params.actionPost };
    end


    % ==================   INITIALIZE VARIABLES   =================

    fprintf('Running a matrix of simulation experiments, varying the following parameters:\n');
    for k = 1:numel(params.propNames)
        fprintf('    %s: ', params.propNames{k});
        for j = 1:numel(params.propValues{k})
            if j > 1;  fprintf(', ');  end
            fprintf('%s', convertToString(params.propValues{k},j));
        end
        fprintf('\n');
    end
    fprintf('\n');

    % turn on sparseness tracking, if requested
    for k = 1:numel(params.captureSparseness)
        cgName      = params.captureSparseness{k};
        measureName = ['sparse_' cgName];
        simParams.model.stats.measure.(measureName).measureType   = 'sparseness';
        simParams.model.stats.measure.(measureName).sourceName    = ['cg_' cgName];
        simParams.model.stats.measure.(measureName).numAvgSamples = params.sparseAvgSamples;
    end

    % initialize 'results' struct
    results = struct();
    results.simParams        = simParams;
    results.experimentParams = params;
    results.propNames        = params.propNames;
    results.propValues       = params.propValues;

    % special case: only one property being varied
    if numel(params.propNames) == 1
        % add a dummy second property to make the nested loop work
        params.propNames{2}  = 'dummy_';     % a dummy property for simParams
        params.propValues{2} = 0;            % a dummy value to be assigned
        results.description  = sprintf('experiment varying %s', params.propNames{1});
    else
        results.description  = sprintf('experiment varying %s and %s', ...
                                     params.propNames{1}, params.propNames{2});
    end


    % ==================   RUN MULTIPLE SIMULATIONS   =================

    startETime = clock();

    for i = 1:numel(params.propValues{1})
        for j = 1:numel(params.propValues{2})

            % set custom parameter values
            eval(['simParams.' params.propNames{1} ' = ' convertToString(params.propValues{1},i) ';'] );
            if ~strcmp(params.propNames{2}, 'dummy_')
                eval(['simParams.' params.propNames{2} ' = ' convertToString(params.propValues{2},j) ';'] );
            end

            % execute any optional caller-provided preprocessing statements
            for k = 1:numel(params.actionPre)
                eval([params.actionPre{k} ';']);
            end

            % run the simulation
            randreset();  model = RunVisionNetworkSimulation(simParams);

            % save an example full parameter set that the model was initialized with.
            % (this includes the framework provided defaults, which are otherwise unrecorded)
            if i == 1 && j == 1
                results.modelInitParams = model.initParams;
            end

            % measure results with learning disabled in a shorter test simulation, higher sample size
            model.learningRate                           = 0;          % turn off learning
            model.stats.measure.resErr_out.numAvgSamples = params.measureAvgSamples;
            model.stats.measure.corr_out.numAvgSamples   = params.measureAvgSamples;
            model = RunVisionNetworkSimulation({'numInputSamples',params.numTestIterations}, model);

            % save results to struct
            results.finalResErr(i,j) = model.stats.measure.resErr_out.rmsResErr;
            results.finalCorr(i,j)   = model.stats.measure.corr_out.rmsCorrelation;

            % capture sparseness information as well, if enabled
            for k = 1:numel(params.captureSparseness)
                cgName      = params.captureSparseness{k};
                measureName = ['sparse_' cgName];
                results.(measureName).time(i,j)     = model.stats.measure.(measureName).timeSparseness;
                results.(measureName).pop(i,j)      = model.stats.measure.(measureName).popSparseness;
                results.(measureName).activity(i,j) = model.stats.measure.(measureName).activitySparseness;
            end

            % execute any optional caller-provided postprocessing statements
            for k = 1:numel(params.actionPost)
                eval([params.actionPost{k} ';']);
            end
        end
    end


    % =========================   WRAP UP   =========================

    elapsedSeconds = etime(clock(), startETime);
    elapsedMinutes = floor(elapsedSeconds / 60);
    elapsedSeconds = elapsedSeconds - elapsedMinutes * 60;
    fprintf('Elapsed time to run experiment matrix: %.0f minute(s), %.2f second(s)\n\n', ...
            elapsedMinutes, elapsedSeconds);
end


function str = convertToString(valueList, idx)

    if iscell(valueList)
        val = valueList{idx};
    else
        val = valueList(idx);
    end

    if ischar(val)
        str = [''' val '''];
    elseif isnumeric(val)
        str = mat2str(val);
    elseif islogical(val)
        if val; fprintf('true'); else fprintf('false');  end
    else
        error('unrecognized type %s for param value', class(val));
    end
end
