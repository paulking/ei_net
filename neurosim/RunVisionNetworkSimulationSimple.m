% RunVisionNetworkSimulationSimple - Run a spiking circuit simulation
%
% Build and train a spiking circuit model that works on natural images.
%
% Note that the printed log of network statistics during training
% is generated from NetModel_Stats, called by NetModel_UpdateFast.
%
% Usage:
%    model = RunVisionNetworkSimulationSimple(params)
%    model = RunVisionNetworkSimulationSimple(params, model)
%
% Inputs:
%    params        = configuration parameters:
%        numInputSamples       = number of image patch samples to train on
%                                (default = 10000)
%        printSummary          = print column headers and time elapsed summary
%                                (default = true)
%        model                 = parameters for constructing model. Required if
%                                generating a new model. Ignored if an already-constructed
%                                model is provided as a second parameter.
%    model         = the already-constructed model to use. Providing this
%                    causes params.model to be ignored. (optional)
%
% Output:
%    model         = the network model ending state
%
% See also:
%    NetModel_Init, CellGroup_Init, CellGroup_AddInput, NetModel_Stats, NetModel_Plot
%
% Created:   10/20/11, Paul King
%--------------------------------------------------------------------------
function model = RunVisionNetworkSimulationSimple( params, model )

    if nargin < 2
        model = [];
    end

    % initialize default parameter values
    defaultValues.numInputSamples       = 10000;
    defaultValues.visualInput           = 'IMAGES.mat';
    defaultValues.printSummary          = true;       % print summary report at end
    defaultValues.model.modelType       = 'EINet';
    defaultValues.model.stats.measure.resErr_out.measureType = 'resError';
    defaultValues.model.stats.measure.resErr_out.sourceName  = 'cg_output.in_input';
    params = ApplyDefaultValues(params, defaultValues);


    % ================  LOAD AND PREPROCESS IMAGE INPUT  ===============

    % load the visual input data set
    if strcmp(params.visualInput, 'IMAGES.mat')
        % load a 512x512 x 10 matrix of 10 images, 512x512 pixels each
        % This is from the support materials for:
        %      Olshausen BA, Field DJ (1997)  Sparse Coding with an Overcomplete
        %          Basis Set: A Strategy Employed by V1?  Vision Research, 37: 3311-3325. 
        load('IMAGES.mat');

        % crop away a 4-pixel margin at the frame edge which contains FFT-induced error
        frameMargin = 4;
        IMAGES = IMAGES(1+frameMargin:end-frameMargin, 1+frameMargin:end-frameMargin, :);

        numImages = size(IMAGES, 3);
        imageSize = size(IMAGES, 1);

    elseif isnumeric(params.visualInput)
        fprintf('using explicitly provided set of %d input samples\n', ...
                size(params.visualInput,2) );

    else
        error('visual input source not supported.');
    end


    % ===================  BUILD THE NETWORK MODEL  ====================

    % build the circuit model if not already built
    if isempty(model)
        if isempty(params.model.modelType)
            error('a modelType must be specified to create and simulate a model');
        end
        switch params.model.modelType
        case 'EINet'
            model = NetModel_InitEINet(params.model);
        case 'V1eV1i'
            model = NetModel_InitV1eV1i(params.model);
        otherwise
            fprintf('WARNING: unrecognized modelType %s\n', params.model.modelType);
            model = NetModel_Init(params.model);
        end
    end

    % reset reporting
    if params.printSummary
        model.stats.areTitlesPrinted             = false;  % force print of column headers
        model.stats.measure.resErr_out.minResErr = Inf;    % reset minResErr tracker
    end

    % get cell group ids and sizes from model
    i_cg_input        = FindByName(model.cellGroup, 'input');
    i_cg_output       = model.outputCellGroupId;
    numCells_output   = model.cellGroup{i_cg_output}.numCells;


    % =====================  RUN THE SIMULATION  =====================

    startETime = clock();

    % iterate through a sequence of sample batches
    sampleCount = 0;
    while sampleCount < params.numInputSamples

        % =================  PREPARE A BATCH OF INPUT SAMPLES  ==================

        if isnumeric(params.visualInput)
            % if explicit input data provided, use that
            sampleIds = 1 + mod( sampleCount:(sampleCount+model.numSamplesPerBatch-1), ...
                    size(params.visualInput,2) );
            inputData = params.visualInput(:, sampleIds);
        else
        
            % otherwise, select a set of random image patches from the image set (normal condition)
            patchDims = model.inputDims;
            inputData = zeros(prod(patchDims), model.numSamplesPerBatch);
            for j = 1:model.numSamplesPerBatch
                image_idx  = ceil( numImages * rand() );
                r          = ceil( (imageSize-patchDims(1)) * rand() );
                c          = ceil( (imageSize-patchDims(2)) * rand() );
                imagePatch = IMAGES(r:r+patchDims(1)-1, c:c+patchDims(2)-1, image_idx);
                inputData(:,j) = imagePatch(:);
            end

            % normalize to zero mean and unit variance
            inputData = bsxfun(@minus, inputData, mean(inputData));
            inputData = bsxfun(@times, inputData, 1 ./ std(inputData));
        end
        
        % store input data in model.stats for use in performance analysis later
        model.snapshot.inputData = inputData;
        sampleCount = sampleCount + model.numSamplesPerBatch;


        % =====================   TRAIN NETWORK   ======================

        % process the image(s) through the spiking network model

        % set the input image
        initialState = struct();
        initialState.y{i_cg_input} = inputData * (model.inputScale * model.simTimeStep);

        % simulate the network in batch
        [model, finalState] = NetModel_UpdateFast(model, model.numIterationsPerSample, initialState);
        %numSpikes = sum(finalState.y_history{i_cg_output}, 3);
    end


    % ========================   WRAP UP   ========================

    % print summary statistics
    if params.printSummary
        fprintf('Minimum residual error achieved (min_resErr) =%7.4f\n', ...
                model.stats.measure.resErr_out.minResErr);

        % summarize elapsed time required
        elapsedSeconds = etime(clock(), startETime);
        elapsedMinutes = floor(elapsedSeconds / 60);
        elapsedSeconds = elapsedSeconds - elapsedMinutes * 60;
        fprintf('Elapsed time to run simulation: %.0f minute(s), %.2f second(s)\n', ...
                elapsedMinutes, elapsedSeconds);
    end
end
