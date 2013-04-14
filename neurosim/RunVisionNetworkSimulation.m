% RunVisionNetworkSimulation - Run a spiking circuit simulation
%
% Build and train a spiking circuit model that works on natural images.
%
% Note that the printed log of network statistics during training
% is generated from NetModel_Stats, called by NetModel_UpdateFast.
%
% Usage:
%    model = RunVisionNetworkSimulation(params)
%    model = RunVisionNetworkSimulation(params, model)
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
function model = RunVisionNetworkSimulation( params, model )

    if nargin < 2
        model = [];
    end

    % initialize default parameter values
    defaultValues.numInputSamples       = 10000;
    defaultValues.visualInput           = 'IMAGES.mat';
    defaultValues.printSummary          = true;       % print summary report at end
    defaultValues.model.modelType       = [];         % no default model type
    defaultValues.model.inputPreprocess.csFilter             = false;
    defaultValues.model.autoAnneal.wait                      = 0;
    defaultValues.model.autoAnneal.minLearningRate           = .01;
    defaultValues.model.autoAnneal.tolerance                 = .001;
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

    elseif strcmp(params.visualInput, 'IMAGES_RAW.mat')
        % load unwhitened (natural) images (same format as above)
        load('IMAGES_RAW.mat');
        IMAGES    = IMAGESr;
        numImages = size(IMAGES, 3);
        imageSize = size(IMAGES, 1);

    elseif isnumeric(params.visualInput)
        fprintf('using explicitly provided set of %d input samples\n', ...
                size(params.visualInput,2) );

    elseif isstruct(params.visualInput)
        error('visual input of type MatLab struct is not (yet) supported.');
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
        case 'SAILnet'
            model = SAILnet_InitNetModel(params.model);
        case 'SAILnet_V1eV1i'
            model = SAILnet_InitModelV1eV1i(params.model);
        case 'V1eV1i'
            model = NetModel_InitV1eV1i(params.model);
        case 'EINet'
            model = NetModel_InitEINet(params.model);
        case 'V1TwoLayer'
            model = NetModel_InitV1TwoLayer(params.model);
        case 'RetV1eV1i'
            model = NetModel_InitRetV1eV1i(params.model);
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

    assert(~model.simNoReset || strcmp(model.simMethod,'cellgroupFastBatch'), ...
            'simNoReset mode is only supported in cellgroupFastBatch mode');

    % prepare for center-surround preprocessing (if enabled)
    % output: csFilterMatrix(numInputs, numPixels)  - filter bank for converting input
    if model.inputPreprocess.csFilter
        csFilterMatrix = GenerateCSFilterArray(model.inputDims, prod(model.inputDims), ...
                model.inputPreprocess)';

        % DEBUG: show center-surround patch filters
        csFilter = reshape(csFilterMatrix', model.inputDims(1), model.inputDims(2), []);
        ShowFigure(5, 'center-surround patch filters');
        ShowImagePatchGrid(csFilter);
    end


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


        % =================  PREPROCESS INPUT (OPTIONAL)  ==================

        % preprocess input: possibly run through center-surround filter
        if model.inputPreprocess.csFilter
            inputData = csFilterMatrix * inputData;
            inputData = bsxfun(@minus, inputData, mean(inputData));
            inputData = bsxfun(@times, inputData, 1./std(inputData));
        end
        
        % preprocess input: split channel into two sign-inverted channels?
        if model.inputPreprocess.splitInvert
            inputData = cat(1, inputData, -inputData);
        end

        % preprocess input: rectify input into the positive range
        % (soft rectification is achieved by adding a Laplace distribution
        % to hard recification, resulting in a smooth differentiable function)
        % note these alternatives to soft rectification:
        %    alpha * exp(...)/2 + ...    (less influence of negative range)
        %    max(0, exp(...)/2 - beta)   (negative domain hits zero axis)
        switch model.inputPreprocess.rectify
        case 'hard'
            inputData = max(inputData,0);                            % hard rectification
        case 'soft'
            inputData = exp(-abs(inputData))/2 + max(inputData,0);   % soft rectification
        case 'sigmoid'
            inputData = 1 ./ (1+exp(-inputData));                    % sigmoid to [0,1] range
        case {'none', ''}
            % do nothing
        otherwise
            error('unknown rectification strategy %s', model.inputPreprocess.rectify);
        end

        % DEBUG: show input image patches
        %imagePatches = reshape(inputData, model.inputDims(1), model.inputDims(2), []);
        %ShowFigure(6, 'input image patches');
        %ShowImagePatchGrid(imagePatches);

        % the number of input data elements must match the number of input cells
        assert(size(inputData,1) == model.cellGroup{i_cg_input}.numCells);


        % =====================   TRAIN NETWORK   ======================

        % process the image(s) through the spiking network model
        switch model.simMethod

        case 'cellgroupFastBatch'
            % execute cellgroup model in batch using optimized NetModel_UpdateFast()

            % set the input image
            initialState = struct();
            if model.simNoReset && sampleCount > model.numSamplesPerBatch
                % special case simNoReset mode: copy state if after the first iteration
                initialState.u     = finalState.u;
                initialState.y_ave = finalState.y_ave;
                initialState.y     = finalState.y;
            end
            initialState.y{i_cg_input} = inputData * (model.inputScale * model.simTimeStep);

            % simulate the network in batch
            [model, finalState] = NetModel_UpdateFast(model, model.numIterationsPerSample, initialState);
            %numSpikes = sum(finalState.y_history{i_cg_output}, 3);
            
        otherwise
            %numSpikes = zeros(numCells_output, model.numSamplesPerBatch);
            for j = 1:model.numSamplesPerBatch
                % reset state of all neurons (although the network works fine without any of this)
                for i = 1:numel(model.cellGroup)
                    model.cellGroup{i}.potential(:)         = 0;
                    model.cellGroup{i}.spikedRecently(:)    = 0;
                    model.cellGroup{i}.spikeRate.instant(:) = 0;
                end

                % set the input image
                % (spikedRecently is used for cell input, in same units as a spike)
                % (spikeRate.instant is used for learning, in units spikes-per-iteration)
                model.cellGroup{i_cg_input}.spikedRecently(:) = inputData(:,j) * (model.inputScale * model.simTimeStep);
                model.cellGroup{i_cg_input}.spikeRate.instant = model.cellGroup{i_cg_input}.spikedRecently;

                % simulate several iterations of the network using the appropriate method
                switch model.simMethod

                case 'cellgroupFast'
                    [model, finalState] = NetModel_UpdateFast(model, model.numIterationsPerSample);
                    %spikeRecord = reshape(finalState.y_history{i_cg_output}, numCells_output, model.numIterationsPerSample);

                case 'cellgroup'
                    %spikeRecord = zeros(numCells_output, model.numIterationsPerSample);
                    for i = 1:model.numIterationsPerSample
                        model = NetModel_Update(model);
                        %spikeRecord(:,i) = (model.cellGroup{i_cg_output}.spikedRecently == 1);
                    end

                otherwise
                    error('unrecognized simulation method ''%s''', model.simMethod);
                end

                %numSpikes(:,j) = sum(spikeRecord, 2);
            end
        end

        % auto-anneal feature: reduce model learning rate incrementally
        if model.autoAnneal.wait > 0
            model.autoAnneal.wait = model.autoAnneal.wait - model.numSamplesPerBatch;
            if model.autoAnneal.wait <= 0
                minDelta = model.autoAnneal.tolerance * model.learningRate;
                if model.stats.measure.dW.ddW > -minDelta && model.stats.measure.dW.ddThresh > -minDelta
                    model.learningRate    = model.learningRate / 2;
                    model.autoAnneal.wait = model.stats.measure.dW.windowSize * 10 / model.numTimeUnitsPerSample;
                    if model.learningRate < model.autoAnneal.minLearningRate
                        break;      % learning rate is too small; stop training
                    end
                else
                    model.autoAnneal.wait = 1;
                end
            end
        end

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
