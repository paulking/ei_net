% NetModel_UpdateFast - Execute many optimized steps of a spiking network simulation
%
% This optimized routine performs the same function as NetModel_Update
% but up to 40x faster. It achieves additional speed in the following way:
%    - performs the full network simulation in local variables
%    - performs many network iterations in a single call (optional)
%    - simulates multiple networks simultaneously in parallel (optional)
%    - can use average firing rates rather than spike-timing plasticity (optional)
%    - calculates weight updates in batch using matrix multiply
%
% The following constants define the size of the simulation:
%    numNetworks   = how many identical networks (but with different inputs
%                    and internal state) to simulate in parallel. When 100
%                    input samples are processed simultaneously, 100 clones
%                    of the full network are created and simulated. This
%                    allows MatLab to do substantial (up to 40x) optimizations
%                    with matrix math and multi-core parallel computation.
%    numIterations = how many simulation time steps to iterate through
%
% Each model simulated has the following basic structure:
%    numCellGroups  = the number of separate neuron populations in the model.
%                     Each neuron population has its own behavior and connectivity
%                     parameters
%    for each CellGroup:
%        numCells       = the number of identical neurons in that cell group
%        numInputBlocks = how many input types (from other cell groups) are received
%        for each InputBlock:
%            numInputs   = how many inputs of that type does each cell receive
%            weight      = matrix of connection weights from all input
%                          cells to all output cells
%            dweight     = accumulated (but not applied) weight matrix changes
%
% If the initialState input parameter is provided, then the network state
% contained in initialState will be used instead of the state contained
% in the model's cell groups. Each state field (u, y, y_ave) can optionally
% contain a matrix of multiple columns, in which case each column represents
% a separate parallel network which will be simulated in parallel and
% in batch for faster performance. Any variables or cell array values that
% are null or missing will be assumed to be initialized to zero.
%
% Note: the spikedRecently feature is not fully accurately simulated.
% Also, there are various restrictions and assumptions.
%
% Usage:
%    model = NetModel_UpdateFast(model, numIterations)
%    model = NetModel_UpdateFast(model, numIterations, initialState)
%    [model, finalState] = NetModel_UpdateFast(model, numIterations, initialState)
%
% Inputs:
%    model         = the simulation model state
%    numIterations = the number of simulation iterations to execute
%    initialState  = struct containing the initial cell group state for batch
%                    processing. any variables or values missing will be filled
%                    in with zeros. If the initialState struct is not provided,
%                    the state in the cell groups will be used instead. (optional)
%        u             = cell array(numCellGroups,1) of matrix(numCells,numNetworks) 
%                        representing initial cell potentials. (optional)
%        y             = cell array(numCellGroups,1) of matrix(numCells,numNetworks)
%                        representing the last spike output of the cells. Alternatively,
%                        y can be a matrix(numCells,numNetworks,numIterations), in
%                        which case y represents the input arriving from an external
%                        cell group to use for each iteration. (optional)
%        y_ave         = cell array(numCellGroups,1) of matrix(numCells,numNetworks) 
%                        representing the most recent spike rate moving average
%                        output of the cells. (optional)
%
% Output:
%    model           = the updated simulation model state
%    finalState      = the final state of the network (optional)
%        u               = matrix(numCells,numNetworks): ending cell potentials
%        y               = matrix(numCells,numNetworks): last spike status
%        y_ave           = matrix(numCells,numNetworks): last iteration
%                          of running average
%        y_history       = matrix(numCells,numNetworks,numIterations):
%                          spike history for all iterations and all networks
%
% Created:   09/03/11, Paul King
%--------------------------------------------------------------------------
function [model, finalState] = NetModel_UpdateFast( model, numIterations, initialState )

    % TODO can't return y_history for external (input) cells because the timing
    % from initialState.y is shifted early by one iteration (no output computed)

    if nargin < 3
        initialState = [];
    end

    cg = model.cellGroup;

    % check all assumptions on which this faster implementation is based
    for i = 1:numel(cg)
        if cg{i}.isExternal
            continue;
        end
        assert(cg{i}.targetSpikeRate > 0);
        for b = 1:numel(cg{i}.inputBlock)
            assert(isempty(cg{i}.inputBlock(b).inputIndices), ...
                'inputIndices not supported; only fully connected models');
        end
    end


    % =============  SLOW DOWN SIMULATION TO TARGET FRAME RATE  ============

    % pause if necessary in order to slow down to the desired frame rate
    % model.etime_offset = number of second behind we are (negative = fast)
    if model.fpsMax > 0
        if isfield(model, 'etime_t0')
            elapsedTime = etime(clock(), model.etime_t0);
            additionalWaitTime = 1/model.fpsMax - elapsedTime - model.etime_offset;
            if additionalWaitTime > 0
                pause(additionalWaitTime);
                elapsedTime = etime(clock(), model.etime_t0);
            end
            model.etime_offset = elapsedTime + model.etime_offset - 1/model.fpsMax;
            model.etime_offset = min(model.etime_offset, 2/model.fpsMax);
        else
            model.etime_offset = 0;         % initialize first time through
        end
        model.etime_t0 = clock();
    end


    % ======================  INITIALIZE NETWORK STATE  ======================

    % variables initialized here:
    %   u0{numCG}(numCells)      = initial cell potential
    %   y0{numCG}(numCells)      = initial cell spike output state (0 or 1)
    %   y_ave0{numCG}(numCells)  = initial moving average spike rate
    % 
    %   numNetworks              = how many networks are being simulated?
    %   totalNumIterations       = total iterations across all networks
    %   fullInputProvided(numCG) = was time-series input provided for this cell group?
    %   fullInputCellGroupIds    = cell group IDs providing full input

    % initialize boolean array for signaling time-series input special case
    % (if true, input for all iterations has been provided for an externally sourced cell group)
    fullInputProvided     = false(numel(cg), 1);
    fullInputCellGroupIds = [];

    if ~isempty(initialState)

        % if initialState was provided, use it to set up u, y, y_ave
        numCellGroups = numel(cg);
        
        % make sure all initialState fields are fully initialized
        if ~isfield(initialState, 'u') || numel(initialState.u) < numCellGroups
            initialState.u{numCellGroups} = [];
        end
        if ~isfield(initialState, 'y') || numel(initialState.y) < numCellGroups
            initialState.y{numCellGroups} = [];
        end
        if ~isfield(initialState, 'y_ave') || numel(initialState.y_ave) < numCellGroups
            initialState.y_ave{numCellGroups} = [];
        end

        % determine numNetworks (the number of network copies being simulated in parallel)
        % This is found by checking the number of columns of the state matrices.
        % Verify that all matrices have the same number of columns
        for i = 1:numCellGroups
            numNetworks_u(i)     = size(initialState.u{i}, 2);
            numNetworks_y(i)     = size(initialState.y{i}, 2);
            numNetworks_y_ave(i) = size(initialState.y_ave{i}, 2);
        end
        sizes = [numNetworks_u; numNetworks_y; numNetworks_y_ave];
        numNetworks = max(sizes(:));
        if ~ all(sizes == 0 | sizes == numNetworks)
            error('all initial state variables must have the same size or be null');
        end
        if numNetworks == 0
            error('if initial state is provided, values must be provided also');
        end

        % initialize state variables from initialState
        u0     = initialState.u;
        y0     = initialState.y;
        y_ave0 = initialState.y_ave;
        for i = 1:numCellGroups
            if numNetworks_u(i) == 0
                u0{i} = zeros(cg{i}.numCells, numNetworks);
                % initialize membrane potentials to random values to reduce time synchronization
                if strcmp(model.simInitPotentials, 'random')
                    u0{i} = cg{i}.spikeThresh * ones(1, numNetworks);
                    u0{i} = u0{i} .* rand(size(u0{i}));
                end
            end
            if numNetworks_y(i) == 0
                y0{i} = zeros(cg{i}.numCells, numNetworks);
            elseif size(y0{i},3) > 1
                % special case: initial value provided for all iterations; use the first one
                assert(size(y0{i},3) == numIterations, ...
                        '3rd dimension of initialState.y must match numIterations');
                y0{i} = y0{i}(:,:,1);
                fullInputProvided(i)         = true;
                fullInputCellGroupIds(end+1) = i;
            end
            if numNetworks_y_ave(i) == 0
                % if no incoming running average provided, initialize to zeros
                y_ave0{i} = zeros(cg{i}.numCells, numNetworks);
            end
        end

    else

        % otherwise set up the initial conditions of the network from the cell groups
        numNetworks = 1;
        for i = 1:numel(cg)
            u0{i}     = cg{i}.potential;
            y0{i}     = cg{i}.spikedRecently;                % TODO should be binary?
            y_ave0{i} = cg{i}.spikeRate.instant;
            if isfield(cg{i}, 'delayLine')
                delayLine{i} = cg{i}.delayLine;
            end
        end
    end

    % threshold and weight change rules are performed on all iterations from all networks
    numTotalIterations = numNetworks * numIterations;


    % ===============  INITIALIZE INNER LOOP OPTIMIZATION VARIABLES  ===========

    % variables initialized here:
    %   isActive(numCells)          = include cell group in the simulation?
    %   activeCellGroupIds(:)       = all active cell group IDs (for iteration)
    %   bw{numCG}(:)                = blockWeight for each input block
    %   thresh{numCG}(numCells)     = spike threshold to be applied during simulation
    %   potDecay(numCG)             = membrane potential decay
    %   y_history{numCG}(:,:,:)     = record of all output spikes
    %   delayLine{numCG}.buffer,idx = delay line buffer for delay line cell groups (if any)
    
    % initialize spike simulation variables
    y_history{numel(cg)} = [];
    delayLine{numel(cg)} = [];
    for i = 1:numel(cg)
 
        isActive(i) = ~cg{i}.isExternal && ~strcmp(cg{i}.cellType, 'disabled');

        % ave_eta is used to calculate the moving average spike rate from the spike train
        ave_eta(i) = 1 - exp(- 1 / cg{i}.spikeRate.instantWindow);

        if ~isActive(i)
            continue;     % the variables below are not needed for inactive groups
        end
        
        % preinitialize history output variable(s)
        y_history{i} = zeros(cg{i}.numCells, numNetworks, numIterations, model.precisionHint);
        if model.stats.keepPotentialHistory
            % save the potential history in single precision to save on space
            model.snapshot.potentialHistory{i} ...
                    = zeros(cg{i}.numCells, numNetworks, numIterations, 'single');
        end

        % replicate spikeThresh across networks as a hack to speed up thresholding
        thresh{i} = cg{i}.spikeThresh * ones(1,numNetworks);

        % determine which input blocks are enabled and their effective block weights (bw)
        bw{i} = [];
        for b = 1:numel(cg{i}.inputBlock)
            bw{i}(b) = cg{i}.inputBlock(b).blockWeight;
            switch cg{i}.inputBlock(b).connectionType
            case 'excitatory'
                bw{i}(b) =  bw{i}(b) * model.spikeSize;
            case 'inhibitory'
                bw{i}(b) = -bw{i}(b) * model.spikeSize;
            case 'continuous'
                bw{i}(b) =  bw{i}(b);
            case 'disabled'
                bw{i}(b) = 0;
            otherwise
                error('unsupported connectionType');
            end
        end
        enabledBlocks{i} = find(bw{i} ~= 0);

        % preallocate any delay line buffers
        if isfield(cg{i}, 'delayLine')
            delayLine{i}.buffer = zeros(cg{i}.numCells, numNetworks, cg{i}.delayLine.len);
            delayLine{i}.idx    = 1;
        end
    end
    activeCellGroupIds = find(isActive);

    % precalculate any fixed inputs from external sources
    for i = activeCellGroupIds
        u_fixedInput{i} = 0;
        for b = enabledBlocks{i}
            k = cg{i}.inputBlock(b).srcId;
            if cg{k}.isExternal && ~fullInputProvided(k)
                W = cg{i}.inputBlock(b).weight;
                % DRAFT INDIRECT: ... y0{k}(cg{i}.inputBlock(b).inputIndices,:));
                u_fixedInput{i} = u_fixedInput{i} + bw{i}(b) * W * y0{k};
                bw{i}(b) = 0;
            end
        end
        enabledBlocks_nonFixed{i} = find(bw{i} ~= 0);
    end
 

    % =======================  SIMULATE THE SPIKING NETWORK(S)  ======================

    % this is the core simulation engine and consumes most of the CPU time.
    % the algorithm works as follows:
    %
    % for each simulation time step, and for all network clones in parallel:
    %    - propagate spikes from outputs to inputs
    %    - copy any externally provided inputs into simulation inputs
    %    - for each cell group:
    %        - for each input block in that cell group:
    %            - process input spikes through connection weights and
    %                  update membrane potentials of all cells in the cell group
    %                  in all networks
    %        - record any membrane potentials that crossed thresholds as spikes
    %        - keep a record of the spike history across all networks and iterations

    % initialize cell state variables
    u     = u0;
    y     = y0;
    y_ave = y_ave0;                   % only used in moving-average-input mode

    % simulate several iterations of the spiking network
    for t = 1:numIterations

        % set y_prev to the previous iteration's output
        y_prev = y;

        % special case: set y_prev to a computed moving average input instead
        if model.simMovingAveInput
            for i = 1:numel(cg)
                y_ave{i}  = (1-ave_eta(i)) * y_ave{i} + ave_eta(i) * y{i};
            end
            y_prev = y_ave;
        end

        % add time-series input to y_prev (time-series input case only)
        for i = fullInputCellGroupIds
            y_prev{i} = initialState.y{i}(:,:,t);
        end

        for i = activeCellGroupIds
            cg_i = cg{i};

            % special case: delay line - shift values over by moving idx pointer
            if ~isempty(delayLine{i})
                y{i} = delayLine{i}.buffer(:,:,delayLine{i}.idx);
                delayLine{i}.buffer(:,:,delayLine{i}.idx) = y_prev{cg_i.delayLine.srcId};
                delayLine{i}.idx = mod(delayLine{i}.idx, cg_i.delayLine.len) + 1;
                continue;
            end

            % compute the new potential and spike output of each cell
            u_i = exp(-model.simTimeStep / cg_i.membraneTC) * u{i} + u_fixedInput{i};
            % if i == 3; u_i = u_i + .2 * rand(size(u_i)) - .1; end % TODO add jitter to V1i cells
            for b = enabledBlocks_nonFixed{i} 
                k  = cg_i.inputBlock(b).srcId;
                W  = cg_i.inputBlock(b).weight;
                if any(y_prev{k}(:))
                    % DRAFT INDIRECT: ... y_prev{k}(cg_i.inputBlock(b).inputIndices,:));
                    u_i = u_i + bw{i}(b) * W * y_prev{k};
                end
            end
            
            % TODO add gaussian noise to potential to simulate extra-network noise
            %u_i = u_i + model.noiseLevelPot * randn(size(u_i));

            y{i}      = (u_i >= thresh{i});           % spike output if potential exceeds threshold
            u_i(y{i}) = 0;                            % reset potential of spiked cells
            u{i}      = u_i;
            y_history{i}(:,:,t) = y{i};
        end
        
        % special case (debugging): store historical membrane potentials
        if model.stats.keepPotentialHistory
            for i = activeCellGroupIds
                model.snapshot.potentialHistory{i}(:,:,t) = u{i};
            end
        end
    end


    % =======================  UPDATE SPIKING THRESHOLDS  ========================

    % update moving average mean spike rate measurements
    % (note: cg.spikeRate.instant is not updated here)
    thisSpikeRate = [];

    % only compute spike rate for active cell groups
    %for i = activeCellGroupIds
    %    thisSpikeRate{i} = mean(reshape(y_history{i},cg{i}.numCells,[]), 2) / model.simTimeStep;

    % compute spike rate for all cell groups, including input
    for i = 1:numel(cg)
        if strcmp(cg{i}.cellType, 'disabled')
            continue;
        elseif fullInputProvided(i)
            thisSpikeRate{i} = mean(reshape(initialState.y{i}, cg{i}.numCells,[]), 2) / model.simTimeStep;
        elseif cg{i}.isExternal
            thisSpikeRate{i} = mean(y0{i}, 2) / model.simTimeStep;
        else
            thisSpikeRate{i} = mean(reshape(y_history{i}, cg{i}.numCells,[]), 2) / model.simTimeStep;
        end

        sr                = cg{i}.spikeRate;
        stats_eta         = 1 - exp(- numTotalIterations / sr.meanWindow);
        sr.meanBiased     = (1-stats_eta) .* sr.meanBiased + stats_eta * thisSpikeRate{i};
        sr.meanBiased_max = (1-stats_eta) .* sr.meanBiased_max + stats_eta;
        sr.mean           = sr.meanBiased ./ sr.meanBiased_max;
        sr.popMean        = mean(sr.mean);
        cg{i}.spikeRate   = sr;
    end

    % update spike thresholds
    for i = activeCellGroupIds
        % update thresholds with Foldiak's rule: keep each neuron firing near target
        % from SAILnet: delta_theta_i = gamma * (n_i - p)
        numTotalTimeUnits  = numTotalIterations * model.simTimeStep;
        lrate              = model.learningRate * model.lrateScale * cg{i}.threshAdaptRate;
        delta_spikeThresh  = lrate * (thisSpikeRate{i} - cg{i}.targetSpikeRate) * numTotalTimeUnits;
        %p = .5;    % TODO alternate experimental adaptation rule
        %delta_spikeThresh  = lrate * (thisSpikeRate{i} - cg{i}.targetSpikeRate*cg{i}.spikeThresh.^p) * numTotalTimeUnits;
        cg{i}.dspikeThresh = cg{i}.dspikeThresh + delta_spikeThresh;

        % only update spike threshold (spikeThresh) occasionally
        cg{i}.updateThreshWait = cg{i}.updateThreshWait + numTotalIterations;
        if cg{i}.updateThreshEvery > 0 ...
                && cg{i}.updateThreshWait >= cg{i}.updateThreshEvery
            cg{i}.spikeThresh      = cg{i}.spikeThresh + cg{i}.dspikeThresh;
            cg{i}.dspikeThresh(:)  = 0;
            cg{i}.updateThreshWait = 0;
        end
    end


    % =======================  UPDATE CONNECTION WEIGHTS  ========================

    % variables initialized here:
    %   y_ave{numCG}(numCells)        = moving average spike rate (spikes/iteration)
    %                                       - ignored in meanRateLearning mode
    %   inputSpikeRate{numCG}(:,:)    = input spike rates (spikes/iteration, pre-spike)
    %   outputSpikeRate{numCG}(:,:)   = output spike rates (spikes/iteration, post-spike)
    %   numSamples                    = number of columns in spikeRate matrices
    %                                   (either numNetworks or numNetworks*numIterations)

    y_ave = y_ave0;                      % note: not updated in meanRateLearning mode
 
    % prepare the input and output spike histories to use for the learning rules
    if model.meanRateLearning

        % collapse spike history across iterations into a single mean value
        for i = 1:numel(cg)
            if isActive(i)
                y_mean = mean(y_history{i}, 3);
            elseif fullInputProvided(i)
                y_mean = mean(initialState.y{i}, 3);
            else
                y_mean = y0{i};
            end
            inputSpikeRate{i}  = y_mean;
            outputSpikeRate{i} = y_mean;
        end
        numSamples = numNetworks;
    else

        % compute moving average history from cell outputs
        % y_ave_history{numCG}(:,:,:) = moving average of output spikes, post-spike
        y_ave_history{numel(cg)}  = [];
        for i = 1:numel(cg)
            if strcmp(cg{i}.cellType, 'disabled')
                continue;
            end
            if cg{i}.isExternal && ~fullInputProvided(i)
                % special case for static input
                if strcmp(model.stdpEnvelopeShape, 'expDecay_continuous_all')
                    % calculate moving average from static input
                    y_ave_history{i} = zeros(cg{i}.numCells, numNetworks, numIterations, model.precisionHint);
                    for t = 1:numIterations
                        y_ave{i} = (1-ave_eta(i))*y_ave{i} + ave_eta(i) * y0{i};
                        y_ave_history{i}(:,:,t) = y_ave{i};
                    end
                else
                    % otherwise, create history by replicating initial value
                    y_ave_history{i}  = repmat(y0{i}, [1,1,numIterations]);
                end
            else
                if fullInputProvided(i)
                    y_history{i} = initialState.y{i};          
                end
                if strncmp(model.stdpEnvelopeShape, 'expDecay_continuous', 19)
                    y_ave_history{i} = zeros(cg{i}.numCells, numNetworks, numIterations, model.precisionHint);
                    for t = 1:numIterations
                        y_ave{i} = (1-ave_eta(i))*y_ave{i} + ave_eta(i) * y_history{i}(:,:,t);
                        y_ave_history{i}(:,:,t) = y_ave{i};
                    end
                else
                    % generate an averaging envelope kernel to compute weighted average
                    % (note: y_ave0 is ignored with averaging envelope calculations)
                    [offset1, offset2] = meshgrid(1:numIterations, 1:numIterations);
                    switch model.stdpEnvelopeShape
                    case 'expDecay'
                        % exponential decay, post-spike only
                        kernel = exppdf(offset1-offset2, cg{i}.spikeRate.instantWindow);
                    case 'expSymmetric'
                        % exponential decay, symmetric pre-post spike
                        kernel = exppdf(abs(offset1-offset2), cg{i}.spikeRate.instantWindow);
                    case 'gaussian'
                        % gaussian averaging window
                        kernel = normpdf(offset1, offset2, cg{i}.spikeRate.instantWindow);
                    end
                    kernel = kernel ./ mean(sum(kernel,1));   % normalize for consistency across kernels
                    y_ave_history{i} = reshape( reshape(y_history{i}, [], numIterations) * kernel, size(y_history{i}));
                end
            end

            % rehape into a 2D matrix that is one large collection of time samples
            outputSpikeRate{i}  = reshape(y_ave_history{i}, [], numTotalIterations);
            if model.stdpTimeOffset
                % for input spike train, use spikes shifted one time step earlier
                y_ave_history0{i} = cat(3, y_ave0{i}, y_ave_history{i}(:,:,1:end-1));
                inputSpikeRate{i} = reshape(y_ave_history0{i}, [], numTotalIterations);
            else
                inputSpikeRate{i} = outputSpikeRate{i};
            end
        end
        numSamples = numTotalIterations;
    end

    % update spike rates used by learning rules (spike rates are in spikes/timeUnit)
    % (TODO use spikes/iteration instead, to avoid need for later conversion?)
    for i = 1:numel(cg)
        if ~ strcmp(cg{i}.cellType, 'disabled')
            sr                 = cg{i}.spikeRate;
            stats_eta          = 1 - exp(- numTotalIterations / sr.lMeanWindow);
            thisSpikeRate      = mean(outputSpikeRate{i}, 2) / model.simTimeStep;
            sr.lMeanBiased     = (1-stats_eta) .* sr.lMeanBiased + stats_eta * thisSpikeRate;
            sr.lMeanBiased_max = (1-stats_eta) .* sr.lMeanBiased_max + stats_eta;
            sr.lMean           = sr.lMeanBiased ./ sr.lMeanBiased_max;
            %thisSpikeRate2     = mean(outputSpikeRate{i}.^2, 2) / model.simTimeStep^2;   % track n^2
            %sr.lMean2Biased    = (1-stats_eta) .* sr.lMean2Biased + stats_eta * thisSpikeRate2;
            %sr.lMean2          = sr.lMean2Biased ./ sr.lMeanBiased_max;
            cg{i}.spikeRate    = sr;
        end
    end

    % GetXYSpiked() - subfunction for generating spike history matrices for STDP rules
    y_history0{numel(cg)} = [];
    function [X_spiked, Y_spiked] = GetXYSpiked(k, i)
        if model.meanRateLearning
            error('spike-timing-dependent learning rules are not supported in mean-rate mode');
        end
        if isempty(y_history0{k})
            y_history0{k} = cat(3, y0{k}, y_history{k}(:,:,1:end-1));
        end
        X_spiked = reshape(y_history0{k}, [], numTotalIterations);
        Y_spiked = reshape(y_history{i},  [], numTotalIterations);
    end

    % calculate weight changes (dweight) using learning rules
    for i = activeCellGroupIds
        for b = enabledBlocks{i}
            k = cg{i}.inputBlock(b).srcId;
            W = cg{i}.inputBlock(b).weight;
            X = inputSpikeRate{k};
            Y = outputSpikeRate{i};

            if cg{i}.inputBlock(b).learningRate == 0 || cg{i}.learningRate == 0 || model.learningRate == 0
                continue;          % learning for this input block is disabled
            end

            % DRAFT INDIRECT: code to make learning rule work with inputIndices
            % (need to iterate over output cells!)
            % if isempty(cg{i}.inputBlock(b).inputIndices)
            %    outputIndices = 0;
            %    X = X_all;
            % else
            %    outputIndices = 1:cg{i}.numCells;
            % end
            % for j = outputIndices
            %    if j > 0
            %        X = X_all(cg{i}.inputBlock(b).inputIndices(j,:), :);
            %    end
            %    switch ...
            %    if j > 0
            %        cg{i}.inputBlock(b).dweight(j,:) = cg{i}.inputBlock(b).dweight(j,:) + dW;
            %    else
            %        cg{i}.inputBlock(b).dweight = cg{i}.inputBlock(b).dweight + dW;
            %    end
            % end

            switch cg{i}.inputBlock(b).learningRule

            case 'hebbian_oja'
                % Oja variant of Hebbian learning rule (dW = y x - y^2 W)
                % (attenuates weight growth and learns linear models)
                dW = Y * X' - diag(sum(Y.^2,2)) * W;

            case 'foldiak'
                % Foldiak learning rule to strengthen correlated connections (dW = y x - <y> <x>)
                % (should use measured spike rates, but SAILnet uses target spike rate)
                %pq = (cg{i}.spikeRate.lMean * cg{k}.spikeRate.lMean') * model.simTimeStep^2;
                p  = cg{i}.targetSpikeRate * model.simTimeStep;
                q  = cg{k}.targetSpikeRate * model.simTimeStep;
                dW = Y * X' - p * q * numSamples;

            case 'correlation_measuring'
                % Foldiak-inspired rule that converges despite correlations (dW = y x - <y> <x> (1 + W))
                % pq is the moving average of the input and output cell spike rates
                % (adjusted to be per-iteration rather than per-time-unit)
                pq = (cg{i}.spikeRate.lMean * cg{k}.spikeRate.lMean') * model.simTimeStep^2;
                dW = (Y * X') - pq .* numSamples .* (1 + W);

            % --------------  experimental learning rules below here  -----------------
                
            case 'hebbian'
                % vanilla hebbian rule (this rule grows forever and requires weight rescaling)
                dW = Y * X';

            case 'anticorrelation'
                % rule that makes weight stronger in proportion to anticorrelation instead of correlation
                pq = (cg{i}.spikeRate.lMean * cg{k}.spikeRate.lMean') * model.simTimeStep^2;
                %dW = (pq  .* exp(-W)) * numSamples - (Y * X');       % exponential version
                dW = pq .* numSamples - (Y * X') .* (1 + W);

            case 'foldiak_bounded_exp'
                % Foldiak-based rule that converges toward equilibrium despite persistent correlation
                % note: the use of a dynamically measured spikeRate means that changing
                % the training batch size will slightly alter network performance.
                %p  = cg{i}.targetSpikeRate * model.simTimeStep;
                %q  = cg{k}.targetSpikeRate * model.simTimeStep;
                p  = mean(cg{i}.spikeRate.lMean) * model.simTimeStep;
                q  = mean(cg{k}.spikeRate.lMean) * model.simTimeStep;
                %p  = mean(Y(:));
                %q  = mean(X(:));
                dW = (Y * X') .* exp(-W) - p * q * numSamples;           % original FB rule

            case 'stdp'
                % spike-timing dependent plasticity (STDP)
                % dW = alpha+ * X * Y_spiked - alpha- * Y * X_spiked
                [X_spiked, Y_spiked] = GetXYSpiked(k, i);
                dW = Y_spiked * X' - Y * X_spiked';

            case 'stdp_attenuated'
                % spike-timing dependent plasticity (STDP) with Maass et al (2010) attenuation
                % dW = alpha+ * X * Y_spiked - alpha- * Y * X_spiked
                [X_spiked, Y_spiked] = GetXYSpiked(k, i);
                dW = (Y_spiked * X' - Y * X_spiked') .* exp(-abs(W));

            case 'none'
                % learning turned off for this input block

            otherwise
                error('learning rule not supported.');

            end

            % store updated dweight back into cellgroup
            cg{i}.inputBlock(b).dweight = cg{i}.inputBlock(b).dweight + dW;
        end
    end

    % calculate learning rate adjustment to convert spikes/iteration to spikes/timeUnit
    % simTimeStep compensates for fact that learning rule is applied k times per time unit
    % (1/simTimeStep)^2 scales spikeRate^2 from (spikes/iteration)^2 to (spikes/timeUnit)^2
    lrateScaleW = model.lrateScale * model.simTimeStep * (1/model.simTimeStep)^2;

    % update weights from dweight every so often
    for i = activeCellGroupIds
        cg{i}.updateWait = cg{i}.updateWait + numTotalIterations;
        if cg{i}.updateEvery > 0 && cg{i}.updateWait >= cg{i}.updateEvery
            lrateScale = model.learningRate * lrateScaleW * cg{i}.learningRate;
            for b = enabledBlocks{i}
                lrate  =  lrateScale * cg{i}.inputBlock(b).learningRate;
                if model.meanRateLearning
                    lrate = lrate * numIterations;      % compensate for collapsed iterations
                end
                if lrate == 0
                    continue;              % learning for this input block is disabled
                end
                weight = cg{i}.inputBlock(b).weight + lrate * cg{i}.inputBlock(b).dweight;
                cg{i}.inputBlock(b).dweight(:) = 0;

                % clamp any requested weights to zero (e.g. self-weights)
                weight(cg{i}.inputBlock(b).clampZero) = 0;

                % constrain weights to be nonnegative
                if strcmp(cg{i}.inputBlock(b).constrainWeights, 'nonneg')
                    weight(weight < 0) = 0;
                end

                cg{i}.inputBlock(b).weight = weight;
            end
            cg{i}.updateWait = 0;
        end
    end


    % ====================  RETURN VALUES   ==================

    % save the final network state values back into cell groups
    % (only if we're not in "initial state" batch mode)
    if isempty(initialState)
        for i = activeCellGroupIds
            cg{i}.potential          = u{i};
            cg{i}.spikedRecently     = y{i};       % TODO partial implementation of "spikedRecently"
            cg{i}.spikeRate.instant  = y_ave{i};
            if ~isempty(delayLine{i})
                cg{i}.delayLine = delayLine{i};
            end
        end
    end

    % save modified cell groups back into model
    model.cellGroup = cg;

    % update performance tracking statistics and figures
    model = NetModel_Stats(model, y_history);

    % return spike history if output argument was provided
    if nargout >= 2
        finalState = struct();
        finalState.y         = y;
        finalState.y_ave     = y_ave;
        finalState.u         = u;
        finalState.y_history = y_history;
        finalState.delayLine = delayLine;
    end
end

