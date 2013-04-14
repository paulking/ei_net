% CellGroup_AddInput - Add a new input block (connection) to a cell group
%
% Add a new input block to a cell group.
%
% Specific connections between cells in the source and target cell groups
% are generated automatically by sampling from the source cell group randomly.
% The random sampling method can be selected with the optional 'connectionPattern'
% parameter.
%
% A side-effect of the random sampling of connections is that connections
% may be duplicated (the same target cell may be assigned multiple
% inputs from the same source cell). While this is statistically reasonable,
% it is computationally redundant. If 'noDuplicates' is true, random connection
% samples that duplicate an existing connection will be resampled until they 
% are unique.  While this succeeds in removing duplicates, it disrupts
% the statistical distribution of the sampling method.  So, for example,
% using the 'geolocal' gaussian connection pattern with noDuplicates will
% result in a non-gaussian distribution if a lot of resampling is done.
%
% Usage:
%    cellGroup = CellGroup_AddInput(cellGroup, params)
%
% Inputs:
%    params    = the parameters defining the input block configuration
%        sourceCellGroup   = the cell group providing the inputs (required)
%        numInputs         = number of inputs per cell (required)
%        connectionPattern = the connection pattern to use. (default = 'full')
%            'full'            = full connectivity.  If numInputs is supplied,
%                                it must equal sourceCellGroup.numCells.
%                                No indirection matrix (inputIndices) is used
%                                for greater speed.
%            'uniform'         = each target cell is connected randomly to a source
%            'fullIndirect'    = full connectivity, but with the usual
%                                indirection matrix (really just for testing).
%            'geolocal'        = connect probabilistically to the closest cells
%                                according to cell location.
%            'geolocalClosest' = connect to the closest cells according to cell location.
%            'none'            = uninitialized; the caller will wire later
%        connectionSigma   = when using the 'geolocal' method, the standard deviation
%                            to use when determining connection probability.
%                            If null is passed in, a value will be calculated
%                            based on the ratio of inputs to source cells.
%        noDuplicates      = if true, duplicate connections (connections between
%                            the same input/output cell pairs) will be thrown out
%                            and resampled. (default = true)
%        noSelfConnections = if true, this input block connects a cell group to itself,
%                            and connections where the input and output cell are the same
%                            will be disallowed. Disallowing self-connections is important
%                            for network stability to prevent self-reinforcing feedback loops
%                            and infinitely scaling connection weights. For fully-connected
%                            input blocks, self-connections will be suppressed by clamping
%                            the diagonal weights, weight(i,i), to zero. (default = false)
%        forceDirect       = force all connections to be direct, with missing connections
%                            marked with connection clamping. This is the fastest
%                            method for execution as long as the input reduction
%                            is not too great (up to 10:1 reduction?).
%        connectionType    = type of connection, if it is desired to override
%                            the default type determined from the input cell type.
%                            (default = input cell type)
%        learningRule      = learning rule to use (default = derive from source cell group).
%                            Partial list:
%            'hebbian_oja'
%            'foldiak'
%            'correlation_measuring'
%            'stdp'
%        learningRate      = learning are to apply. Inhibitory connections should
%                            use a negative learning rate. This is later multiplied
%                            by cellGroup.learningRate and model.learningRate
%                            to arrive at a final value.  (default = +/- .01)
%        blockWeight       = how much to weight this block of inputs relative
%                            to the other blocks. (default = 1)
%        initWeights       = method to use for initializing weights (default = 'uniform')
%                            weight blocks are normalized to unit vector length input
%                            to each cell (unless the weights are zero). Unconstrained weights
%                            will have a zero mean. 'nonneg' weights will all be >= 0.
%            'uniform'         = uniform random values (default)
%            'gaussian'        = gausian distribution, or positive half-gaussian
%                                if constrainWeights = 'nonneg'.
%            'zero'            = initialize all weights to zero
%            scalar value      = weights are initialized to that value
%            matrix            = weight are initialized to the supplied weight matrix
%        constrainWeights  = constraint to apply to weights (default = 'nonneg').
%            'nonneg'          = weights cannot go below zero (default)
%            'none'            = no constraint to apply
%        clampZero         = array(n,1) or logical array(numCells,numInputs) indicating
%                            which weights, if any, to clamp to zero. (default = []).
%                            Clamping weights to zero can be used to simulate
%                            partial connectivity. When a cell group is connected
%                            to itself, all self-weights (weight(i,i) along the diagonal)
%                            should be clamped to zero.
%
% Output:
%    cellGroup  = the cell group with the new input block added
%
% See also:
%    CellGroup_Init, CellGroup_Update
%
% ISSUES
%    - should params.sourceCellGroup actually be the cell group id
%        - would require passing in a model, which might be a little circular!
%    - save connection pattern name in inputBlock?
%        - in case there is a desire to use the pattern to allow connection rewiring
%        - i.e. structural plasticity
%
% Created:   11/15/09, Paul King
%--------------------------------------------------------------------------
function cellGroup = CellGroup_AddInput( cellGroup, params )

    % initialize and apply default parameter values
    if nargin < 2
        params = struct();
    end
    defaultValues.sourceCellGroup   = [];             % must be provided
    defaultValues.numInputs         = [];             % default computed (below)
    defaultValues.connectionPattern = 'full';
    defaultValues.connectionSigma   = [];             % default computed (below)
    defaultValues.noDuplicates      = true;
    defaultValues.noSelfConnections = false;
    defaultValues.forceDirect       = cellGroup.defForceDirect;
    defaultValues.connectionType    = [];             % default computed (below)
    defaultValues.learningRule      = [];             % default computed (below)
    defaultValues.learningRate      = [];             % default computed (below)
    defaultValues.blockWeight       = 1;
    defaultValues.initWeights       = 'uniform';
    defaultValues.initWeightsNorm   = 'unitLength';
    defaultValues.constrainWeights  = 'nonneg';
    params = ApplyDefaultValues(params, defaultValues);
    
    assert( ~isempty(params.sourceCellGroup) );       % sourceCellGroup is a required parameter

    % calculate computed defaults
    if isempty(params.numInputs) && any(strcmp(params.connectionPattern, {'full','fullIndirect'}))
        params.numInputs = params.sourceCellGroup.numCells;
    end
    if isempty(params.connectionType)
        params.connectionType = params.sourceCellGroup.cellType;
    end
    if isempty(params.learningRate)
        if strcmpi(params.connectionType, 'inhibitory')
            params.learningRate = -1;
        else
            params.learningRate = 1;
        end
    end
    if isempty(params.learningRule)
        params.learningRule = params.sourceCellGroup.defLearningRule;
    end

    assert(any( strcmp(params.constrainWeights, {'nonneg','none'}) ));

    if params.noSelfConnections && params.sourceCellGroup.numCells ~= cellGroup.numCells
        error('noSelfConnections is true, but input and output cell populations differ.');
    end


    % ======================   INITIALIZE STRUCT   =======================

    % the number of inputs for each element of the input block
    numInputs      = params.numInputs;
    numSourceCells = params.sourceCellGroup.numCells;

    % initialize the input block data structure
    inputBlock = struct();
    inputBlock.connectionType    = params.connectionType;
    inputBlock.numSourceInputs   = numSourceCells;
    inputBlock.name              = params.sourceCellGroup.name;
    inputBlock.srcId             = 0;                            % filled in later
    inputBlock.numInputs         = numInputs;
    inputBlock.blockWeight       = params.blockWeight;
    inputBlock.learningRule      = params.learningRule;
    inputBlock.learningRate      = params.learningRate;
    inputBlock.constrainWeights  = params.constrainWeights;
    inputBlock.clampZero         = [];
    inputBlock.inputIndices      = [];
    inputBlock.spikedRecently    = zeros(numSourceCells, 1);
    inputBlock.spikeRate         = struct();                     % filled in later
    inputBlock.causalFeedback    = [];
    inputBlock.weight            = [];
    inputBlock.dweight           = [];                           % accumulated delta to weight

    % initialize the causal feedback state if appropriate
    if ~ isempty(params.sourceCellGroup.causalFeedback)
        inputBlock.causalFeedback = zeros(numSourceCells, 1);
    end


    % ======================   ADD CONNECTIONS   =======================

    % add the wiring connections
    switch params.connectionPattern
    case 'geolocal'
        % wire connections according to a geometric local connection probability
        inputBlock = WireConnections_Geolocal(inputBlock, cellGroup, params);

    case 'geolocalClosest'
        % wire connections according to a geometric local connection probability
        inputBlock = WireConnections_GeolocalClosest(inputBlock, cellGroup, params);

    case 'uniform'
        % wire connections using uniform connectivity
        inputBlock = WireConnections_Uniform(inputBlock, cellGroup, params);

    case 'full'
        % use direct full connectivity with no inputIndices indirection matrix
        if numInputs ~= numSourceCells
            error('full connectivity requires that numInputs equal the number of source cells.');
        end
        inputBlock.inputIndices = [];

    case 'fullIndirect'
        % wire every source cell to every target cell
        if numInputs ~= numSourceCells
            error('full connectivity requires that numInputs equal the number of source cells.');
        end
        inputBlock.inputIndices = ones(cellGroup.numCells,1) * (1:params.sourceCellGroup.numCells);

    case 'none'
        % leave the indices uninitialized; The caller must initialize later
        inputBlock.inputIndices = zeros(cellGroup.numCells, numInputs);

    otherwise
        error('unknown connection pattern');
    end
    inputBlock.inputIndices = uint16(inputBlock.inputIndices);    % speeds up later access

    % convert connections to direct connectivity with clamping
    % (this is faster for execution as long as reduction ratio is not too high... 10:1 max?)
    if params.forceDirect && ~isempty(inputBlock.inputIndices)
        if ~isempty(inputBlock.clampZero)
            error('forceDirect is not compatible with clamping input parameter!');
        end
        if numInputs < .1 * numSourceCells
            fprintf('WARNING: forceDirect used with large connection reduction ratio.\n');
        end
        % create new logical array and combine connections with pre-existing clamping
        inputBlock.clampZero = true(cellGroup.numCells, numSourceCells);
        for i = 1:cellGroup.numCells
            % open holes in the connection-supression matrix for each active connection
            inputBlock.clampZero(i, inputBlock.inputIndices(i,:)) = false;
        end
        inputBlock.inputIndices = [];
        inputBlock.numInputs    = numSourceCells;
        numInputs               = numSourceCells;
    end

    % prevent self-connections by clamping any self-weights to zero
    if params.noSelfConnections
        if ~isempty(inputBlock.inputIndices) && (numInputs ~= numSourceCells)
            error('can''t disable self-connection for this connectionPattern');
        end
        N = numSourceCells;
        selfWeightIdx = (0:N-1)*(N + 1) + 1;
        if isempty(inputBlock.clampZero)
            inputBlock.clampZero = selfWeightIdx;
        elseif islogical(inputBlock.clampZero)
            inputBlock.clampZero(selfWeightIdx) = true;
        else isnumeric(inputBlock.clampZero)
            inputBlock.clampZero = union(inputBlock.clampZero, selfWeightIdx);
        end
    end


    % ===================   INITIALIZE WEIGHTS   ====================

    weightDims           = [cellGroup.numCells, numInputs];
    constrainNonNegative = strcmp(params.constrainWeights, 'nonneg');
    
    if isnumeric(params.initWeights)
        % initialize weights to a numeric value, or possibly an explicit matrix
        if numel(params.initWeights) == 1
            W = params.initWeights * ones(weightDims);
        elseif all(size(params.initWeights) == weightDims)
            W = params.initWeights;
        else
            error('invalid dimensions for initWeights matrix');
        end
    else
        % initialize weights according to a named initialization scheme
        switch params.initWeights
        case 'uniform'
            if constrainNonNegative
                W = rand(weightDims);
            else
                W = rand(weightDims)*2 - 1;
            end
        case 'gaussian'
            if constrainNonNegative
                W = abs(randn(weightDims));
            else
                W = randn(weightDims);
            end
        case 'zero'
            W = zeros(weightDims);
        otherwise
            error('unrecognized weight initialization scheme');
        end
    end
    W(inputBlock.clampZero) = 0;
    if strcmp(params.initWeightsNorm, 'unitLength') && ~strcmp(params.initWeights, 'zero')
        % optionally normalize weights so the weight vector input to each cell has unit length
        rescale = 1 ./ sqrt(sum(W.^2,2));
        rescale(~isfinite(rescale)) = 1;       % just in case any weights are all zero
        W = bsxfun(@times, W, rescale);
    end
    inputBlock.weight  = W;
    inputBlock.dweight = zeros(weightDims);


    % ==========================   WRAP UP   ===========================

    % add the input block to the cell group
    if isempty(cellGroup.inputBlock)
        cellGroup.inputBlock = inputBlock;
    else
        cellGroup.inputBlock(end+1) = inputBlock;
    end
end


% WireConnections_Geolocal
%
% Wire connections according to euclidean distance of the cell locations.
% Connections are established according to a 2D gaussian probability
% defined by the parameter connectionSigma. If connectionSigma == [],
% then one will be calculated based on the number of inputs relative
% to the number of source inputs.
function inputBlock = WireConnections_Geolocal(inputBlock, cellGroup, params )

    % source and target cell groups must have location information
    % available in order to use this method.
    if isempty(cellGroup.location) || isempty(params.sourceCellGroup.location)
        error(['both source and target cell groups must have locations in order\n' ...
                'to use ''geolocal'' connectivity.' ]);
    end

    % collect the key metrics used for wiring the cells
    numInputs       = params.numInputs;
    numSourceCells  = params.sourceCellGroup.numCells;
    connectionSigma = params.connectionSigma;

    % if no connectionSigma is provided, calculate a reasonable one
    % (sigma such that approximately numInputs will be within two sigma distances)
    if isempty(connectionSigma)
        connectionSigma  = (1/2) * sqrt(numInputs / numSourceCells);
    end
    
    % the number of times to try to get a valid connection before giving up
    maxAttempts  = 100;
    totalSamples = 0;

    % assign a set of connections for each cell in the cell group
    inputBlock.inputIndices = zeros(cellGroup.numCells, numInputs);
    for i = 1:cellGroup.numCells

        % calculate the connection probability density to all source cell locations
        % (connection probability is based on euclidean distance)
        delta1       = params.sourceCellGroup.location(:,1) - cellGroup.location(i,1);
        delta2       = params.sourceCellGroup.location(:,2) - cellGroup.location(i,2);
        dist_squared = delta1.^2 + delta2.^2;
        pd           = 1/(2*pi*connectionSigma^2) * exp(-.5 * dist_squared / connectionSigma^2);
        %pd          = mvnpdf(sqrt(dist_squared)*[1 0], 0, connectionSigma^2*[1 1]);

        % if we need to avoid duplicates, create a master list of available cells
        availableSourceCells = [];
        if params.noDuplicates
            availableSourceCells = 1:numSourceCells;
        end

        % generate random input connections for each cell according to
        % a gaussian distribution centered on the receiving neuron's location
        numConnected = 0;
        numAttempts  = 0;
        while numConnected < numInputs && numAttempts < maxAttempts
            numAttempts        = numAttempts + 1;
            numInputsRemaining = numInputs - numConnected;

            % pick a set of candidate input cells
            if params.noDuplicates && ~isempty(availableSourceCells)
                candidateIds = availableSourceCells;
            else
                candidateIds = 1:numSourceCells;
            end
            totalSamples = totalSamples + numel(candidateIds);

            % calculate the probability of connecting to each candidate input cell, given
            % the probability density and the number of additional connections needed
            pd_candidates = pd(candidateIds);
            p_connection  = numInputsRemaining * (pd_candidates / sum(pd_candidates));
            
            % randomly test the connection probabilities to find the cells that succeed
            selected      = rand(size(p_connection)) <= p_connection;
            selectedIds   = candidateIds(selected);

            % disallow duplicate connections
            if params.noDuplicates
                if ~isempty(availableSourceCells)
                    % duplicate check already made, just remove selected cells from list
                    availableSourceCells = availableSourceCells(~selected);
                else
                    selectedIds = setdiff(unique(selectedIds), inputBlock.inputIndices(i,1:numConnected));
                end
            end

            % disallow self-weights (connections from a cell to itself)
            if params.noSelfConnections
                selectedIds = selectedIds(selectedIds ~= i);
            end

            % add the selected cells to the input block
            if numel(selectedIds) > 0
                while numel(selectedIds) > numInputsRemaining
                    % picked too many, delete a random subset
                    numExcess     = numel(selectedIds) - numInputsRemaining;
                    deleteIndices = ceil( rand(numExcess,1) * numel(selectedIds) );
                    keep          = true(size(selectedIds));
                    keep(deleteIndices) = false;
                    selectedIds = selectedIds(keep);
                end
                % this works, but is biased rather than random
                %if numel(selectedIds) > numInputsRemaining
                %    % picked too many, keep the most probable connections
                %    [~, indices] = sort(pd(selectedIds),'descend');
                %    selectedIds  = selectedIds(indices(1:numInputsRemaining));
                %end
                slots = (numConnected + 1) : (numConnected + numel(selectedIds));
                inputBlock.inputIndices(i,slots) = selectedIds;
                numConnected = numConnected + numel(selectedIds);
            end
        end

        % exit with an error if we couldn't find enough valid connections
        if numConnected < numInputs
            if params.noDuplicates
                disp('disallowing duplicates may have caused all connection options to be exhausted.');
            end
            error('exceeded maxAttempts without successfully completing connections');
        end
    end
    
    % DEBUG - to see how efficient the algorithm is
    avgSamplesPerConnection = totalSamples / (cellGroup.numCells * numInputs);
    if avgSamplesPerConnection > 10
        fprintf('warning: high number of samples per connection (%.1f)\n', ...
                avgSamplesPerConnection);
    end
end


% WireConnections_GeolocalClosest
%
% Wire connections according to euclidean distance of the cell locations.
% Rather than connecting probabilistically, as above, the closest connections
% are made. There will never be any duplicates with this scheme, regardless
% of the setting of the noDuplicates field.
%
function inputBlock = WireConnections_GeolocalClosest(inputBlock, cellGroup, params )

    % source and target cell groups must have location information
    % available in order to use this method.
    if isempty(cellGroup.location) || isempty(params.sourceCellGroup.location)
        error(['both source and target cell groups must have locations in order\n' ...
                'to use ''geolocalClosest'' connectivity.' ]);
    end

    if params.numInputs > (params.sourceCellGroup.numCells - params.noSelfConnections)
        error('Too many inputs per cell to satisfy all constraints.');
    end

    % assign a set of connections for each cell in the cell group
    inputBlock.inputIndices = zeros(cellGroup.numCells, params.numInputs);
    for i = 1:cellGroup.numCells

        % calculate distance to all source cells
        delta1  = params.sourceCellGroup.location(:,1) - cellGroup.location(i,1);
        delta2  = params.sourceCellGroup.location(:,2) - cellGroup.location(i,2);
        dist    = sqrt(delta1.^2 + delta2.^2);

        % disallow self-weights (connections from a cell to itself)
        if params.noSelfConnections
            dist(i) == Inf;        % make the connection too far away to connect to
        end

        % sort source cells by distance and pick the N closest
        [distSorted, indicesSorted] = sort(dist);
        inputBlock.inputIndices(i,:) = indicesSorted(1:params.numInputs);
    end
end


% WireConnections_Uniform
%
function inputBlock = WireConnections_Uniform(inputBlock, cellGroup, params )

    % collect the key metrics used for wiring the cells
    numInputs      = params.numInputs;
    numSourceCells = params.sourceCellGroup.numCells;

    % number of times to try making a connection before giving up
    % (only an issue when noDuplicates == true)
    maxAttempts = 20;

    % assign a set of connections for each cell in the cell group
    inputBlock.inputIndices = zeros(cellGroup.numCells, numInputs);
    for i = 1:cellGroup.numCells

        % special case: relatively small nujmber of source inputs
        % (select randomly from a list to select connections.
        if params.noDuplicates && numSourceCells < 5*numInputs
            % build a list of all available source cells
            availableSourceCells = 1:numSourceCells;
            if params.noSelfConnections
                availableSourceCells = [1:(i-1), i+1:numSourceCells];
            end

            % randomly select inputs one at a time
            for j = 1:numInputs
                k_idx = floor(rand() * numel(availableSourceCells)) + 1;
                inputBlock.inputIndices(i,j) = availableSourceCells(k_idx);
                availableSourceCells = availableSourceCells([1:k_idx-1, k_idx+1:end]);
            end
            continue;
        end

        % generate uniform random input connections for each cell
        for j = 1:numInputs

            % generate a random connection source
            for attempts = 1:maxAttempts

                sourceIndex = round(rand() * numSourceCells) + 1;

                % don't allow duplicate connections (if requested)
                if params.noDuplicates && any(inputBlock.inputIndices(i,1:j-1) == sourceIndex)
                    continue;
                end
                
                % disallow self-weights (connections from a cell to itself)
                if params.noSelfConnections && sourceIndex == i
                    continue;
                end
                
                % (note: this loop will never iterate if noDuplicates == false)
                break;
            end

            % add the connection to the input block
            inputBlock.inputIndices(i,j) = sourceIndex;
        end
    end
end

