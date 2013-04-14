% NetModel_Init - Initialize a spiking network simulation model
%
% Initialize a spiking network simulation model.
%
% Notes:
%    self-connections:
%       If a cell group is connected to itself, the input block parameter
%       'noSelfConnections' will be set to true by default. For fully
%       connected networks, the diagonal weights, W(i,i), will be clamped
%       to zero.
%    learning rates:
%        If a cell group is inhibitory, its learning rate should probably
%        be negative unless the cell group connects to itself, in which
%        case the learning rate should be positive.
%
% Usage:
%    model = NetModel_Init(params)
%
% Inputs:
%    params      = the network model configuration parameters:
%        refNumCells       = a reference number of cells for each cell group
%                            The actual number of cells for a cell group can then
%                            be specified as a fractional percentage of this number.
%        fpsMax            = maximum allowable iteration rate (frames per second).
%                            if set, this will slow down the simulation to a certain
%                            maximum speed
%        learningRate      = learning rate to apply to the model as a whole. If 0,
%                            learning is disabled for the whole model. (default = 1).
%        lrateScale        = scaling factor to make learning rates more human-readable
%                            (default = .001).
%        simTimeUnitSize   = size of simulation time unit in seconds (default = .001 = 1 ms)
%        simTimeStep       = size of discrete simulation time step in simualtion time
%                            units. The default interpretation is that each simulation
%                            time unit is 1 ms and the time step is 1. (default = 1)
%        simMovingAveInput = simulate synaptic neurotransmitter accumulation
%        simInitPotentials = how to initialize cell membrane potentials before network
%                             simulation. Either 'zero' or 'random' (default = 'zero')
%        simNoReset        = reset potentials before each training sample? (default = true)
%        spikeSize         = how much does a single spike contribute to a target
%                            cell's potential? (default = 1)
%        meanRateLearning  = use mean-rate learning instead of spike-timing learning
%                            (default = false)
%        stdpEnvelopeShape = shape of weighting envelope to use for STDP moving average calculation
%            'expDecay_continuous' = exponential decay, using y_ave0 to start
%            'expDecay'            = exponential decay (typical moving average)
%            'expDecaySymmetric'   = exponential decay but both forward and backward in time
%            'gaussian'            = symmetric gaussian envelope
%        stdpTimeOffset    = when doing STDP learning (non-mean-rate learning), shift
%                            input moving averages to be one time step earlier in time.
%        precisionHint     = floating point precision to favor, either 'single' or 'double'.
%                            'single' can be up to 25% faster; results are 99.9% the same.
%                            (default = 'double')
%        outputCellGroupName = name of cell group defined as representing the output
%                              of the model as a whole. (default = [])
%        cgDefaults        = parameter struct containing default values to apply
%                            to all cell groups in model (default = [])
%        cg_<name>         = parameters for cell group with name <name>
%            numCells         = number of cells in the cell group
%            cellType         = the type of cell (e.g. 'excitatory' or 'inhibitory')
%            refNumInputs     = the target number of inputs to this cell group.  This
%                               number, if provided, is used to calculate the input block
%                               input count when fractional ratios are provided.
%                               This parameter provides a way to experiment
%                               with different numbers of total inputs without
%                               changing the synapse type ratios. (optional)
%            in_<srcName>     = describes a set of inputs to this cell.
%                               <srcName> is the name of the cell group providing
%                               the input.  The value of the parameter can simply
%                               be a number, in which case it is taken to be
%                               the number of inputs to this cell from the source
%                               cell group.  If the number is in the range (0,1)
%                               then it is assumed to be a percentage of
%                               cg_xxx.refNumInputs.
%                numInputs         = if in_<srcName> is a property struct, then
%                                    numInputs contains the value of the number
%                                    of inputs, identical to the description above.
%                connectionType    = type of connection, if it is desired to override
%                                    the default type determined from the input cell type.
%                                    (default = input cell type)
%                connectionPattern = how to wire the source cells to this group.
%                                    options include: 'geolocal', 'uniform', 'full',
%                                    and 'none'. See CellGroup_AddInput for more details.
%                connectionSigma   = gaussian sigma to use for 'geolocal' connection
%                                    method. (see CellGroup_AddInput for details)
%                learningRule      = learning rule to use (default = derive from source cell group)
%                                    (See CellGroup_Init, CellGroup_Update for details)
%                learningRate      = learning are to apply. This is later multiplied
%                                    by cellGroup.learningRate and model.learningRate
%                                    to arrive at a final value.  (default = +/- .01)
%                [See CellGroup_AddInput() for more in_xxx parameters.]
%        cgDefaults        = parameter struct containing default values to apply
%        displayTimeUnits  = time units to use for interpreting display/print frequencies
%            'iteration'       = iterations
%            'sample'          = learning samples
%            'simTimeUnit'     = simulation time units
%        debugDisplayEvery = number of displayTimeUnits between debug output events.
%                            0 = never. Only used for single-sample linear simulation mode.
%                            (default = 0)
%        stats             = struct of params to configure figures and statistics
%                            gathering. (See NetModel_Stats for details.)
%            measure              = measurements to take
%            figure               = figures to draw
%            printStatsEvery      = how often to print a line of statistics. Units determined
%                                   by displayTimeUnits. 0 means don't print statistics.
%                                   (default = 1000)
%            updateFiguresEvery   = how often to update display figures. Units determined
%                                   by displayTimeUnits. Value of [] will copy value from
%                                   printStatsEvery. 0 means don't show figures. (default = [])
%            keepSpikeHistory     = retain record of spikes from simulation? Useful
%                                   for debugging and charting, can take up a lot of memory.
%                                   History is kept in model.snapshot.spikeHistory.
%                                   (default = false)
%            keepPotentialHistory = retain record of membrane potentials from simulation?
%                                   History is kept in model.snapshot.potentialHistory.
%                                   (default = false)
%
% Output:
%    model       = the initialized network model
%        [all properties above, plus...]
%        cellGroup       = cell array(numCellGroups) of cell group structs
%        snapshot        = struct containing temporarily saved information from
%                          one "batch" of simulations; can be deleted at any time
%                          to free up memory.
%            inputData        = matrix(numInputs,numSamplesPerBatch,numIterationsPerSample)
%                               of inputData.
%            spikeHistory     = cell array(numCellGroups) of matrix(numCells,
%                               numSamplesPerBatch,numIterationsPerSample) containing
%                               retained spike history for later analysis.
%            potentialHistory = cell array(numCellGroups) of matrix(numCells,
%                               numSamplesPerBatch,numIterationsPerSample) containing
%                               retained cell potentials for later analysis.
%
% See also:
%    CellGroup_AddInput, CellGroup_Init, V1NetModel_Init, NetModel_Update, NetModel_Stats
%
% Created:   11/13/09, Paul King
%--------------------------------------------------------------------------
function model = NetModel_Init( params )

    % ===================   SET DEFAULT MODEL PARAMETERS   ===================

    % initialize default parameter values
    if nargin < 1
        params = struct();
    end
    defaultValues.fpsMax                 = 0;           % don't limit frame rate
    defaultValues.learningRate           = 1;
    defaultValues.lrateScale             = .001;
    defaultValues.simTimeUnitSize        = .001;        % default simulation time unit is 1 ms
    defaultValues.simTimeStep            = 1;
    defaultValues.simMovingAveInput      = false;
    defaultValues.simInitPotentials      = 'zero';
    defaultValues.simNoReset             = false;
    defaultValues.numSamplesPerBatch     = 100;
    defaultValues.numIterationsPerSample = 20;
    defaultValues.spikeSize              = 1;
    defaultValues.meanRateLearning       = false;
    defaultValues.stdpEnvelopeShape      = 'expDecay_continuous';
    defaultValues.stdpTimeOffset         = false;
    defaultValues.precisionHint          = 'double';
    defaultValues.outputCellGroupName    = [];
    defaultValues.refNumCells            = [];          % must be supplied if pct cell count is used
    defaultValues.cgDefaults             = struct();
    defaultValues.displayTimeUnits       = 'iteration'; % time units for display interval calculation
    defaultValues.debugDisplayEvery      = 0;           % TODO obsolete; need to phase out
    defaultValues.stats.printStatsEvery      = 1000;
    defaultValues.stats.updateFiguresEvery   = [];      % default set from printStatsEvery
    defaultValues.stats.numSTASamples        = 400;
    defaultValues.stats.keepSpikeHistory     = false;   % retain spike history? (for debugging)
    defaultValues.stats.keepPotentialHistory = false;   % retain membrane potential history?
    defaultValues.stats.measure              = struct();
    defaultValues.stats.print                = struct();
    defaultValues.stats.figure               = [];

    % apply the default model parameters to params
    params = ApplyDefaultValues(params, defaultValues);

    % default for updateFiguresEvery is printStatsEvery
    if isempty(params.stats.updateFiguresEvery)
        params.stats.updateFiguresEvery = params.stats.printStatsEvery;
    end


    % ===================   INITIALIZE MODEL STRUCT   ======================

    % copy parameters into the simulation model data structure
    % This copies all parameters into the model, including parameters
    % that are only used during initialization like cg_<name> and cg_<name>.in_<srcName>.
    model = params;
    model.initParams              = params;  % save copy of initialization params
    model.cellGroup               = [];
    model.displayTimeScale        = 1;       % # iterations per time unit for conversions
    model.debugDisplayWait        = 0;       % # displayTimeUnits since last stats display
    model.stats.areTitlesPrinted    = false;
    model.stats.areFiguresShown     = false;
    model.stats.printStatsWait      = 0;
    model.stats.updateFiguresWait   = 0;
    model.stats.trainingSampleCount = 0;

    % select time scaling factor for display intervals (params.displayTimeScale)
    switch model.displayTimeUnits
    case 'iteration'
        % to match SAILnet's default of 5 time constants per image * 100 images/batch
        model.displayTimeScale = 1;

    case 'sample'
        % applies to simulations that use input sample units, won't work for continuous time
        model.displayTimeScale = model.numIterationsPerSample;

    case 'simTimeUnit'
        % time units convert between simulated real-time and iterations
        error('simulation time units are not implemented yet');    % TODO need to implement

    otherwise
        error('unsupported value for displayTimeUnits');
    end

    % get a list of all 'cg_' (cell group) parameters
    paramNames   = fieldnames(params);
    cgParamNames = paramNames( strncmp(paramNames, 'cg_', 3) );


    % ===================   INITIALIZE CELL GROUPS   ===================

    % iterate through 'cg_' parameters to create cell groups (pass 1)
    for i = 1:numel(cgParamNames)
    
        % get the sub-struct representing the parameters for this cell group
        cgParams = params.(cgParamNames{i});
        if ~isstruct(cgParams)
            error(['non-struct encountered for parameter ' cgParamNames{i} ]);
        end

        % initialize the cell group
        if isfield(cgParams, 'numCells') && cgParams.numCells > 0;

            % initialize the cell group parameter struct
            cgParams      = ApplyDefaultValues(cgParams, params.cgDefaults, 'underlay');
            cgParams.name = cgParamNames{i}(4:end);

            % special case: numCells is a ratio relative to params.refNumCells
            if cgParams.numCells < 1
                if isempty(params.refNumCells)
                    error(['''' cgParams.name ''' specifies a numCells ratio, ' ...
                            'but params.refNumCells was not specified.']);
                end
                cgParams.numCells = max(round(cgParams.numCells * params.refNumCells), 1);
            end

            % add a cell group to the model
            model.cellGroup{end+1} = CellGroup_Init(cgParams);
        end
    end


    % ===================   INITIALIZE INPUT BLOCKS   ===================

    % iterate through 'cg_' parameters to create connections (pass 2)
    for i = 1:numel(cgParamNames)
    
        % get the sub-struct containing the parameters for this cell group
        cgParams = ApplyDefaultValues(params.(cgParamNames{i}), params.cgDefaults);

        % get the target cell group index
        cgName = cgParamNames{i}(4:end);
        i_cg   = FindByName(model.cellGroup, cgName);
        if isempty(i_cg)
            continue;                     % properties for missing cell groups are ignored
        end

        % process the input block descriptors - param name = cg_<name>.in_<srcName>
        paramNames   = fieldnames(cgParams);
        ibParamNames = paramNames( strncmp(paramNames, 'in_', 3) );
        for j = 1:numel(ibParamNames)

            % get the input block parameters (or input count)
            ibParams = cgParams.(ibParamNames{j});
            
            % get the index of the source cell group. The name if the source cell group
            % is derived from the last part of the property name, or can optionally
            % be overridden by specifying in_xxx.sourceName.
            if isfield(ibParams, 'sourceName')
                srcName = ibParams.sourceName;         % explicit sourceName overrides input name
            else
                srcName = ibParamNames{j}(4:end);      % default is to derive source name from input name
            end
            srcIndex = FindByName(model.cellGroup, srcName);
            if isempty(srcIndex)                
                % can't find the source cell group, so just skip this input descriptor.
                % If the cell group was not explicitly set to 0, then print an error.
                if isfield(params, ['cg_' srcName]) ...
                        && isfield(params.(['cg_' srcName]), 'numCells') ...
                        && params.(['cg_' srcName]).numCells
                    fprintf('Input cell group ''%s'' not found for target ''%s''. Skipping...\n', ...
                            srcName, cgName);
                end
                continue;             % inputs from missing cell groups are ignored
            end

            % determine the number of inputs per cell. If the input block descriptor
            % is a single number, then that is used. If it is a struct, then the
            % subparameter 'numInputs' is used. If missing but fully connected,
            % don't require. If a value < 1, calculate as a percentage of either the
            % total available inputs or refNumInputs.
            if isnumeric(ibParams)
                ibParams = struct('numInputs', ibParams);
            elseif ~isfield(ibParams, 'numInputs') || isempty(ibParams.numInputs)
                if ~isfield(ibParams, 'connectionPattern')
                    ibParams.connectionPattern = 'full';
                end
                if any(strcmpi(ibParams.connectionPattern, {'full','fullIndirect'}))
                    ibParams.numInputs = model.cellGroup{srcIndex}.numCells;
                else
                    error(['numInputs is required for ' cgParamNames{i} '.' ibParamNames{j} ]);
                end
            end
            if ibParams.numInputs == 0
                fprintf('WARNING: cellgroup %s, input block %s specified with numInputs=0; ignoring.\n', ...
                        cgName, srcName);
                continue;                       % no inputs. skip.
            elseif ibParams.numInputs > 0 && ibParams.numInputs < 1
                if isfield(cgParams, 'refNumInputs') && ~isempty(cgParams.refNumInputs)
                    refNumInputs = cgParams.refNumInputs;               % use cg_xxx.refNumInputs
                else
                    refNumInputs = model.cellGroup{srcIndex}.numCells;  % use max available inputs
                end                    
                ibParams.numInputs = max(round(ibParams.numInputs * refNumInputs), 1);
            elseif ibParams.numInputs >= 1
                ibParams.numInputs = round(ibParams.numInputs);
            else
                error('the number of inputs per cell cannot be negative');
            end

            % clamp self-to-self weights to zero if this connection is from a cell group to itself
            if ~isfield(ibParams, 'noSelfConnections')
                ibParams.noSelfConnections = (srcIndex == i_cg);
            elseif ibParams.noSelfConnections && (srcIndex ~= i_cg)
                error('can only disable self-connections for self-connected cell groups');
            end

            % create the new input block and add it to the cell group
            ibParams.sourceCellGroup = model.cellGroup{srcIndex};
            model.cellGroup{i_cg}    = CellGroup_AddInput(model.cellGroup{i_cg}, ibParams);
        end
    end


    % =============   POST-NETWORK MODEL STATE INITIALIZATION   ===============

    % fill in inputBlock cell group pointers ('srcId' field for fast lookup later)
    for i = 1:numel(model.cellGroup)
        for b = 1:numel(model.cellGroup{i}.inputBlock)
            srcName = model.cellGroup{i}.inputBlock(b).name;
            model.cellGroup{i}.inputBlock(b).srcId = FindByName(model.cellGroup, srcName);
        end
        
        % delayLine cell group pointer
        if isfield(model.cellGroup{i}, 'delayLine')
            srcName = model.cellGroup{i}.delayLine.srcName;
            model.cellGroup{i}.delayLine.srcId = NetModel_FindElement(model, srcName);
        end
    end

    % look up outputCellGroupName (if provided)
    model.outputCellGroupId = 0;
    if ~isempty(model.outputCellGroupName)
        model.outputCellGroupId = FindByName(model.cellGroup, model.outputCellGroupName);
    end

    % calculate the number of targets of each source cell
    % (this number is needed when interpreting the causalFeedback signal)
    model = UpdateTargetCounts(model);


    % ====================   DEBUG OUTPUT   =====================

    fprintf('Created network model with %d cell groups:\n', numel(model.cellGroup));
    for i = 1:numel(model.cellGroup)
        cg = model.cellGroup{i};
        if strcmp(cg.cellType, 'disabled')
            flag = '*';
        else
            flag = '';
        end
        %fprintf('    %d) %s%s -- ', i, cg.name, flag);
        fprintf('    %d) %s%s (n=%d) -- ', i, cg.name, flag, cg.numCells);
        if cg.isExternal
            fprintf('no inputs (external source)');
        elseif isfield(cg, 'delayLine')
            fprintf('delay line (%d iterations) from (%d)%s', ...
                    cg.delayLine.len, cg.delayLine.srcId, model.cellGroup{cg.delayLine.srcId}.name);
        elseif numel(cg.inputBlock) < 1
            fprintf('no inputs (possible error?)');
        else
            fprintf('inputs from: ');
            for b = 1:numel(cg.inputBlock)
                if strcmp(cg.inputBlock(b).connectionType, 'disabled')
                    flag = '*';
                else
                    flag = '';
                end
                if b > 1
                    fprintf(', ')
                end
                fprintf('(%d%s)%s', b, flag, cg.inputBlock(b).name);
            end
        end
        fprintf('\n');
    end
end


% UpdateTargetCounts - Update target counts for all cell groups
%
% Update target counts for each cell group in the sim model.  This
% will leave each cell group with an updated count of the number of
% cells it connects to.  The results are stored in the 'targetCount'
% property of each cell group.
%
% Usage:
%    model = UpdateTargetCounts(model)
%
% Inputs:
%    model    = the simulation model state
%
% Output:
%    model    = the updated simulation model state
%--------------------------------------------------------------------------
function model = UpdateTargetCounts( model )

    % iterate through each source cell group
    for i = 1:numel(model.cellGroup)

        srcName     = model.cellGroup{i}.name;
        targetCount = zeros(model.cellGroup{i}.numCells, 1);

        % iterate through each potential target cell group
        for j = 1:numel(model.cellGroup)

            % add up the target counts for each input block
            for b = 1:numel(model.cellGroup{j}.inputBlock)

                % skip if this input block is not receiving from the source cell group
                if ~strcmp(model.cellGroup{j}.inputBlock(b).name, srcName)
                    continue;
                end

                % add the number of targets for each source cell to the running total
                inputIndices = model.cellGroup{j}.inputBlock(b).inputIndices;
                if isempty(inputIndices)
                    targetCount = targetCount + model.cellGroup{j}.inputBlock(b).numInputs;
                else
                    for k = 1:numel(targetCount)
                        targetCount(k) = targetCount(k) + sum( inputIndices(:) == k );
                    end
                end
            end
        end

        % save the accumulated results into the source cell group
        model.cellGroup{i}.targetCount = targetCount;
    end
end

