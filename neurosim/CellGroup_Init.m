% CellGroup_Init - Initialize a group of computational units (neurons)
%
% Create a set of cells for wiring into networks.
%
% The cellType indicates the type of cell being created.  Supported
% types are 'excitatory' and 'inhibitory'.  Uniformly distributed (x,y)
% coordinates in the range (0,1) are assigned to each cell, which can
% be used to implement locally-connected wiring.
%
% For details on the struct properties of cell group, see CellGroup_Update.
%
% Usage:
%    cellGroup = CellGroup_Init(params)
%
% Inputs:
%    params      = the cell group configuration parameters:
%        numCells            = the number of cells to allocate
%        name                = name to give this cell group
%        cellType            = cell type (default = 'excitatory')
%            'excitatory'        = input has a positive effect on potential
%            'inhibitory'        = input has a negative effect on potential
%            'disabled'          = cellgroup is ignored in all calculations
%        location            = method for allocating cell location. (default = 'tile2D')
%            'uniform2D'         = assign each cell two random coordinates in the range (0,1)
%            'tile2D'            = evenly tile the 2D unit square
%            matrix              = if a numeric matrix(numCells,2) is provided,
%                                  use those locations
%        defLearningRule     = default learning rule for outputs from this cell group.
%                              Only used during wiring phase. (default = [])
%        defForceDirect      = whether or not to use direct connections by default
%                              (only used during wiring phase)
%        displayColor        = RGB color for graphs and plots (default = [1 1 1] = white)
%        membraneTC          = membrane time constant determining potential decay rate
%                              in time units (default = 10)
%        spikeDecayRate      = "spiked recently" decay rate in simulation time units
%                              (default = .2)
%        learningRate        = learning rate to apply to all connection weights in the cell group.
%                              suggested values are between 0 and 1, where 0 = no weight updates
%                              and 1 = recommended 'fast' update rate. spike threshold updates
%                              are not affected (see threshAdaptRate). (default = 1).
%        threshAdaptRate     = learning rate to apply to thresholds; weights are not affected
%        initSpikeThresh     = initial threshold for spiking (default = 1)
%        targetSpikeRate     = target spike rate (spikes/interation) for spike rate autoranging.
%                              0 disables autoranging. (default = 0)
%        trackCausalFeedback = whether or not to track causal feedback (default = false)
%        updateThreshEvery   = how often (# iterations) up update spike thresholds (default = 100)
%        updateEvery         = how often (# iterations) to update weights (default = 100)
%        rescaleEvery        = how often to rescale weights; 0 = never (default = 0)
%        spikeRate           = struct with various fields for spike rate tracking
%            meanWindow          = sampling window (#iterations) to use to compute
%                                  long-term mean rate of spiking for network monitoring
%                                  and reporting. This value has no effect on network
%                                  learning behavior. (default = 10000)
%            lmeanWindow         = sampling window (#iterations) to use to compute long-term
%                                  mean rate of spiking used by learning rules. Changing
%                                  this can alter network behavior. (default = 5000)
%            instantWindow       = sampling window (#iterations) to use for measurement
%                                  of instantaneous spike rate. (default = 5)
%        delayLine           = substruct of delay line params
%            len                 = length of delay (# iterations)
%
% Output:
%    cellGroup    = the allocated inhibitory cell population
%                   (See CellGroup_Update for struct property details.)
%
% See also:
%    CellGroup_Update
%
% Created:   11/15/09, Paul King
%--------------------------------------------------------------------------
function cellGroup = CellGroup_Init( params )

    % initialize and apply default parameter values
    if nargin < 1
        params = struct();
    end
    defaultValues.cellType            = 'excitatory';
    defaultValues.location            = 'tile2D';
    defaultValues.isExternal          = false;
    defaultValues.name                = '';
    defaultValues.defLearningRule     = [];
    defaultValues.defForceDirect      = true;
    defaultValues.displayColor        = [1 1 1];
    defaultValues.membraneTC          = 10;
    defaultValues.spikeDecayRate      = .2;
    defaultValues.learningRate        = 1;
    defaultValues.initSpikeThresh     = 1;
    defaultValues.targetSpikeRate     = .20;
    defaultValues.threshAdaptRate     = 1;               % relative to model.lrateScale
    defaultValues.trackCausalFeedback = false;
    defaultValues.updateEvery         = 40;              % update weights every n cycles
    defaultValues.updateThreshEvery   = 10;              % update thresholds every n cycles
    defaultValues.rescaleEvery        = 0;               % rescale weights every n cycles (disabled)
    defaultValues.spikeRate.instantWindow = 5;
    defaultValues.spikeRate.meanWindow    = 10000;
    defaultValues.spikeRate.lMeanWindow   = 5000;
    params = ApplyDefaultValues(params, defaultValues);


    % ===================   INITIALIZE CELL GROUP STRUCT   ======================

    % initialize the cell group data structure
    % (see CellGroup_Update for details on cell group structure elements)
    cellGroup = struct();
    cellGroup.name             = params.name;
    cellGroup.numCells         = params.numCells;
    cellGroup.cellType         = params.cellType;
    cellGroup.isExternal       = params.isExternal;
    cellGroup.displayColor     = params.displayColor;
    cellGroup.defLearningRule  = params.defLearningRule;
    cellGroup.defForceDirect   = params.defForceDirect;
    cellGroup.location         = [];
    cellGroup.potential        = [];
    cellGroup.spikeThresh      = params.initSpikeThresh * ones(params.numCells,1);
    cellGroup.dspikeThresh     = zeros(params.numCells,1);   % accumulated delta to spikeThresh
    cellGroup.spikedRecently   = zeros(params.numCells,1);
    cellGroup.causalFeedback   = [];
    cellGroup.targetCount      = [];
    cellGroup.inputBlock       = [];
    cellGroup.membraneTC       = params.membraneTC;
    cellGroup.spikeDecayRate   = params.spikeDecayRate;
    cellGroup.learningRate     = params.learningRate;
    cellGroup.targetSpikeRate  = params.targetSpikeRate;
    cellGroup.threshAdaptRate  = params.threshAdaptRate;
    cellGroup.updateThreshEvery  = params.updateThreshEvery;
    cellGroup.updateThreshWait   = 0;                    % # iteration ticks since last thresh rescaling
    cellGroup.updateEvery        = params.updateEvery;
    cellGroup.updateWait         = 0;
    cellGroup.rescaleEvery       = params.rescaleEvery;
    cellGroup.rescaleWait        = 0;                    % # iteration ticks since last weight rescaling
    cellGroup.spikeRate.meanWindow        = params.spikeRate.meanWindow;
    cellGroup.spikeRate.meanBiased        = zeros(params.numCells, 1);
    cellGroup.spikeRate.meanBiased_max    = 0;
    cellGroup.spikeRate.mean              = zeros(params.numCells, 1);
    cellGroup.spikeRate.lMeanWindow       = params.spikeRate.lMeanWindow;
    cellGroup.spikeRate.lMeanBiased       = zeros(params.numCells, 1);
    cellGroup.spikeRate.lMeanBiased_max   = 0;
    cellGroup.spikeRate.lMean             = zeros(params.numCells, 1);
    %cellGroup.spikeRate.lMean2Biased      = zeros(params.numCells, 1);
    %cellGroup.spikeRate.lMean2            = zeros(params.numCells, 1);    % mean(spikeRate^2)
    cellGroup.spikeRate.popMean           = 0;
    cellGroup.spikeRate.instantWindow     = params.spikeRate.instantWindow;
    cellGroup.spikeRate.instant           = zeros(params.numCells, 1);

    % initialize delayLine mode
    if isfield(params, 'delayLine')
        cellGroup.delayLine = params.delayLine;
        cellGroup.delayline.buffer = zeros(params.numCells, params.delayLine.len);
        cellGroup.delayline.idx    = 1;
    end

    % initialize causalFeedback store if this feature is turned on
    if params.trackCausalFeedback
        cellGroup.causalFeedback = zeros(params.numCells,1);
    end

    % assign euclidean positions in unit 2D space to the cells
    if isnumeric(params.location)
        % locations were passed in -- use these
        if ~isempty(params.location) && any( size(params.location) ~= [params.numCells, 2] )
            error('invalid dimensions for location matrix parameter');
        end
        cellGroup.location = params.location;
    else
        switch params.location

        case 'uniform2D'
            % assign random uniformly distributed locations
            cellGroup.location = rand(cellGroup.numCells,2);

        case 'tile2D'
            % assign geometrically tiled locations
            % (this will be a grid if numCells is square, and staggered otherwise)
            numCols = round(sqrt(cellGroup.numCells));
            idx     = (1:cellGroup.numCells)';
            col     = floor((idx - .5) * numCols / cellGroup.numCells);  % row # (zero-based)
            X       = (col + .5) / numCols;
            Y       = (idx - .5) * numCols / cellGroup.numCells - col;
            Y       = Y + (1 - max(Y) - min(Y)) / 2;              % adjust to equalize margins
            cellGroup.location = [Y, X];

        otherwise
            error('unknown location allocation strategy ''%s''', params.location);
        end
    end

    % initialize the internal state variables (unless this is a "dummy" cell group)
    if ~params.isExternal
        cellGroup.potential = zeros(params.numCells,1);
    end
end

