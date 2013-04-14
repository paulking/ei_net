% NetModel_InitEINet - Initialize a network model with E and I populations
%
% Initialize a basic V1 network circuit model that uses two populations
% of neurons:
%    E - excitatory cells with input from image, I, E (E->E disabled by default)
%    I - inhibitory cells with input from E and I
%
% The connections from E -> E are normally disabled. When enabled, these
% are inhibitory connections to match the lateral inhibition model of SAILnet.
%
% Usage:
%    model = NetModel_InitEINet(params)
%
% Inputs:
%    params      = the configuration parameters:
%        modelSubtype        = which variation of the EI network is desired:
%            'jneuroscience'    = model used in J Neuroscience paper (CM rule, STDP)
%            'SAILnet'          = SAILnet - only one neuron layer, F rule
%            'fast'             = use mean rate learning on a smaller network
%                                 for speed as a demonstration (default)
%        learningRate        = learning rate to apply to the model as a whole. If 0,
%                              learning is disabled for the whole model. (default = 1).
%        stdpEnvelopeShape   = how to average spikes for input to learning; ignored
%                              if meanRateLearning = true. (default = 'expDecay_continuous')
%            'expDecay_continuous'     = use exponential moving average
%            'expDecay_continuous_all' = use exponential moving average, including for static inputs
%            'gaussian'                = use a temporal gaussian weighting
%        [see NetModel_Init for more model parameters]
%
% Output:
%    model       = the initialized network model
%
% See also:
%    NetModel_Init
%
% Created:   12/15/11, Paul King
%--------------------------------------------------------------------------
function model = NetModel_InitEINet( params )

    % initialize and apply default parameter values
    if nargin < 1
        params = struct();
    end
    defaultValues = struct();

    % general default values
    defaultValues.modelType                  = 'EINet';
    defaultValues.modelSubtype               = 'fast';
    defaultValues.numSamplesPerBatch         = 100;
    defaultValues.numIterationsPerSample     = 50;
    defaultValues.inputDims                  = [10,10];
    defaultValues.inputPreprocess.splitInvert  = false;
    defaultValues.inputPreprocess.rectify      = 'none';
    defaultValues.inputPreprocess.poissonScale = 0;          % poisson spiking (spikes-per-timeUnit)
    defaultValues.inputScale                 = (1/5);        % for compatibility with SAILnet (# spikes / 5 timeUnits)
    defaultValues.outputCellGroupName        = 'E';
    defaultValues.simMethod                  = 'cellgroupFastBatch';
    defaultValues.learningRate               = .5;
    defaultValues.learningRateUnits          = 'per100Samples';
    defaultValues.lrateScale                 = [];           % [derived from learningRateUnits]
    defaultValues.simTimeStep                = .1;           % to match SAILnet paper
    defaultValues.simTimeUnitSize            = .01;          % default iteration duration is 1 ms
    defaultValues.spikeSize                  = 1;            % SAILnet uses .1
    defaultValues.meanRateLearning           = false;
    defaultValues.stdpEnvelopeShape          = 'expDecay_continuous';
    defaultValues.displayTimeUnits           = 'sample';     % time units for display frequency
    defaultValues.autoStats                  = {'weights'};
    defaultValues.autoFigures                = [];

    % default values shared by all cell groups
    defaultValues.cgDefaults.location        = 'tile2D';     % (not used)
    defaultValues.cgDefaults.membraneTC      = 1;
    defaultValues.cgDefaults.spikeDecayRate  = 2;
    defaultValues.cgDefaults.initSpikeThresh = 2;
    defaultValues.cgDefaults.threshAdaptRate = 1;
    defaultValues.cgDefaults.targetSpikeRate = .05;          % spikes-per-timeUnit
    defaultValues.cgDefaults.spikeRate.instantWindow = 5;    % # iterations
    defaultValues.cgDefaults.spikeRate.meanWindow    = 5000; % # iterations

    % define external input cells
    defaultValues.cg_input.numCells          = [];           % derived from inputDims (below)
    defaultValues.cg_input.isExternal        = true;         % source is input image
    defaultValues.cg_input.cellType          = 'continuous';

    % define E excitatory cells
    defaultValues.cg_E.numCells              = 121;
    defaultValues.cg_E.cellType              = 'excitatory';
    defaultValues.cg_E.displayColor          = [0 1 0];      % green
    defaultValues.cg_E.initSpikeThresh       = 1;

    % define I inhibitory cells
    defaultValues.cg_I.numCells              = 36;
    defaultValues.cg_I.cellType              = 'inhibitory';
    defaultValues.cg_I.displayColor          = [1 0 0];      % red

    % add connections input->E (continuous) - image input
    defaultValues.cg_E.in_input.connectionType    = 'continuous';
    defaultValues.cg_E.in_input.connectionPattern = 'full';
    defaultValues.cg_E.in_input.learningRule      = 'hebbian_oja';
    defaultValues.cg_E.in_input.initWeights       = 'uniform';
    defaultValues.cg_E.in_input.constrainWeights  = 'none';
    defaultValues.cg_E.in_input.learningRate      = .2;

    % add connections E->E (inhibitory) - to decorrelate E cells directly (DISABLED)
    defaultValues.cg_E.in_E.connectionType        = 'disabled';  % 'inhibitory' to enable
    defaultValues.cg_E.in_E.connectionPattern     = 'full';
    defaultValues.cg_E.in_E.learningRule          = 'foldiak';
    defaultValues.cg_E.in_E.initWeights           = 'zero';
    defaultValues.cg_E.in_E.constrainWeights      = 'nonneg';
    defaultValues.cg_E.in_E.noSelfConnections     = true;
    defaultValues.cg_E.in_E.learningRate          = 2;

    % add connections E->I (excitatory) - to allow I to detect E correlations
    defaultValues.cg_I.in_E.connectionType        = 'excitatory';
    defaultValues.cg_I.in_E.connectionPattern     = 'full';
    defaultValues.cg_I.in_E.learningRule          = 'correlation_measuring';
    defaultValues.cg_I.in_E.initWeights           = 'uniform';
    defaultValues.cg_I.in_E.constrainWeights      = 'nonneg';
    defaultValues.cg_I.in_E.learningRate          = .7;

    % add connections I->E (inhibitory) - to decorrelate E cells using I input
    defaultValues.cg_E.in_I.connectionType        = 'inhibitory';
    defaultValues.cg_E.in_I.connectionPattern     = 'full';
    defaultValues.cg_E.in_I.learningRule          = 'correlation_measuring';
    defaultValues.cg_E.in_I.initWeights           = 'zero';
    defaultValues.cg_E.in_I.constrainWeights      = 'nonneg';
    defaultValues.cg_E.in_I.learningRate          = .7;

    % add connections I->I (inhibitory) - to decorrelate I cells from each other
    defaultValues.cg_I.in_I.connectionType        = 'inhibitory';
    defaultValues.cg_I.in_I.connectionPattern     = 'full';
    defaultValues.cg_I.in_I.learningRule          = 'correlation_measuring';
    defaultValues.cg_I.in_I.initWeights           = 'zero';
    defaultValues.cg_I.in_I.constrainWeights      = 'nonneg';
    defaultValues.cg_I.in_I.noSelfConnections     = true;
    defaultValues.cg_I.in_I.learningRate          = 1.5;


    % --------------------  DEFINE STATS MEASUREMENTS  ---------------------

    defaultValues.stats.measure.resErr_out.measureType    = 'resError';
    defaultValues.stats.measure.resErr_out.sourceName     = 'cg_output.in_input';

    defaultValues.stats.measure.corr_out.measureType      = 'correlation';
    defaultValues.stats.measure.corr_out.sourceName       = 'cg_output';

    defaultValues.stats.measure.spikeRate_out.measureType = 'spikeRate';
    defaultValues.stats.measure.spikeRate_out.sourceName  = 'cg_output';


    % --------------------  DEFINE STATS PRINT COLUMNS  ---------------------

    defaultValues.stats.printStatsEvery              = 2000;

    defaultValues.stats.print.numImages.position     = 1;
    defaultValues.stats.print.numImages.title        = '#images';
    defaultValues.stats.print.numImages.builtin      = 'sampleCount';
    defaultValues.stats.print.numImages.formatTitle  = '%7s';
    defaultValues.stats.print.numImages.formatData   = '%6d ';

    defaultValues.stats.print.spikeRateE.position    = 2;
    defaultValues.stats.print.spikeRateE.title       = 'sRateE';
    defaultValues.stats.print.spikeRateE.measureName = 'spikeRate_out';

    defaultValues.stats.print.threshE.position       = 3;
    defaultValues.stats.print.threshE.title          = 'threshE';
    defaultValues.stats.print.threshE.matlabExpr     = 'mean(cg_E.spikeThresh)';

    defaultValues.stats.print.threshI.position       = 3.1;
    defaultValues.stats.print.threshI.title          = 'threshI';
    defaultValues.stats.print.threshI.matlabExpr     = 'mean(cg_I.spikeThresh)';

    defaultValues.stats.print.ecorr.position         = 5;
    defaultValues.stats.print.ecorr.title            = 'E-corr';
    defaultValues.stats.print.ecorr.measureName      = 'corr_out';

    defaultValues.stats.measure.dW.measureType       = 'deltaWeight';
    defaultValues.stats.print.dW.position            = 6;
    defaultValues.stats.print.dW.title               = 'dW_rms';
    defaultValues.stats.print.dW.measureName         = 'dW';

    defaultValues.stats.print.resErr.position        = 7;
    defaultValues.stats.print.resErr.title           = 'resErr';
    defaultValues.stats.print.resErr.measureName     = 'resErr_out';


    % --------------------  SECOND PASS DEFAULT INITIALIZATION  ---------------------

    % determine selected model subtype
    modelSubtype = defaultValues.modelSubtype;
    if isfield(params, 'modelSubtype')
        modelSubtype = params.modelSubtype;
    end

    % set second-pass defaults according to model subtype
    switch modelSubtype

    case 'jneuroscience'
        % parameters for default network presented in J Neuroscience paper
        defaultValues.learningRate                   = .4;
        defaultValues.meanRateLearning               = false;    % use STDP
        defaultValues.stdpEnvelopeShape              = 'expDecay_continuous';
        defaultValues.autoStats                      = {'weights', 'corr_EI', 'corr_I', 'spikeRate_I'};
        defaultValues.cgDefaults.spikeRate.instantWindow = 10;
        defaultValues.cg_E.numCells                  = 400;
        defaultValues.cg_E.targetSpikeRate           = .02;
        defaultValues.cg_E.membraneTC                = 1;
        defaultValues.cg_E.initSpikeThresh           = 2;        % needed for higher input blockWeight
        defaultValues.cg_E.in_input.learningRule     = 'hebbian_oja';
        defaultValues.cg_E.in_input.blockWeight      = 5;        % for compatibility with SAILnet input scaling
        defaultValues.cg_I.numCells                  = 49;       % use 49 for .5x overcomplete, 100 for 1x
        defaultValues.cg_I.targetSpikeRate           = .04;      % I is 2x faster than E
        defaultValues.cg_I.membraneTC                = .5;
        defaultValues.stats.numSTASamples            = 2000;
        defaultValues.stats.printStatsEvery          = 10000;    % less frequent figure updates
        defaultValues.stats.measure.dW.deltaInterval = 10000;
        defaultValues.meanRateLearning               = false;

    case 'SAILnet'
        % override params to disable I cells and enable lateral E cell inhibition
        defaultValues.meanRateLearning             = true;
        defaultValues.cg_E.targetSpikeRate         = .025;
        defaultValues.cg_E.initSpikeThresh         = 2;            % needed for higher input blockWeight
        defaultValues.cg_E.membraneTC              = -.1/log(1-.1);% .9491 to match exactly SAILnet "TC" = 1
        defaultValues.cg_E.in_input.blockWeight    = 5;            % for compatibility with SAILnet input scaling
        defaultValues.cg_E.in_I.connectionType     = 'disabled';   % disable connection: I->E
        defaultValues.cg_E.in_E.connectionType     = 'inhibitory'; % ENABLE lateral inhibition: E->E
        defaultValues.cg_I.cellType                = 'disabled';   % disable I cells
        defaultValues.stats.print.threshI.position = -1;         % disable I thresh print column

    case 'fast'
        % faster network that uses fewer iterations per sample and mean rate learning
        defaultValues.numIterationsPerSample       = 20;
        defaultValues.meanRateLearning             = true;
        defaultValues.learningRate                 = .2;
        defaultValues.autoStats                    = {'weights', 'corr_EI', 'corr_I', 'spikeRate_I'};
    end


    params = ApplyDefaultValues(params, defaultValues);

    assert(any( strcmp(params.simMethod, {'cellgroup','cellgroupFast', 'cellgroupFastBatch'}) ));


    % ====================   PRE-INITIALIZE MODEL PARAMETERS   ===================

    % set input and output cell counts
    numPixels = prod(params.inputDims);
    params.cg_input.numCells = numPixels;

    % generate locations of the input pixels for geolocal connectivity later
    [X_input,Y_input] = meshgrid( ((1:params.inputDims(2)) - .5) / params.inputDims(2), ...
                      ((1:params.inputDims(1)) - .5) / params.inputDims(1) );
    location_input    = cat(2, Y_input(:), X_input(:));

    % special case: factored input
    if params.inputPreprocess.splitInvert
        params.cg_input.numCells                = numPixels * 2;  % double pixels for pos/neg split input
        location_input = cat(1, location_input, location_input);
    end

    % a derived constant for convenience (TODO move to NetModel_Init?)
    params.numTimeUnitsPerSample = params.numIterationsPerSample * params.simTimeStep;

    % select scaling factor for learning rates (params.lrateScale)
    % (note: per X timeUnits calculation isn't right, because per100Samples assumes
    % 2 time units per sample)
    switch params.learningRateUnits
    case 'per1000TimeUnits'
        % a round number, assuming simTimeStep = 1 (as CellGroup default spiking networks)
        params.lrateScale = 1 / 1000;

    case 'per100TimeUnits'
        % a round number, assuming simTimeStep = .1 (as in SAILnet)
        params.lrateScale = 1 / 100;

    case 'per100Samples'
        % compatible with SAILnet paper, independent of batch size grouping
        params.lrateScale = 1 / 100;

    case 'perBatch'
        % exactly matches SAILnet, but blows up for small batches
        params.lrateScale = 1 / params.numSamplesPerBatch;
        if params.numSamplesPerBatch < 100
            fprintf('WARNING: Using ''perBatch'' learning rate units with small batch size!\n');
        end

    otherwise
        error('unsupported value for learningRateUnits');
    end

    % fill in secondary (derived) default values: spike threshold & rate, learning rates
    numIterationsPerBatch = params.numIterationsPerSample * params.numSamplesPerBatch;
    defaultValues2 = struct();
    defaultValues2.cgDefaults.updateEvery       = numIterationsPerBatch;
    defaultValues2.cgDefaults.updateThreshEvery = numIterationsPerBatch;
    defaultValues2.cg_input.location            = location_input;
    params = ApplyDefaultValues(params, defaultValues2, 'underlay');


    % ======================   INITIALIZE THE MODEL   ======================

    fprintf('Using simulation method "%s"\n', params.simMethod);

    % initialize the base network simulation model
    model = NetModel_Init(params);


    % ==================   ADD AUTO-GENERATED PRINT COLUMNS   ==================

    % add auto-generated statistics columns
    for k = 1:numel(model.autoStats)
        cg = model.cellGroup;
        switch model.autoStats{k}

        case 'weights'
            % add columns for each set of connection weights (each active input block)
            for i = 1:numel(cg)
                if cg{i}.isExternal || strcmp(cg{i}.cellType, 'disabled')
                    continue;
                end
                for j = 1:numel(cg{i}.inputBlock)
                    ib = cg{i}.inputBlock(j);
                    if strcmp(ib.connectionType, 'disabled')
                        continue;
                    end

                    if strcmp(ib.name, 'input')
                        inName = 'in';
                    else
                        inName = ib.name;
                    end
                    meanName = sprintf('Wmean_%d_%d', i, j);
                    srcName  = ['cg_' cg{i}.name '.in_' ib.name];
                    model.stats.print.(meanName).position    = 4 + i/10 + j/100;
                    model.stats.print.(meanName).title       = [inName '->' cg{i}.name];
                    model.stats.print.(meanName).matlabExpr  = ['mean(abs(' srcName '.weight(:)))'];
                    model.stats.print.(meanName).formatTitle = '%9s ';
                    model.stats.print.(meanName).formatData  = '%9.4f ';
                end
            end

        case 'spikeRate_I'
            % add print column for I cell spikeRate
            model.stats.measure.spikeRateI.measureType = 'spikeRate';    % TODO just use built-in popMean?
            model.stats.measure.spikeRateI.sourceName  = 'cg_I';
            model.stats.print.spikeRateI.position      = 2.1;
            model.stats.print.spikeRateI.title         = 'sRateI';
            model.stats.print.spikeRateI.measureName   = 'spikeRateI';

        case 'corr_EI'
            % add print column for EI_corr (correlation between E and I cells)
            model.stats.measure.corr_Ei.measureType = 'correlation';
            model.stats.measure.corr_Ei.sourceName  = {'cg_E', 'cg_I'};
            model.stats.print.eicorr.position         = 5.1;
            model.stats.print.eicorr.title            = 'EI-corr';
            model.stats.print.eicorr.measureName      = 'corr_Ei';

        case 'corr_I'
            % add print column for I_corr (self-correlation within I population)
            model.stats.measure.corr_I.measureType  = 'correlation';
            model.stats.measure.corr_I.sourceName   = 'cg_I';
            model.stats.print.icorr.position          = 5.2;
            model.stats.print.icorr.title             = 'I-corr';
            model.stats.print.icorr.measureName       = 'corr_I';

        case 'resErr_I'
            % add print column for resErr_I (does I reconstruct E?)
            model.stats.measure.resErr_I.measureType = 'resError';
            model.stats.measure.resErr_I.sourceName  = 'cg_I.in_E';
            model.stats.measure.resErr_I.normalize   = true;     % because of spiking input
            model.stats.print.resErr_I.position      = 6.9;
            model.stats.print.resErr_I.title         = 'resErr_I';
            model.stats.print.resErr_I.measureName   = 'resErr_I';
            model.stats.print.resErr_I.formatTitle   = '%11s';
            model.stats.print.resErr_I.formatData    = '%9.4f  ';
        end
    end

    % add auto-generated figures
    % (static figure numbers are used so that the same figure type displays in the same place)
    for k = 1:numel(model.autoFigures)
        switch model.autoFigures{k}

        case 'STA_E'
            model.stats.figure{end+1}.figureNum  = 1001;
            model.stats.figure{end}.plotType     = 'STA';
            model.stats.figure{end}.sourceName   = 'cg_E';

        case 'STA_I'
            i_cg = NetModel_FindElement(model, 'cg_I');
            if strcmp(model.cellGroup{i_cg}.cellType, 'disabled')
                continue;
            end
            model.stats.figure{end+1}.figureNum  = 1002;
            model.stats.figure{end}.plotType     = 'STA';
            model.stats.figure{end}.sourceName   = 'cg_I';

        case 'spikeRaster'
            model.stats.figure{end+1}.figureNum  = 1003;
            model.stats.figure{end}.plotType     = 'spikeRasterAllImage';

        case 'networkDynamics'
            model.stats.figure{end+1}.figureNum  = 1004;
            model.stats.figure{end}.plotType     = 'networkDynamics';
            model.stats.figure{end}.YMin         = -2;
            model.stats.figure{end}.markerType   = 's';
            if ~ model.stats.keepPotentialHistory
                fprintf('turning on model.stats.keepPotentialHistory to support networkDynamics plot\n');
                model.stats.keepPotentialHistory = true;
            end

        case 'weightDist'
            model.stats.figure{end+1}.figureNum  = 1005;
            model.stats.figure{end}.plotType     = 'weightDistributionAll';

        case 'weightCorr'
            model.stats.figure{end+1}.figureNum  = 1006;
            model.stats.figure{end}.plotType     = 'weightCorrelation';
            model.stats.figure{end}.sourceName   = {'cg_I.in_E', 'cg_E.in_I'};
        end
    end
end
