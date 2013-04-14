% AnalyzeV1eV1iPSTH - Analyze PSTH for E-I Net (VeV1i Network)
%
% Generate a set of histogram metrics for computing PSTH histograms of the
% relationship between E and I cells.
%
% Note: histMeanStrong + histMeanWeak does not necessarily equal sum(histMean,1), because
% middle deciles may be excluded from averaging.
%
% Inputs:
%    model        = the model (read-only)
%    params       = configuration parameters
%        windowDev        = amount +/- of temporal window in iterations (default=15)
%        numStrongDeciles = number of deciles to average as the "strong targets" (default = 2)
%        numWeakDeciles   = number of deciles to average as the "weak targets" (default = 5)
%
% Outputs:
%    results      = struct containing output results
%        timeBinVal     = the delta-t location of each time bin
%        strengthBinVal = the mean connection strength of search connection strength bin
%        histMean       = matrix(numStrengthBins,numTimeBins) representing I cell
%                         spike triggered E cell spike rate histogram.
%                         Each element is a histogram bin indicating average E cell
%                         spike rate as a multiple of the target E cell spike rate
%                         (1 = target rate). numStrengthBins == 7, where the first 6 bins
%                         are the top 6 deciles and the 7th bin aggregates the bottom
%                         4 deciles (which are often all zero).
%        histMeanStrong = matrix(1,numTimeBins) of histogram bins for the E cells represented
%                         by the 20% strongest I->E connections.
%        histMeanWeak   = matrix(1,numTimeBins) of histogram bins for the E cells represented
%                         by the 80% weakest I->E connections.
%        histInhibMean  = array(1,numTimeBins) representing the total inhibitory feedback
%                         to E cells I->E (on average) for each E-spike-triggered
%                         relative time bin.
%        histIInhibMean = array(1,numTimeBins) representing the total inhibitory feedback
%                         to I cells I->I (on average) for each I-spike-triggered
%                         relative time bin.
%
% Created:   11/11/2012, Paul King
%--------------------------------------------------------------------------
function results = AnalyzeV1eV1iPSTH( model, params )

    % initialize and apply default parameter values
    if nargin < 2
        params = struct();
    end
    defaultValues = struct();

    % general default values
    defaultValues.sourceName              = 'cg_V1e.in_V1i';
    defaultValues.sourceNameII            = 'cg_V1i.in_V1i';
    defaultValues.windowDev               = 15;
    defaultValues.numStrongDeciles        = 2;    % # top deciles to include in 'strong' metric
    defaultValues.numWeakDeciles          = 5;    % # bottom deciles to include in 'weak' metric
    params = ApplyDefaultValues(params, defaultValues);


    % ==================   PRECALCULATE HISTOGRAM INFO   =================

    % get I->E connection weights
    [i_cg_E, i_ib_EI] = NetModel_FindElement(model, params.sourceName);
    bw_I_to_E         = model.cellGroup{i_cg_E}.inputBlock(i_ib_EI).blockWeight;
    weight_I_to_E     = bw_I_to_E * model.cellGroup{i_cg_E}.inputBlock(i_ib_EI).weight;
    [numE, numI]      = size(weight_I_to_E);

    % bin by deciles globally (last bin has 4 deciles = 40%)
    weightDecile    = prctile(weight_I_to_E(:), 10:10:100);
    if weightDecile(1) == 0 && min(weight_I_to_E(:)) == 0
        numStrengthBins = 11 - find(weightDecile, 1);
    else
        numStrengthBins = 10;
    end
    strengthBinIdx  = zeros(size(weight_I_to_E));
    for i = 1:numStrengthBins
        strengthBinIdx(weight_I_to_E <= weightDecile(11-i)) = i;
    end

    % number of iterations on either side of spike trigger to evaluate
    windowDev = params.windowDev;

    % allocate histogram bins
    numTimeBins = windowDev*2 + 1;
    histSum     = zeros(numStrengthBins, numTimeBins);
    histCount   = zeros(numStrengthBins, numTimeBins);

    %{
    % experimental: bin by deciles locally
    % (arguably, this is more correct; but should we sort within E->I or I->E)
    % (here, we sort by E->I, i.e. we sort for each I cell, all E inputs)
    strengthBinIdx = zeros(size(weight_I_to_E));
    for k = 1:numI
        weightDecile = prctile(weight_I_to_E(:,k), 10:10:100);
        for i = 1:numStrengthBins
            strengthBinIdx(weight_I_to_E(:,k) <= weightDecile(11-i), k) = i;
        end
    end
    %}


    % =======  GENERATE 2D I-TRIGGERED PSTH OF E SPIKERATE BINNED BY I->E CONNECTION STRENGTH  ======

    % average across samples to get the average spikes for each cell at each time step
    i_cg_E      = NetModel_FindElement(model, 'cg_V1e');
    i_cg_I      = NetModel_FindElement(model, 'cg_V1i');
    spikeHist   = model.snapshot.spikeHistory;
    spikeHist_E = spikeHist{i_cg_E};
    spikeHist_I = spikeHist{i_cg_I};
    [numCellsE, numSamples, numIterations] = size(spikeHist_E);

    % debug: plot average spike rate over time for E and I cells (to show oscillations)
    %figure(11); plot(1:numIterations, sum(spikeMean_E,1));  set(11, 'Name','E sumSpikes');
    %figure(12); plot(1:numIterations, sum(spikeMean_I,1));  set(12, 'Name','I sumSpikes');

%
    % for each shifted window center, accumulate histogram bin counts
    for windowCenter = (1+windowDev) : (numIterations-windowDev)
        spikeHist_Ebuffer   = spikeHist_E(:, :, (windowCenter-windowDev):(windowCenter+windowDev));
        spikeHist_Itrigger  = spikeHist_I(:, :, windowCenter);   % size = numI x numSamples

        % accumulate histogram bins according to E-I strength
        spikeHist_Ereshaped = reshape(spikeHist_Ebuffer, [], numTimeBins);
        for idxBin = 1:numStrengthBins
            selectedConnections = (strengthBinIdx == idxBin);   % size = numE x numI
            sumItriggers        = selectedConnections * spikeHist_Itrigger; % size = numE x numSamples
            histSum(idxBin,:)   = histSum(idxBin,:)   + sumItriggers(:)' * spikeHist_Ereshaped;
            histCount(idxBin,:) = histCount(idxBin,:) + sum(sumItriggers(:));
        end
    end
%}
    
%{
    % alternate implementation, 5x slower
    numCellsI = size(spikeHist_I, 1);
    for k = 1:numTimeBins
        i0 = k;
        i1 = numIterations - numTimeBins + k;
        s0 = 1 + windowDev;
        s1 = numIterations - windowDev;
        subSpikeHistE = reshape(spikeHist_E(:,:,i0:i1), numCellsE, []); % size = numE x numSamples*numSubIt
        subItrigger   = reshape(spikeHist_I(:,:,s0:s1), numCellsI, []); % size = numI x numSamples*numSubIt
        for idxBin = 1:numStrengthBins
            selectedConnections = (strengthBinIdx == idxBin);   % size = numE x numI
            sumItriggers        = selectedConnections * subItrigger; % size = numE x numSamples*numSubIt
            histSum(idxBin,k)   = sumItriggers(:)' * subSpikeHistE(:);
            histCount(idxBin,k) = sum(sumItriggers(:));
        end
    end
%}

%{
    % alternate implementation, only 5% faster
    % precalculate array variables
    spikeHist_Ereshaped = reshape(spikeHist{i_cg_E}, numCellsE*numSamples, numIterations);
    selectedConnections = [];
    for idxBin = 1:numStrengthBins
        selectedConnections{idxBin} = (strengthBinIdx == idxBin);   % size = numE x numI x numBins
    end

    % for each shifted window center, accumulate histogram bin counts
    for windowCenter = (1+windowDev) : (numIterations-windowDev)
        winRange            = (windowCenter-windowDev):(windowCenter+windowDev);
        spikeHist_Itrigger  = spikeHist_I(:, :, windowCenter);   % size = numI x numSamples
        spikeHist_Ebuffer   = spikeHist_Ereshaped(:,winRange);   % size = numE*numSamples x numBins

        % accumulate histogram bins according to E-I strength
        for idxBin = 1:numStrengthBins
            sumItriggers        = selectedConnections{idxBin} * spikeHist_Itrigger; % size = numE x numSamples
            histSum(idxBin,:)   = histSum(idxBin,:)   + sumItriggers(:)' * spikeHist_Ebuffer;
            histCount(idxBin,:) = histCount(idxBin,:) + sum(sumItriggers(:));
        end
    end
%}

    % normalize histogram by dividing by count (unless count = 0)
    histMean = histSum ./ histCount ./ model.simTimeStep;
    histMean(histCount == 0) = 0;

    % calculate bin values
    timeBinVal     = (1:numTimeBins) - windowDev - 1;
    strengthBinVal = zeros(numStrengthBins,1);
    for i = 1:numStrengthBins
        strengthBinVal(i) = mean(weight_I_to_E(strengthBinIdx == i));
    end

    % experiment: collapse histMean into bins containing two iterations
    % (to smooth how sawtoothed spike rate graph;  doesn't help much)
    %{
    histMean    = cat(2, histMean(:,1), ...
                   ( histMean(:,2:2:size(histMean,2)) + histMean(:,3:2:size(histMean,2)) )/2 );
    windowDev   = windowDev / 2;
    numTimeBins = round(numTimeBins / 2);
    %}

    % generate histogram of E cells connected by weak vs. strong connections
    sd = min(params.numStrongDeciles, numStrengthBins-1);
    wd = params.numWeakDeciles;
    histMeanStrong = mean(histMean(1:sd,:), 1);
    if wd > 11 - numStrengthBins
        histMeanWeak = (mean(histMean((11-wd):end-1,:), 1) ...
                + (11-numStrengthBins)*histMean(end,:)) / wd;     % last bin is larger
    else
        histMeanWeak = histMean(end,:);              % use bottom bin
    end


    % ==========  GENERATE E-TRIGGERED PSTH OF TOTAL I INPUT TO E CELL ==========

    % compute E-spike-triggered total I input to E cells
    histInhibSum   = zeros(1,numTimeBins);
    histInhibCount = 0;
    numSamples     = size(spikeHist_E, 2);
    inhibInput_E   = weight_I_to_E * reshape(spikeHist_I, numI, []);
    inhibInput_E   = reshape(inhibInput_E, [numE, numSamples, numIterations] );
    for windowCenter = (1+windowDev) : (numIterations-windowDev)
        inhib_Ebuffer       = inhibInput_E(:, :, (windowCenter-windowDev):(windowCenter+windowDev));
        inhib_Ereshaped     = reshape(inhib_Ebuffer, [], numTimeBins);  % size = numE*numSamples x numTimeBins
        spikeHist_Etrigger  = spikeHist_E(:, :, windowCenter);   % size = numE x numSamples
        histInhibSum        = histInhibSum   + spikeHist_Etrigger(:)' * inhib_Ereshaped;
        histInhibCount      = histInhibCount + sum(spikeHist_Etrigger(:));
    end
    histInhibMean = histInhibSum ./ histInhibCount;
    histInhibMean(histInhibCount == 0) = 0;


    % ==========  GENERATE I-TRIGGERED PSTH OF TOTAL I INPUT TO I CELL ==========

    % experimental
    if ~isempty(params.sourceNameII)
        [i_cg_I, i_ib_II] = NetModel_FindElement(model, params.sourceNameII);
        bw_I_to_I         = model.cellGroup{i_cg_I}.inputBlock(i_ib_II).blockWeight;
        weight_I_to_I     = bw_I_to_I * model.cellGroup{i_cg_I}.inputBlock(i_ib_II).weight;
        spikeHist_I       = spikeHist{i_cg_I};

        % compute I-spike-triggered total I input to I cells
        histInhibSum   = zeros(1,numTimeBins);
        histInhibCount = 0;
        numI           = model.cellGroup{i_cg_I}.numCells;
        numSamples     = size(spikeHist_I, 2);
        inhibInput_I   = weight_I_to_I * reshape(spikeHist_I, numI, []);
        inhibInput_I   = reshape(inhibInput_I, [numI, numSamples, numIterations] );
        for windowCenter = (1+windowDev) : (numIterations-windowDev)
            inhib_Ibuffer       = inhibInput_I(:, :, (windowCenter-windowDev):(windowCenter+windowDev));
            inhib_Ireshaped     = reshape(inhib_Ibuffer, [], numTimeBins);  % size = numI*numSamples x numTimeBins
            spikeHist_Itrigger  = spikeHist_I(:, :, windowCenter);   % size = numI x numSamples
            histInhibSum        = histInhibSum   + spikeHist_Itrigger(:)' * inhib_Ireshaped;
            histInhibCount      = histInhibCount + sum(spikeHist_Itrigger(:));
        end
        histIInhibMean = histInhibSum ./ histInhibCount;
        histIInhibMean(histInhibCount == 0) = 0;
    else
        histIInhibMean = [];
    end


    % =======================  RETURN RESULTS  =======================

    results = struct();
    results.timeBinVal     = timeBinVal;
    results.strengthBinVal = strengthBinVal;
    results.histMean       = histMean;
    results.histMeanStrong = histMeanStrong;
    results.histMeanWeak   = histMeanWeak;
    results.histInhibMean  = histInhibMean;
    results.histIInhibMean = histIInhibMean;
end

