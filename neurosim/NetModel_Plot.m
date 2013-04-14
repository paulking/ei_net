% NetModel_Plot - Draw a chart based on model state (internal use)
%
% Figure supports the following subfields:
%    figureNum      = figure number to use. (default, assigned automatically
%                     starting at 1001)
%    title          = title to use to label the figure
%    sourceName     = a name identifying the signal source for the chart.
%                     Examples: 'cg_V1e' or 'cg_V1e.in_V1i'
%    plotType      = type of chart to draw in this figure (see plotTypes below)
%       'STA'                = spike-triggered average
%       'STA_movie'          = spike-triggered movie
%       'inputImagePatches'  = batch of input training samples (either still image patches
%                              or image patch movie clips)
%       'weightPatches'      = grid of 2D image patch basis functions
%       'weightMatrix'       = image showing 2D weight strengths, with scale
%       'weightDistribution' = histogram showing distribution of weight values
%       'weightDistributionAll' = histogram showing weight distributions for all input blocks
%       'weightCorrelation'  = scatter plot comparing recurrent A->B and B->A weights
%       'spikeRateHistogram' = histogram of spike rates across training samples
%                              or spike rates across cells
%       'populationActivity' = line plot of mean spike rate over time, all cell groups
%       'connectionDensities'= shows where connections occur geospatially
%       'STAWRatioHistogram' = bar chart showing distribution of STA/weight ratio
%       'rStdHistogram'      = histogram of reconstructed sample standard deviations
%       'spikeRaster'        = raster plot of spikes over time
%       'spikeRasterAll'     = raster plot of spikes over time for all cell groups
%       'spikeRasterAllImage'= same as spikeRasterAll, but with faster image drawing
%       'spikeRasterAllMulti'= a grid of multiple raster plots
%       'potential'          = timeline plot of the evolving potentials of active cells
%       'networkDynamics'    = plot cell potentials and spike raster for whole network
%       'timelineStats'      = long-term timeline ticker of tracked metrics
%                              (e.g. spikeRate, threshold, ...)
%       'cellGroupVars'      = plot cell group activity metrics
%       'custom'             = a user-defined figure drawn via callback function
%       'composite'          = a figure composed of other figures
%
% Special fields for particular plotTypes:
%    'weightPatches':
%        maxDisplay      = max number of patches to display.
%        combinePlusMinus= subtract second half of weights from first to recreate
%                          receptive fields in factored plus-minus scenario
%    'STA':
%        sourceName      = cells to use for spike trigger (e.g. 'cg_V1e')
%        measureName     = name of an STA measure to use as source for plot
%                          (use only one of sourceName or measureName)
%    'STA_movie':
%        numFrames       = how many frames to display in the movie (default = 5)
%        playDuration    = duration of clip playback in seconds (default = 1)
%        minImageSize    = minimum image size so it isn't too small to see (default = 300)
%    'weightDistribution':
%        sourceName      = which weights to display (e.g. 'cg_V1e.in_input')
%        useLogScale     = if true, use a log scale on the X (weight size) axis
%        minWeight       = if provided, discard all weights with magnitude smaller than minWeight
%        numBins         = if provided, overrides the auto-calculated number of bins
%        plotStyle       = either 'line' or 'bar'. (default = 'line')
%    'weightDistributionAll':
%        [same as 'weightDistribution']
%    'weightCorrelation'
%        sourceName      = two bi-directed input blocks, e.g. {'cg_V1i.in_V1e', 'cg_V1e.in_V1i'}
%    'spikeRateHistogram':
%        dimension       = which spike rate measurements to use, either 'perSample' or 'perCell'
%        numBins         = number of bins to use when drawing histogram (optional). If
%                          numBins = 'tabulate', then perform discrete value tabulation.
%    'STAWRatioHistogram':
%        sourceName      = which weights to display (e.g. 'cg_V1e.in_input')
%        measureName     = which STA measure to use (optional)
%    'spikeRasterAll':
%        bgColor         = background color for MatLab scatter plot (default = [0 0 0] = black)
%        markerSize      = area size to use for spike markers (default = 11)
%    'spikeRasterAllImage':
%        networkId       = which network (sample) to display when for multi-sample
%                          simulations (default = 1)
%    'potential':
%        sourceName      = which cell group to display (e.g. 'cg_V1e')
%        networkId       = which network / sample to plot (default = 1)
%        normalize       = rescale each cell's membrane potential to be relative
%                          to its spike threshold, so that the cell spikes at
%                          potential >= 1. (default = true)
%        maxLines        = how many lines to draw (default = 40)
%        colormap        = colormap to use, from 'help colormap'. 'Lines'
%                          and 'Winter' are good. Can also be an explicit RGB table.
%                          'cellColor' will use cell-by-cell colors in the
%                          cellGroup.cellColor field. 'colorRamp' will generate
%                          ramp based on the cell group's displayColor. absColorRamp
%                          does the same thing, but always assigns a given cell
%                          the same color.
%                          (default = 'Lines')
%        cellIds         = which cells to plot and in what order
%                          (overrides automatic ranking based on cell activity)
%    'networkDynamics':
%        networkId       = which network / sample to plot (default = 1)
%    'timelineStats':
%        metricExpr     = a cell array of MatLab expressions to track
%        colormap       = matrix(N,3) of RGB colors to use for lines.
%                         (default = colormap('Lines'))
%        legendText     = cell array of labels for each line plotted
%        legendLocation = where to place the legend (default = 'SouthWest').
%    'cellGroupVars':
%        varNames       = cell array(N) of variables to plot. The first variable
%                         becomes the dependent variable on X axis, and all
%                         remaining variables are plotted on the Y axis against
%                         the X variable. (required)
%        showPDF        = draw probability distribution?  (default = false)
%        normalize      = normalize Y values to the mean (default = false
%                         unless multiple lines are drawn)
%        colormap       = colormap to use, which is a matrix(N,3) of RGB values.
%                         [] indicates default to colormap('Lines'). (default = [])
%        legendText     = cell array of names for each variable (optional)
%        lineWidth      = line with to use for line plotting (default = 1.5)
%        fontSize       = font size for text (default = 13)
%    'custom':
%        figureFn        = function pointer of form fig = fn(fig, model, spikeHistory)
%                          to plot into a drawing context that has already been established.
%    'composite':
%        subFigure{}     = subfigures (only one level of nesting is allowed)
%        layoutDims      = dimensions [rows,cols] for plot arrangement (optional)
%
% Inputs:
%    fig          = the figure to draw (see subfields described above)
%    model        = the model (read-only)
%    spikeHistory = the spike history (used by spikeRaster plot)
%
% Created:   2/5/2012, Paul King
%--------------------------------------------------------------------------
function fig = NetModel_Plot( fig, model, spikeHistory )

    defaultChartTitle = '';

    switch fig.plotType

    case 'weightPatches'
        % show image-style weight patches
        if ~isfield(fig, 'maxDisplay')
            fig.maxDisplay = 256;
        end
        [i_cg, i_ib] = NetModel_FindElement(model, fig.sourceName);
        W  = model.cellGroup{i_cg}.inputBlock(i_ib).weight;
        if isfield(fig, 'combinePlusMinus') && fig.combinePlusMinus;
            half = size(W,2) / 2;
            W    = W(:,1:half) - W(:,half+1:end);
        end
        if size(W,1) > fig.maxDisplay
            W = W(1:fig.maxDisplay,:);
        end
        if strcmp(model.cellGroup{i_cg}.inputBlock(i_ib).name, 'input')
            inputDims = model.inputDims;
        elseif rem(size(W,2), 1) == 0
            inputDims = [1 1] * sqrt(size(W,2));
        else
            error('can''t display weight patches for non-square number of cells');
        end
        bases = reshape(W', inputDims(1), inputDims(2), []);
        defaultChartTitle = ['Weight Patches (' model.cellGroup{i_cg}.name ')'];
        ShowImagePatchGrid(bases);

    case 'connectionDensities'
        % show 2D connection densities, e.g. from geolocal or uniform partial wiring
        if ~isfield(fig, 'initialized') || ~fig.initialized
            [i_cg, i_ib] = NetModel_FindElement(model, fig.sourceName);
            ib = model.cellGroup{i_cg}.inputBlock(i_ib);
            defaultChartTitle = ['Connection Densities (' ib.name ' -> ' model.cellGroup{i_cg}.name ')'];
            if ~isempty(ib.inputIndices)
                error('connection density for indirectly wired input blocks not yet supported.');
            end
            if strcmp(ib.name, 'input')
                inputDims = model.inputDims;
            elseif rem(ib.numInputs, 1) == 0
                inputDims = [1 1] * sqrt(ib.numInputs);
            else
                error('can''t display connection densities for non-square number of cells');
            end
            connectionRFs = reshape(~ib.clampZero', inputDims(1), inputDims(2), []);
            ShowImagePatchGrid(connectionRFs);
            fig.initialized = true;
        end

    case 'STA'
        % show spike-triggered average for a cell group
        if ~isfield(fig, 'maxDisplay')
            fig.maxDisplay = 256;
        end
        m          = model.stats.measure.(fig.measureName);
        i_cg       = m.cellGroupId;
        numDisplay = min(model.cellGroup{i_cg}.numCells, fig.maxDisplay);
        defaultChartTitle = ['Spike-Triggered Average (' model.cellGroup{i_cg}.name ')'];
        imagePatches = reshape(m.STA(:,1:numDisplay), model.inputDims(1), model.inputDims(2), []);
        ShowImagePatchGrid(imagePatches);

    case 'STA_movie'
        % show spike-triggered movie of time-series input for a cell group
        if ~isfield(fig, 'initialized') || ~fig.initialized
            defaultValues.minImageSize = 300;
            defaultValues.playDuration = .5;    % seconds
            fig = ApplyDefaultValues(fig, defaultValues);
            fig.initialized = true;
        end
        if ~isfield(fig, 'enabled') || fig.enabled
            m            = model.stats.measure.(fig.measureName);
            i_cg         = m.cellGroupId;
            numDisplay   = min(model.cellGroup{i_cg}.numCells, 256);
            defaultChartTitle = ['Spike-Triggered Average Movie (' model.cellGroup{i_cg}.name ')'];
            mov = struct('colormap', [], 'cdata', []);
            for i = 1:m.numFrames
                k = ceil( (i-.5) / m.numFrames * size(m.STA,3) );
                imagePatches = reshape(m.STA(:,1:numDisplay,k), model.inputDims(1), model.inputDims(2), []);
                data         = ShowImagePatchGrid(imagePatches);
                if max(size(data)) < fig.minImageSize
                    data = imresize(data, fig.minImageSize/max(size(data)), 'nearest');
                end
                mov(i).cdata = uint8(data * 255);
            end
            fps = m.numFrames / fig.playDuration;        % calculate playback speed
            PlayMovie(mov, fps);
        end

    case 'inputImagePatches'
        % show input image patches (or movie patches)
        if ~isfield(fig, 'maxDisplay')
            fig.maxDisplay = 256;
        end
        numDisplay = min(model.numSamplesPerBatch, fig.maxDisplay);
        inputData  = model.snapshot.inputData(:,1:numDisplay,:);
        
        if size(inputData, 3) == 1
            % show image patches
            defaultChartTitle = 'input image patches';
            imagePatches      = reshape(inputData, model.inputDims(1), model.inputDims(2), []);
            ShowImagePatchGrid(imagePatches);
        else
            % show movie patches
            if ~isfield(fig, 'initialized') || ~fig.initialized
                defaultValues.minImageSize = 300;
                defaultValues.playDuration = .5;    % seconds
                defaultValues.numFrames    = min(size(inputData, 3), 5);   % up to 5 frames
                fig = ApplyDefaultValues(fig, defaultValues);
                fig.initialized = true;
            end
            defaultChartTitle = 'input movie patches';
            mov = struct('colormap', [], 'cdata', []);
            numIterations = size(inputData, 3);
            for i = 1:fig.numFrames
                k = ceil( (i-.5) / fig.numFrames * numIterations );
                imagePatches = reshape(inputData(:,:,k), model.inputDims(1), model.inputDims(2), []);
                data         = ShowImagePatchGrid(imagePatches);
                if max(size(data)) < fig.minImageSize
                    data = imresize(data, fig.minImageSize/max(size(data)), 'nearest');
                end
                mov(i).cdata = uint8(data * 255);
            end
            fps = fig.numFrames / fig.playDuration;        % calculate playback speed
            PlayMovie(mov, fps);
        end

    case 'weightMatrix'
        % show weight matrix with weight level represented as pixel brightness
        [i_cg, i_ib] = NetModel_FindElement(model, fig.sourceName);
        ib = model.cellGroup{i_cg}.inputBlock(i_ib);
        defaultChartTitle = ['Weight Matrix (' ib.name ' -> ' model.cellGroup{i_cg}.name ')'];
        set(gca,'units','points');    % TODO prevents plot-shrinking bug, but breaks resizing
        imagesc(ib.weight), colorbar, axis image

    case 'weightDistribution'
        % show probability distribution of weights
        if ~isfield(fig, 'useLogScale')
            fig.useLogScale = false;          % default
        end
        [i_cg, i_ib]      = NetModel_FindElement(model, fig.sourceName);
        ib                = model.cellGroup{i_cg}.inputBlock(i_ib);
        W                 = ib.weight;
        defaultChartTitle = ['Weight Distribution (' ib.name ' -> ' model.cellGroup{i_cg}.name ')'];

        % remove disabled weights
        clampZero = ib.clampZero;
        if ~isempty(clampZero)
            keep            = true(size(W));
            keep(clampZero) = false;
            W               = W(keep);
        end

        % remove weights <= 0 if in nonneg mode or using a log scale
        if strcmp(ib.constrainWeights, 'nonneg') || fig.useLogScale
            W = W(W > 0);
        end

        % remove small weights if requested
        if isfield(fig, 'minWeight') && fig.minWeight > 0
            W = W(abs(W) >= fig.minWeight);
        end

        % calculate probability density histogram (possibly using log scale histogram bins)
        numPoints = numel(W);
        if isfield(fig, 'numBins') && fig.numBins > 0
            numBins = fig.numBins;
        else
            numBins = ceil(numPoints^.5);
        end
        if fig.useLogScale
            [numPointsPerBin, logBinLocation] = hist(log(W(:)), numBins);
            binLocation = exp(logBinLocation);
        else
            [numPointsPerBin, binLocation] = hist(W(:), numBins);
        end
        binDensity = numPointsPerBin / numPoints;

        % generate plot, either line graph (default) or bar chart
        if isfield(fig, 'plotStyle') && strcmp(fig.plotStyle, 'bar')
            bar(binLocation, binDensity);
        else
            plot(binLocation, binDensity);
        end
        
        % set log axis if we're using log scale bins
        if fig.useLogScale
            set(gca, 'XScale', 'log');
        end

    case 'weightDistributionAll'
        % show weight distribution for all active weights
        cg = model.cellGroup;
        defaultChartTitle = 'Weight Distribution - ALL';
        
        % initialize: generate a list of all active/enabled input blocks
        if ~isfield(fig, 'initialized') || ~fig.initialized
            clf();
            fig.ibIds = zeros(0,2);
            for i = 1:numel(cg)
                if cg{i}.isExternal || strcmp(cg{i}.cellType, 'disabled')
                    continue;
                end
                ibIds = find( ~strcmp({cg{i}.inputBlock.connectionType}, 'disabled') );
                for k = ibIds
                    fig.ibIds(end+1, :) = [i, k];
                end
            end
            fig.initialized = true;
        end
        
        % draw each input block as a row in a composite figure
        numSubplots = size(fig.ibIds,1);
        for i = 1:numSubplots
            subplot(numSubplots, 1, i);
            subFig = fig;
            subFig.plotType   = 'weightDistribution';
            subFig.title      = [];
            subFig.sourceName = [ 'cg_' cg{fig.ibIds(i,1)}.name '.in_' ...
                    cg{fig.ibIds(i,1)}.inputBlock(fig.ibIds(i,2)).name ];
            subFig = NetModel_Plot(subFig, model, spikeHistory);
            title(subFig.title);
        end

    case 'weightCorrelation'
        % compare a recurrent weight set pair
        [i_cg_1, i_ib_2] = NetModel_FindElement(model, fig.sourceName{1});
        [i_cg_2, i_ib_1] = NetModel_FindElement(model, fig.sourceName{2});
        cg               = model.cellGroup;

        defaultChartTitle = ['Weight Correlation (' cg{i_cg_1}.name ' - ' cg{i_cg_2}.name ')'];

        % make sure weights connect the two cell groups
        assert( cg{i_cg_1}.inputBlock(i_ib_2).srcId == i_cg_2);
        assert( cg{i_cg_2}.inputBlock(i_ib_1).srcId == i_cg_1);

        % get the respective weight blocks
        W_1to2     = cg{i_cg_1}.inputBlock(i_ib_2).weight';
        W_2to1     = cg{i_cg_2}.inputBlock(i_ib_1).weight';
        W_1to2_inv = W_2to1';

        % plot 1->2 weights on X axis against 2->1 weights on Y axis:
        scatter(W_1to2(:), W_1to2_inv(:), 7, 'o');

    case 'spikeRateHistogram'
        % show histogram of spike rates across training samples for a cell group
        m  = model.stats.measure.(fig.measureName);

        % determine which spike rate measurements to use, either 'perSample' or 'perCell'
        if isfield(fig, 'dimension') && strcmp(fig.dimension, 'perCell')
            dimName  = 'Cell';
            srPoints = m.spikeRate;
        else
            dimName  = 'Sample';
            srPoints = m.srPoints(:,1:m.srCount);
            srPoints = srPoints(:);       % reshape to vector
        end
        defaultChartTitle = ['Mean Spike Rate per ' dimName ' (' model.cellGroup{m.cellGroupId}.name ')'];
        
        % compute histogram of mean spike rates across samples
        numPoints = numel(srPoints);
        if isfield(fig, 'numBins') && strcmp(fig.numBins, 'tabulate')
            table           = tabulate(srPoints);
            binLocation     = table(:,1);
            numPointsPerBin = table(:,2);
        else
            if isfield(fig, 'numBins')
                numBins = fig.numBins;
            else
                numBins = ceil(min(numPoints^.4, numPoints^.8/2));
            end
            [numPointsPerBin, binLocation] = hist(srPoints, numBins);
        end
        binDensity = numPointsPerBin / numPoints;
        
        % generate plot, either line graph (default) or bar chart
        if isfield(fig, 'plotStyle') && strcmp(fig.plotStyle, 'bar')
            bar(binLocation, binDensity, 'BarWidth',1, 'EdgeColor','none');
            xlim([0 max(binLocation)]);
        else
            plot(binLocation, binDensity);
        end

    case 'populationActivity'
        % show graph of population activity over time (iterations)
        defaultChartTitle = 'Population Activity (spikes/timeUnit)';

        % first-time variable initialization
        if ~isfield(fig, 'initialized') || ~fig.initialized
            fig.biasedStats_max     = 0;
            fig.biasedMeanSpikeRate = [];
            if ~isfield(fig, 'numInputSamples')
                fig.numInputSamples = 1000;            % default # input samples in moving average
            end
            for j = 1:numel(model.cellGroup)
                fig.biasedMeanSpikeRate{j} = zeros(1, model.numIterationsPerSample);
            end
            fig.initialized = true;
        end

        stats_eta           = 1 - exp(- model.numSamplesPerBatch / fig.numInputSamples);
        fig.biasedStats_max = (1-stats_eta) .* fig.biasedStats_max + stats_eta;

        newplot();
        hold('on');
        for j = 1:numel(model.cellGroup)
            if model.cellGroup{j}.isExternal || strcmp(model.cellGroup{j}.cellType, 'disabled')
                continue;               % skip external and disabled cell groups
            end
            [numCells, numNetworks, numIterations] = size(spikeHistory{j});
            assert(numNetworks == model.numSamplesPerBatch);
            assert(numIterations == model.numIterationsPerSample);

            meanSpikes = reshape(mean(mean(spikeHistory{j}, 1), 2), 1, []) / model.simTimeStep;
            fig.biasedMeanSpikeRate{j} = (1 - stats_eta) * fig.biasedMeanSpikeRate{j} + stats_eta * meanSpikes;
            meanSpikeRate = fig.biasedMeanSpikeRate{j}  / fig.biasedStats_max;
            % xAxisUnits = (1:numIterations) * model.simTimeStep;
            xAxisUnits = 1:numIterations;
            plot(xAxisUnits, meanSpikeRate, 'Color',model.cellGroup{j}.displayColor);
        end
        hold('off');

    case 'STAWRatioHistogram'
        % show probability distribution of std(STA) / std(weight) ratio across cells
        [i_cg, i_ib]      = NetModel_FindElement(model, fig.sourceName);
        cg                = model.cellGroup{i_cg};
        defaultChartTitle = ['STA/W Distribution (' cg.inputBlock(i_ib).name ' -> ' cg.name ')'];
        STA               = model.stats.measure.(fig.measureName).STA;
        W                 = cg.inputBlock(i_ib).weight;

        % calculate ratio coefficients
        if isfield(fig, 'combinePlusMinus') && fig.combinePlusMinus;
            half = size(W,2) / 2;
            W    = W(:,1:half) - W(:,half+1:end);
        end
        coeffs = std(STA, 0, 1) ./ std(W, 0, 2)';

        % plot probability density histogram of coeffs 
        numPoints = numel(coeffs);
        numBins   = ceil(sqrt(numPoints));
        [numPointsPerBin, binLocation] = hist(coeffs, numBins);
        binDensity = numPointsPerBin / numPoints;
        bar(binLocation, binDensity, 'blue');

    case 'rStdHistogram'
        % show histogram of reconstructed sample standard deviations
        m  = model.stats.measure.(fig.measureName);
        cg = model.cellGroup{m.src_i_cg};
        defaultChartTitle = ['Reconstructed Sample Std Dev Histogram (' ...
                cg.inputBlock(m.src_i_ib).name ' -> ' cg.name ')'];

        % plot probability density histogram of standard deviations 
        numPoints = numel(m.stdPoints);
        numBins   = ceil(sqrt(numPoints));
        [numPointsPerBin, binLocation] = hist(m.stdPoints, numBins);
        binDensity = numPointsPerBin / numPoints;
        bar(binLocation, binDensity, 'blue');

    case 'spikeRaster'
        % show a spike raster plot of the spike history (using MatLab scatter plot)
        if ~isfield(fig, 'networkId')
            fig.networkId = 1;             % default is to display the first network / sample
        end
        i_cg = NetModel_FindElement(model, fig.sourceName);
        defaultChartTitle = ['Spike Raster Plot (' model.cellGroup{i_cg}.name ')'];
        if ~model.stats.areFiguresShown
            colormap(gray);
        end
        spikeRecord = reshape(spikeHistory{i_cg}(:,fig.networkId,:), model.cellGroup{i_cg}.numCells, []);
        [r,c] = find(spikeRecord == 1);
        scatter(c,r,11,'diamond', 'filled');   % filled diamond markers
        axis('ij');                            % put the origin in the top-left corner
        axis([1 model.numIterationsPerSample 1 model.cellGroup{i_cg}.numCells]);  % force full axes
        xlabel('time (iteration #)');
        ylabel('cell ID');
        set(gca, 'YGrid', 'on');
        %grid('on');

    case 'spikeRasterAll'
        % show a spike raster plot for all cell groups (using MatLab scatter plot)
        defaultChartTitle = 'Spike Raster Plot (all cell groups)';
        assert(~isempty(spikeHistory), 'spikeHistory is required to geneate spike raster plot');
        if ~isfield(fig, 'networkId')
            fig.networkId = 1;
        end
        if ~isfield(fig, 'bgColor')
            fig.bgColor = 'black';
        end
        if ~isfield(fig, 'markerSize')
            fig.markerSize = 11;
        end
        if ~isfield(fig, 'markerType')
            % 'o' (disc) looks nice but is slow; 'd' (diamond) and 's' (square) are faster
            fig.markerType = 'o';
        end
        if ~isfield(fig, 'fontSize') || isempty(fig.fontSize)
            fig.fontSize = 11;
        end
        if ~isfield(fig, 'activeOnly')
            fig.activeOnly = false;
        end
        rows = [];  cols = []; colorRGB = [];
        totalNumCells = 0;
        cla();
        cg = model.cellGroup;
        for j = 1:numel(cg)
            if cg{j}.isExternal || strcmp(cg{j}.cellType, 'disabled')
                continue;               % skip external and disabled cell groups
            end
            spikeRecord   = reshape(spikeHistory{j}(:,fig.networkId,:), cg{j}.numCells, []);
            if fig.activeOnly
                activeCellIds = find(sum(spikeRecord, 2) > 0);
                spikeRecord   = spikeRecord(activeCellIds,:);
            end
            [r,c]         = find(spikeRecord == 1);
            if isfield(cg{j}, 'cellColor') && ~isempty(cg{j}.cellColor)
                rgb = cg{j}.cellColor(activeCellIds(r),:);
            else
                rgb = ones(size(r)) * cg{j}.displayColor;
            end
            rows          = [rows;     r(:)+totalNumCells];    % offset row as we iterate
            cols          = [cols;     c(:)];
            colorRGB      = [colorRGB; rgb];
            %alt version: plot points as group; does not activate Painter's mode, so looks worse
            %hold on; scatter(c, r+totalNumCells, fig.markerSize, model.cellGroup{j}.displayColor, 'o', 'filled');  hold off;
            totalNumCells = totalNumCells + size(spikeRecord, 1);
        end
        if ~isempty(rows)
            scatter(cols, rows, fig.markerSize, colorRGB, fig.markerType, 'filled');
        end
        axis('ij');                                         % put origin at top-left corner
        axis([1 model.numIterationsPerSample 0 max(totalNumCells+1,5)]);      % force full axes
        xlabel('time (iteration #)', 'FontSize',fig.fontSize);
        if fig.activeOnly
            ylabel('active cells', 'FontSize',fig.fontSize);
        else
            ylabel('cell ID', 'FontSize',fig.fontSize);
        end
        set(gca, 'color',fig.bgColor, 'FontSize',fig.fontSize); % set background color & font
        if totalNumCells >= 10
            set(gca, 'YGrid', 'on');                            % not visible with black background!
        end

    case 'spikeRasterAllImage'
        % show a spike raster plot for all cell groups (using pixelated image)
        if ~isfield(fig, 'networkId')
            fig.networkId = 1;
        end
        defaultChartTitle = 'Spike Raster Plot (all cell groups)';
        
        % generate display image
        rasterImage = [];
        for j = 1:numel(model.cellGroup)
            if model.cellGroup{j}.isExternal || isempty(model.cellGroup{j}.inputBlock) || strcmp(model.cellGroup{j}.cellType, 'disabled')
                continue;               % skip external & disabled cell groups (and delayLines)
            end
            spikeRecord   = reshape(spikeHistory{j}(:,fig.networkId,:), model.cellGroup{j}.numCells, []);
            cgRasterImage = zeros( [size(spikeRecord), 3] );
            for k = 1:3
                cgRasterImage(:,:,k) = (spikeRecord == 1) .* model.cellGroup{j}.displayColor(k);
            end
            rasterImage = cat(1, rasterImage, cgRasterImage);
        end
        
        % draw the spike raster
        image(rasterImage);
        axis('ij');                            % put origin at top-left corner
        xlabel('Time (iteration #)');
        ylabel('Cell ID');

    case 'spikeRasterAllMulti'
        % draw multiple spike raster plots
        if ~isfield(fig, 'layoutDims')
            fig.layoutDims = [1 4];
        end

        % draw each 'spikeRasterAllImage' plot in a subplot array described by layoutDims
        for i = 1:prod(fig.layoutDims)
            h = subplot(fig.layoutDims(1), fig.layoutDims(2), i);
            subFig = fig;
            subFig.plotType   = 'spikeRasterAllImage';
            subFig.title      = [];
            subFig.networkId  = i;
            subFig = NetModel_Plot(subFig, model, spikeHistory);
        end

    case 'potential'
       % plot timeline of the evolving potentials of active cells
       fig = Plot_Potential(fig, model, spikeHistory);

    case 'networkDynamics'
        % show the membrane potential dynamics for several (all?) cell groups
        if ~isfield(fig, 'networkId')
            fig.networkId = 1;
        end
        if ~isfield(fig, 'activeOnly')
            fig.activeOnly = true;
        end
        cg = model.cellGroup;
        defaultChartTitle = 'Network Dynamics';
        
        % initialize
        if ~isfield(fig, 'initialized') || ~fig.initialized
            clf();
            if isfield(fig, 'sourceName')
                fig.srcIds = NetModel_FindElement(model, sourceName);
            else
                % generate a list of all cell groups
                fig.srcIds = [];
                for i = 1:numel(cg)
                    if ~cg{i}.isExternal && ~strcmp(cg{i}.cellType, 'disabled')
                        fig.srcIds(end+1) = i;
                    end
                end
            end
            fig.initialized = true;
        end
        
        % plot potentials for each cell group in a subplot
        numCGs = numel(fig.srcIds);
        for i = 1:numCGs
            ii = fig.srcIds(i);
            subplot(numCGs+1, 1, i);
            subFig            = fig;                  % inherit attributes
            subFig.plotType   = 'potential';
            subFig.sourceName = [ 'cg_' cg{ii}.name ];
            subFig.title      = [];
            subFig.baseTitle  = [cg{ii}.name ' cells'];
            if isfield(cg{ii}, 'cellColor')
                subFig.colormap = 'cellColor';
            else
                subFig.colormap   = 'absColorRamp';    % either 'colorRamp' or 'absColorRamp'
            end
            subFig = Plot_Potential(subFig, model, spikeHistory);
            title(subFig.title);

            % assign plotted colors temporarily to cellGroup
            if ~isfield(model.cellGroup{ii}, 'cellColor')
                model.cellGroup{ii}.cellColor = zeros(model.cellGroup{ii}.numCells,3);
                for k = 1:numel(subFig.plottedCellIds)
                    cellId = subFig.plottedCellIds(k);
                    model.cellGroup{ii}.cellColor(cellId,:) = subFig.plottedCellColors(k,:);
                end
            end
        end

        % plot the raster summary in the last subplot
        subplot(numCGs+1, 1, numCGs+1);
        rasterFig            = fig;                  % inherit attributes
        rasterFig.plotType   = 'spikeRasterAll';
        rasterFig.networkId  = fig.networkId;
        rasterFig.activeOnly = fig.activeOnly;
        rasterFig.title      = '';
        rasterFig.bgColor    = [1 1 1];
        rasterFig = NetModel_Plot(rasterFig, model, spikeHistory);

    case 'timelineStats'
       fig = Plot_TimelineStats(fig, model, spikeHistory);

    case 'cellGroupVars'
       fig = Plot_CellGroupVars(fig, model, spikeHistory);

    case 'composite'
        % show a composite figure containing subfigures
        if ~isfield(fig, 'layoutDims')
            numSubplots         = numel(fig.subFigure);
            fig.layoutDims(1,2) = ceil(sqrt(numSubplots));
            fig.layoutDims(1,1) = ceil(numSubplots / fig.layoutDims(2));
        end
        for i = 1:numel(fig.subFigure)
            if ~isempty(fig.subFigure{i})
                subplot(fig.layoutDims(1), fig.layoutDims(2), i);
                fig.subFigure{i} = NetModel_Plot(fig.subFigure{i}, model, spikeHistory);
                title(fig.subFigure{i}.title);
            end
        end

    case 'custom'
        % execute custom plot drawing function
        fig = feval(fig.figureFn, fig, model, spikeHistory);

    otherwise
        error('Unknown figure plot type "%s"', fig.plotType);
    end


    % set the figure title if it has not been set yet
    if ~isfield(fig, 'title') || ~ischar(fig.title)
        fig.title = defaultChartTitle;
    end
end


% Plot_Potential
%
% Plot the (membrane) potentials for a particular cell group over time.
% The most "important" cells will be plotted, according to a heuristic
% based on variation in the potential and number of spikes emitted.
%
% Inputs:
%    fig          = the figure descriptor (model.stats.figure{i})
%        sourceName   = which cell group to display (e.g. 'cg_V1e')
%        networkId    = which network to plot (default = 1)
%        baseTitle    = root title to which stats will be added (preferred method)
%        title        = title to display (overrides baseTitle)
%        colormap     = colormap to use, from 'help colormap'. 'Lines'
%                       and 'Winter' are good. Can also be an explicit RGB table.
%                       (default = 'Lines')
%        cellIds      = which cellIds to plot and in what order
%                       (overrides automatic ranking based on cell activity)
%        activeOnly   = only display active cells
%        maxLines     = maximum number of lines to draw (default = 40)
%        normalize    = rescale each cell's membrane potential to be relative
%                       to its spike threshold, so that the cell spikes at
%                       potential >= 1. (default = true)
%        overlay      = if true, draw on top of existing plot (default = false)
%        ylabel       = alternate label to use for y axis (optional)
%
% Outputs:
%    fig
%        plottedCellIds
%        plottedCellColors
%--------------------------------------------------------------------------
function fig = Plot_Potential( fig, model, spikeHistory )

    % initialize figure parameters
    i_cg = NetModel_FindElement(model, fig.sourceName);
    if ~isfield(fig, 'networkId')
        fig.networkId = 1;             % default is to display the first network / sample
    end
    if ~isfield(fig, 'maxLines')
        fig.maxLines = 40;             % max number of lines to draw
    end
    if ~isfield(fig, 'colormap')
        fig.colormap = 'Lines';        % colors to use for lines (see 'help colormap')
    end
    if ~isfield(fig, 'normalize')
        fig.normalize = true;          % normalize membrane potential to spike threshold
    end
    if ~isfield(fig, 'YMin')
        fig.YMin = [];
    end
    if ~isfield(fig, 'lineWidth') || isempty(fig.lineWidth)
        fig.lineWidth = 1;
    end
    if ~isfield(fig, 'fontSize') || isempty(fig.fontSize)
        fig.fontSize = 11;
    end
    if ~isfield(fig, 'title') || ~ischar(fig.title)
        if isfield(fig, 'baseTitle') && ~isempty(fig.baseTitle)
            fig.title = fig.baseTitle;
        else
            fig.title = [model.cellGroup{i_cg}.name ' Membrane Potentials (selected cells)'];
        end
    end

    % initialize variables
    if ~model.stats.keepPotentialHistory
        cla(); set(gca,'Visible','off');
        text(.2, .4, 'model.stats.keepPotentialHistory must be turned on first');
        error('model.stats.stats.keepPotentialHistory must be turned on first');
    end
    spikeRecord       = reshape(spikeHistory{i_cg}(:,fig.networkId,:), ...
                            model.cellGroup{i_cg}.numCells, []);
    potentialHist     = reshape(model.snapshot.potentialHistory{i_cg}(:,fig.networkId,:), ...
                            model.cellGroup{i_cg}.numCells, []);
    spikeThresh       = model.cellGroup{i_cg}.spikeThresh;
    potentialHistNorm = bsxfun(@times, potentialHist, 1 ./ spikeThresh);
    potentialHistNorm(spikeRecord > 0) = 2;

    % select which cells to plot and in what order
    if isfield(fig, 'cellIds') && ~isempty(fig.cellIds)
        cellIds = fig.cellIds;            % fig.cellIds overrides our calculations
    else
        % calculate an importance score and sort based on that
        %importanceScore = var(potentialHist, 1, 2);            % potential variance
        importanceScore  = var(max(potentialHist,0), 1, 2);     % rectified potential variance
        [~, cellIds] = sort(importanceScore, 'descend');
    end
    if isfield(fig, 'activeOnly') && fig.activeOnly
        % filter out any non-active cells
        activeCellIds = find(sum(spikeRecord, 2) > 0);
        cellIds       = cellIds(ismember(cellIds, activeCellIds));
        if isfield(fig, 'baseTitle') && ~isempty(fig.baseTitle)
            fig.title = sprintf('%s (%d active, %d spikes)', fig.baseTitle, ...
                    numel(activeCellIds), sum(spikeRecord(:)));
        end
    end
    numPlot  = min(fig.maxLines, numel(cellIds));
    cellIds  = cellIds(1:numPlot);

    % select color to use for each line
    if strcmp(fig.colormap, 'cellColor')
        if isfield(model.cellGroup{i_cg}, 'cellColor') && ~isempty(model.cellGroup{i_cg}.cellColor)
            linecolor = model.cellGroup{i_cg}.cellColor(cellIds,:);
        else
            error('colormap = ''cellColor'' requires cellColor to be set in cellGroup');
        end
    elseif strcmp(fig.colormap, 'colorRamp') || strcmp(fig.colormap, 'absColorRamp')
        % generate a color ramp based on the cell's display color
        if strcmp(fig.colormap, 'absColorRamp')
            rampVal = (cellIds - .5) / model.cellGroup{i_cg}.numCells;
        else
            rampVal = (.5:numPlot) / numPlot;
        end
        hsv          = rgb2hsv(model.cellGroup{i_cg}.displayColor);
        hue          = hsv(1);
        cellHSV      = zeros(numPlot,3);
        cellHSV(:,1) = hue;
        cellHSV(:,2) = min(.3 + (1-rampVal), 1);
        cellHSV(:,3) = min(.4 + .8*rampVal, 1);
        linecolor    = hsv2rgb(cellHSV);
        
    elseif ischar(fig.colormap)
        linecolor = colormap(fig.colormap);       % look up color map name
    elseif isnumeric(fig.colormap) && size(fig.colormap, 2) == 3
        linecolor = fig.colormap;                 % use color map explicitly
    else
        error('unknown type for property "colormap"');
    end
    if ~ strcmpi(fig.colormap, 'Lines') && numel(cellIds) ~= size(linecolor,1)
        % for gradated colormaps, space lines across whole color space
        colorIds  = ceil((1:numel(cellIds)) / numel(cellIds) * size(linecolor,1));
        linecolor = linecolor(colorIds, :);
    end

    % draw plot
    if ~isfield(fig, 'overlay') || ~fig.overlay
        cla();              % erase previous lines
    end
    X = 1:model.numIterationsPerSample;
    Y = zeros(1, model.numIterationsPerSample);
    for i = 1:numel(cellIds)
        k = cellIds(i);
        if fig.normalize
            Y(:) = potentialHistNorm(k,:);
        else
            Y(:) = potentialHist(k,:);
            Y(spikeHistory(k,:) > 0) = max(spikeThresh) * 2;
        end
        hold on; plot(X, Y, 'Color',linecolor(i,:), 'LineWidth',fig.lineWidth); hold off;
    end
    set(gca, 'box','off', 'FontSize',fig.fontSize);
    xlabel('time (iteration #)','FontSize',fig.fontSize);
    if isfield(fig, 'ylabel')
        ylabel(fig.ylabel,'FontSize',fig.fontSize);
    elseif fig.normalize
        ylabel('norm potential','FontSize',fig.fontSize);
    else
        ylabel('potential','FontSize',11);
    end
    xlim([0 model.numIterationsPerSample]);       % force X axis, even if no spikes
    if fig.normalize
        if ~isempty(fig.YMin)
            YMin = fig.YMin;
        elseif ~isempty(cellIds)
            YMin = min(min(potentialHistNorm(cellIds,:)));
        else
            YMin = -2;
        end
        ylim([min(YMin,-2) 2]);                   % prevent awkward Y axis auto-scaling
    end
    
    % save plot info in figure for possible future reference
    fig.plottedCellIds    = cellIds;
    fig.plottedCellColors = linecolor;
end


% Plot_TimelineStats
%
% Plot a timeline ticker showing the evolution of various metrics
% up to the present.
%
% Inputs:
%    fig          = the figure descriptor (model.stats.figure{i})
%        metricExpr     = a cell array of MatLab expressions for each metric
%                         to track.  See NetModel_EvalExpression for details
%                         on automatic variables in the evaluation scope.
%        colormap       = colormap to use, which is a matrix(N,3) of RGB values
%                         [] indicates default to colormap('Lines'). (default = [])
%        legendText     = cell array of labels for each line plotted
%        legendLocation = where to place the legend (default = 'SouthWest').
%                         (see 'help legend')
%--------------------------------------------------------------------------
function fig = Plot_TimelineStats( fig, model, spikeHistory )

    % get measure containing the historical data
    if ~isfield(fig, 'measureName')
        fig.measureName = 'timelineStats';
    end
    m = model.stats.measure.(fig.measureName);
    
    if ~isfield(fig, 'title') || ~ischar(fig.title)
        fig.title = 'Network Statistics Over Time';
    end
    if isfield(fig, 'colormap')
        colors = fig.colormap;
    else
        colors = colormap('Lines');
    end
    
    % generate plot
    cla();
    if m.bufferWrapped,
        datapoints = [m.currentPtr+1:m.historySize, 1:m.currentPtr];
    else
        datapoints = 1:m.currentPtr;
    end
    % datapoints = datapoints(end:-1:1);   % TODO reverse so plot is L-to-R
    X = m.sampleNum(datapoints);
    for i = 1:numel(m.metricExpr)
        sd = std(m.history(i,datapoints));
        Y = m.history(i,datapoints) / sd;
        if sd == 0
            Y(:) = sign(Y(1));     % if all values are the same, normalize to -1, 0, or 1
        end
        hold on; plot(X, Y, 'Color', colors(i,:));
    end
    
    % add legend
    if isfield(fig, 'legendText')
        if isfield(fig, 'legendLocation')
            legend(fig.legendText, 'Location',fig.legendLocation);
        else
            legend(fig.legendText, 'Location','SouthWest');
        end
    end
end


% Plot_CellGroupVars
%
% Plot activity variables related to cells in a particular cell group.
% The variables that can be plotted are the spike rate, and the mean
% input current from different input channels. Any variable can be
% selected as the dependent variable for the X axis, and other variables
% can be plotted on the Y axis. A probability distribution representing
% the occurance frequency for the variable along the X axis can also
% be plotted. Variable names must be an input block, e.g. 'ib_xxx'
% or the string 'spikeRate'.
%
% Inputs:
%    fig        = the figure descriptor (model.stats.figure{i})
%        varNames    = cell array(N) of variables to plot. The first variable
%                      becomes the dependent variable on X axis, and all
%                      remaining variables are plotted on the Y axis against
%                      the X variable. (required)
%        showPDF     = draw probability distribution? If this is a scalar value,
%                      scale the PDF by that value for easier viewing. (default = false)
%        normalize   = normalize Y values to the mean (default = false
%                      unless multiple lines are drawn)
%        colormap    = colormap to use, which is a matrix(N,3) of RGB values.
%                      [] indicates default to colormap('Lines'). (default = [])
%        legendText  = cell array of names for each variable (optional)
%        lineWidth   = line with to use for line plotting (default = 1.5)
%        fontSize    = font size for text (default = 13)
% Outputs:
%    fig         = updated figure descriptor
%        varInfo     = struct array of the variables used in the plot, including
%                      recently plotted values.
%        pdf         = array(1,numBins) of the probability density at each bin
%--------------------------------------------------------------------------
function fig = Plot_CellGroupVars( fig, model, spikeHistory )

    if ~isfield(fig, 'initialized') || ~fig.initialized
        fig.srcId = NetModel_FindElement(model, fig.sourceName);

        if ~isfield(fig, 'showPDF') || isempty(fig.showPDF)
            fig.showPDF = true;
        end
        if ~isfield(fig, 'lineWidth') || isempty(fig.lineWidth)
            fig.lineWidth = 1.5;
        end
        if ~isfield(fig, 'fontSize') || isempty(fig.fontSize)
            fig.fontSize = 11;
        end

        % look up variable names in model
        fig.varInfo = [];
        cg = model.cellGroup{fig.srcId};
        for k = 1:numel(fig.varNames)
            name = fig.varNames{k};
            fig.varInfo(k).sourceName = name;      % TODO not needed?
            if strcmp(name, 'spikeRate')
                fig.varInfo(k).ibId   = 0;
                fig.varInfo(k).targetSpikeRate = cg.targetSpikeRate;
                fig.varInfo(k).legend = [cg.name ' spike rate'];
            elseif numel(name) > 3 && strcmp(name(1:3), 'ib_')
                fig.varInfo(k).ibId = FindByName(model.cellGroup{fig.srcId}.inputBlock, name(4:end));
                if isempty(fig.varInfo(k).ibId)
                    error('input block "%s" not found', name);
                end
                fig.varInfo(k).targetSpikeRate = model.cellGroup{fig.srcId}.targetSpikeRate;
                fig.varInfo(k).legend = [name(4:end) '\rightarrow' cg.name ' input'];
            else
                error(['unknown cell variable "' name '"']);
            end
        end
        if ~isfield(fig, 'title') || ~ischar(fig.title)
            fig.title = [cg.name ' cell activity metrics'];
        end
        if ~isfield(fig, 'normalize') || isempty(fig.normalize)
            fig.normalize = (numel(fig.varNames) > 2);
        end
        fig.initialized = true;
    end

    assert(~isempty(spikeHistory), 'spikeHistory is required to geneate cell group variables plot');

    % compute all variables for the average cell in the cell group
    vals = [];
    for k = 1:numel(fig.varInfo)
        if fig.varInfo(k).ibId == 0
            vals{k} = mean(spikeHistory{fig.srcId}, 3) / model.simTimeStep;
            continue;
        end
        ib = model.cellGroup{fig.srcId}.inputBlock(fig.varInfo(k).ibId);
        W  = ib.weight;
        bw = ib.blockWeight;
        if strcmp(model.cellGroup{ib.srcId}.name, 'input')
            netInput = W * mean(model.snapshot.inputData, 3);
        else
            netInput = W * mean(spikeHistory{ib.srcId}, 3);
        end
        vals{k} = bw * netInput;
    end

    % bin according to the dependent variable, varInfo(1)
    if ~isfield(fig, 'numBins') || isempty(fig.numBins) || fig.numBins == 0
        numBins = min(50, sqrt(numel(vals{1})));
    else
        numBins = fig.numBins;
    end
    [~, sortedIdx] = sort(vals{1}(:));
    for k = 1:numel(fig.varInfo)
        fig.varInfo(k).val = zeros(1,numBins);
    end
    lastEnd  = 0;
    N        = numel(sortedIdx);
    for i = 1:numBins
        newEnd     = round(N*i/numBins);
        dataPoints = sortedIdx(lastEnd+1:newEnd);
        for k = 1:numel(fig.varInfo)
            fig.varInfo(k).val(i) = mean(vals{k}(dataPoints));
        end
        lastEnd = newEnd;
    end

    % calculate occurance histogram
    % (normalize so that the area under the pdf curve sums to 1)
    pdfHeight = zeros(1,numBins);
    XVar = vals{1};                 % TODO name?
    X    = fig.varInfo(1).val;
    for i = 2:numBins-1
        pdfHeight(i) = 1 / (X(i+1) - X(i-1));
    end
    pdfHeight(1)       = .5 / (XVar(sortedIdx(round(N/numBins))) - XVar(sortedIdx(1)));
    pdfHeight(numBins) = .5 / (XVar(sortedIdx(N)) - XVar(sortedIdx(N - round(N/numBins))));
    pdfHeight(~isfinite(pdfHeight)) = NaN;
    fig.pdf            = pdfHeight * (fig.varInfo(1).val(end) - fig.varInfo(1).val(1)) / numBins;

    % prepare plotting properties
    LW      = fig.lineWidth;                % line width for data lines
    LW_mean = min(fig.lineWidth, 1.5);      % line width for average marker
    if isfield(fig, 'colormap')
        colors = fig.colormap;
    else
        colors = colormap('Lines');
        colors = colors(2:end,:);        % skip blue because it is used for the PDF
    end
    
    % generate plot
    cla();
    X = fig.varInfo(1).val;
    for k = 2:numel(fig.varInfo)
        Y = fig.varInfo(k).val;
        if fig.normalize
            Y = Y / mean(abs(Y(:)));
        end
        hold on;
        plot(X, Y, '-', 'LineWidth',LW, 'MarkerSize',3, 'Color',colors(k-1,:));
        hold off;
    end
    xlim([min(X), max(X)]);          % force range to not exceed the points plotted
    set(gca, 'box','off', 'FontSize',fig.fontSize, 'YTick',[0]);
    if isfield(fig, 'legendText') && ~isempty(fig.legendText)
        legendText = fig.legendText;
    else
        legendText = { fig.varInfo.legend };
    end
    if fig.showPDF
        % add occurance frequency histogram
        ylimSave = ylim();
        hold on;  plot(X(isfinite(fig.pdf)), fig.pdf(isfinite(fig.pdf))*fig.showPDF, ...
                '--b', 'LineWidth',LW);  hold off;
        ylim(ylimSave);    % prevent large values from disrupting the Y axis scale
        legend([legendText(2:end), 'occurrence frequency'], 'Location','NorthWest');
    elseif numel(fig.varInfo) > 2
        legend(legendText(2:end), 'Location','NorthWest');
    end
    if fig.normalize
        % plot Y = normalized mean value = 1
        hold on; plot(xlim(), [1 1], '--k', 'LineWidth',LW_mean);  hold off;
    elseif numel(fig.varInfo) == 2 && fig.varInfo(2).ibId == 0
        % plot Y = target spike rate
        hold on; plot(xlim(), [1 1]*fig.varInfo(2).targetSpikeRate, '--k', 'LineWidth',LW_mean);  hold off;
    end
    hold on; plot(xlim(), [0 0], '-k', 'LineWidth',LW/3);  hold off;  % plot X axis (Y = 0)
    xlabel(legendText{1},'FontSize',fig.fontSize);
    if fig.normalize
        ylabel('normalized value','FontSize',fig.fontSize);
    elseif numel(fig.varInfo) == 2
        ylabel(legendText{2},'FontSize',fig.fontSize);
    end
end
