% train3_EICount - Train E vs. I cell count data for surface plot in Figure 4
%
% Run a matrix of training simulations to produce the 3D surface plot
% of Figure 4. 7 values for the E cell count and 7 values for the I
% cell count are used, for a total of 49 full training simulations.
% Training is brief (100,000 input samples) for speed.
%
% Reference:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.


% =================  RUN A MATRIX OF EXPERIMENTS: E vs. I CELL COUNT  =================

% run a matrix of training experiments on networks with different numbers of E and I cells.
% The results of each experiment are collected in a MatLab struct "results" and saved
% to a file. The x 7 x matrix of experiments (7 E cell counts times 7 I cell counts)
% results in 49 total simulations. The results produce the surface plots in Figure 4.
% (total time to run training simulations = 150 minutes)

% Figure 4 (J Neuroscience 2013): surface plot of ecorr & resErr for varying # E and I cells

% Initialize MatLab configuration
global neurosim_initialized; if isempty(neurosim_initialized), init_neurosim; end


% base E-I network to vary (no figures)
simParams = struct();
simParams.numInputSamples              = 100000;
simParams.model.modelType              = 'V1eV1i';
simParams.model.modelSubtype           = 'jneuroscience';
simParams.model.inputDims              = [10,10];        % needed to calc numPixels, below

% add I correlation tracking
simParams.model.stats.measure.corr_V1i.measureType = 'correlation';
simParams.model.stats.measure.corr_V1i.sourceName  = 'cg_V1i';

% run simulation multiple times to find optimal E/I ratio
overcomplete_E = [1, 1.5, 2, 3, 4, 5, 6];
overcomplete_I = [.05 .1, .2, .3, .5, .75, 1];
numPixels      = prod(simParams.model.inputDims);
experimentParams = struct();
experimentParams.propNames         = { 'model.cg_V1e.numCells', 'model.cg_V1i.numCells' };
experimentParams.propValues        = { round(overcomplete_E*numPixels), round(overcomplete_I*numPixels) };
experimentParams.captureSparseness = {'V1e', 'V1i'};
% EXPERIMENT BWvary: adjust blockweight to counteract effect of adding E cells
%    experimentParams.actionPre = 'simParams.model.cg_V1i.in_V1e.blockWeight = 600/simParams.model.cg_V1e.numCells';
experimentParams.actionPost        = 'results.finalCorrI(i,j) = model.stats.measure.corr_V1i.rmsCorrelation';
results = RunParameterizedExperiment(simParams, experimentParams);
results.overcomplete_E = overcomplete_E;
results.overcomplete_I = overcomplete_I;
results.numPixels      = numPixels;
save experiment3_EICount results

