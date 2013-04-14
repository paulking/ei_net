% train4_IFastSpiking - Train Fast Spiking data for surface plot in Figure 7
%
% Run a matrix of training simulations to produce the 3D surface plot
% of Figure 7. 8 values for the I cell spike rate and 6 values
% for the I connection multiplier are used, for a total of 48 full
% training simulations. Training is brief (100,000 input samples)
% for speed.
%
% Reference:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.


% ======  RUN A MATRIX OF EXPERIMENTS: I CELL FAST-SPIKING vs. CONNECTION WEIGHT MULTIPLIER  =======

% run a matrix of training experiments on networks with varying I cell spike rates
% (cg_V1i.targetSpikeRate) and varying I->E connection weight multipliers
% (cg_V1e.in_V1i.blockWeight). The results of each experiment are collected
% in a MatLab struct "results" and saved to a file. The 8 x 6 matrix of experiments
% (8 spike rates, 6 connection weight multipliers) results in 48 total simulations.
% The results produce the surface plots in figure 7.
% (total time to run training simulations = 130 minutes)

% Figure 7 (J Neuroscience): surface plot of recorr and resErr for FS I cells

% Initialize MatLab configuration
global neurosim_initialized; if isempty(neurosim_initialized), init_neurosim; end

% base E-I network to vary (no figures)
simParams = struct();
simParams.numInputSamples              = 100000;
simParams.model.modelType              = 'V1eV1i';
simParams.model.modelSubtype           = 'jneuroscience';

% run simulation multiple times to find optimal blockWeight to spikeRate ratio
% 6*8 = 48 runs (190 minutes)
spikeRate    = [.5  1  1.5  2  2.5  3  3.5  4] * .02;     % multiple of E cell default spike rate
blockWeight  = [1 2 3 4 5 6];
experimentParams = struct();
experimentParams.propNames   = { 'model.cg_V1i.targetSpikeRate', 'model.cg_V1e.in_V1i.blockWeight' };
experimentParams.propValues  = { spikeRate, blockWeight };
results = RunParameterizedExperiment(simParams, experimentParams);
results.spikeRate   = spikeRate;
results.blockWeight = blockWeight;
save experiment4_FSBW results
