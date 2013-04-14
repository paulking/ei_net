% demo1_SmallNetwork - Run a small E-I Net network that trains quickly
%
% Train a small E-I Net network at fast rate (15 seconds). The network has
% 121 E cells and 36 I cells (square numbers for better display layout).
% Smaller-sized 8x8 image patches are used for speed. Several real-time
% plots are displayed to show the progress of learning.
%
% Reference:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.


% Initialize MatLab configuration
global neurosim_initialized; if isempty(neurosim_initialized), init_neurosim; end


% demo: fast training with figure display (30 seconds)
% (this is fast due to smaller 8x8 image patches, fewer neurons, and meanRateLearning=true)
simParams = struct();
simParams.numInputSamples              = 50000;
simParams.model.modelType              = 'V1eV1i';
simParams.model.inputDims              = [8,8];
simParams.model.cg_V1e.numCells        = 121;
simParams.model.cg_V1i.numCells        = 36;
simParams.model.meanRateLearning       = true;       % speeds up training by 2x, results mostly the same
simParams.model.autoFigures            = {'STA_E', 'STA_I', 'weightDist', 'spikeRaster', 'networkDynamics'};
simParams.model.stats.printStatsEvery  = 2000;
randreset(); model = RunVisionNetworkSimulation(simParams);

