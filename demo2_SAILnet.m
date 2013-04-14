% demo4_SAILnet - Run a small SAILnet network that trains quickly
%
% Train a small E-I Network at fast rate (15 seconds). The network has
% 121 E cells and 36 I cells (square numbers for better display layout).
% Smaller-sized 8x8 image patches are used for speed. Several real-time
% plots are displayed to show the progress of learning.
%
% Reference:
%    Joel Zylberberg, Jason Timothy Murphy, Michael Robert DeWeese (2011).
%        Beyond receptive fields: consequences of learning a sparse code
%        for natural images with spiking neurons and synaptically local
%        plasticity rules. PLoS Computational Biology.
%
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.


% Initialize MatLab configuration
global neurosim_initialized; if isempty(neurosim_initialized), init_neurosim; end


% demo: fast training with figure display (20 seconds)
simParams = struct();
simParams.numInputSamples              = 50000;
simParams.model.modelType              = 'V1eV1i';
simParams.model.modelSubtype           = 'SAILnet';
simParams.model.cg_V1e.numCells        = 121;
simParams.model.autoFigures            = {'STA_E', 'weightDist', 'networkDynamics'};
simParams.model.stats.printStatsEvery  = 2000;
randreset(); model = RunVisionNetworkSimulation(simParams);

