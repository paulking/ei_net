% demo2_ReferenceEINet - Run the full EI Net with short training duration and display
%
% Train the full E-I Net quickly as a demonstration (3 minutes). The network
% has 400 E cells and 49 I cells.
%
% Reference:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.


% Initialize MatLab configuration
global neurosim_initialized; if isempty(neurosim_initialized), init_neurosim; end


% demo: training with figure display (3 minutes)
simParams = struct();
simParams.numInputSamples              = 100000;
simParams.model.modelType              = 'V1eV1i';
simParams.model.modelSubtype           = 'jneuroscience';
simParams.model.autoFigures            = {'STA_E', 'STA_I', 'weightDist', 'spikeRaster', 'networkDynamics'};
simParams.model.stats.printStatsEvery  = 2000;
randreset(); model = RunVisionNetworkSimulation(simParams);
