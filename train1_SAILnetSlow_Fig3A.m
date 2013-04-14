% train1_SAILnetSlow - Train SAILnet slowly, for plot in Figure 3A
%
% Train a network using the SAILnet algorithm (400 neurons). Training
% uses a very low (slow) learning rate in order to prevent the network
% from becoming unstable during learning.
%
% Reference:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.


% Initialize MatLab configuration
global neurosim_initialized; if isempty(neurosim_initialized), init_neurosim; end

% train a SAILnet network (190 minutes total)
% (a long training with slow learning is used to prevent network instability during training)
simParams = struct();
simParams.numInputSamples              = 5000000;
simParams.model.modelType              = 'V1eV1i';
simParams.model.modelSubtype           = 'SAILnet';
simParams.model.inputDims              = [10,10];
simParams.model.numIterationsPerSample = 50;
simParams.model.learningRate           = .1;
simParams.model.autoFigures            = {'STA_E'};
simParams.model.cg_V1e.numCells        = 400;
simParams.model.stats.numSTASamples    = 5000;            % average more samples when computing STA
simParams.model.stats.printStatsEvery  = 10000;           % log stats, update figures every 10000 samples
randreset(); model = RunVisionNetworkSimulation(simParams); % (100 minutes)
model.learningRate                     = .02;             % lower rate to smooth RFs
model = RunVisionNetworkSimulation({'numInputSamples',1000000}, model);   % (20 minutes)
model.learningRate                     = 0;               % compute STAs with learning frozen
model = RunVisionNetworkSimulation({'numInputSamples',10000}, model);
save experiment1_SAILnet model

