% demo3_TrainEINetFast - Train EI Net quickly using shortcuts
%
% Train E-I Net using a fast learning rate in order to train quickly.
% This takes around 20 minutes to train and produces nearly identical
% results to the slow training in train2_EINetSlow_Fig3B.
%
% Reference:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.


% Initialize MatLab configuration
global neurosim_initialized; if isempty(neurosim_initialized), init_neurosim; end


% train an E-I Net network using faster learning rate (20 minutes)
% (training the network at higher speed achieves nearly identical results)
simParams = struct();
simParams.numInputSamples              = 500000;
simParams.model.modelType              = 'V1eV1i';
simParams.model.modelSubtype           = 'jneuroscience';
simParams.model.learningRate           = .4;
simParams.model.autoFigures            = {'STA_E', 'STA_I'};
simParams.model.stats.numSTASamples    = 5000;            % average 5000 samples when computing STA
randreset(); model = RunVisionNetworkSimulation(simParams); % (15 minutes)
model.learningRate                     = .1;
model = RunVisionNetworkSimulation({'numInputSamples',100000}, model);  % (3 minutes)
model.learningRate                     = .02;
model = RunVisionNetworkSimulation({'numInputSamples',50000}, model);   % (1 minute)
model.learningRate                     = 0;
model = RunVisionNetworkSimulation({'numInputSamples',10000}, model);
save experiment2_EINetFast model
