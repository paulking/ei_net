% train2_EINetSlow - Train EI_Net slowly, for plot in Figure 3B
%
% Train E-I Net using a slow learning rate for maximally stable
% network results. The results of this training appear in Figure 3B.
%
% Reference:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.


% Initialize MatLab configuration
global neurosim_initialized; if isempty(neurosim_initialized), init_neurosim; end

% train an E-I Net network (110 minutes total)
simParams = struct();
simParams.numInputSamples              = 2000000;
simParams.model.modelType              = 'V1eV1i';
simParams.model.modelSubtype           = 'jneuroscience';
simParams.model.learningRate           = .1;
simParams.model.autoFigures            = {'STA_E', 'STA_I'};
simParams.model.stats.numSTASamples    = 5000;            % average more samples when computing STA
randreset(); model = RunVisionNetworkSimulation(simParams); % (60 minutes)
model.learningRate                     = .02;             % lower rate to smooth RFs
model = RunVisionNetworkSimulation({'numInputSamples',1000000}, model);   % (25 minutes)
model.learningRate                     = 0;               % compute STAs with learning frozen
model = RunVisionNetworkSimulation({'numInputSamples',10000}, model);
save experiment2_EINet model
