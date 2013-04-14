% draw_figures_extra - Draw extra figures not included in the paper
%
% Before this script is run, certain "train" scripts must first be run
%
% Reference:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.



% ==============  SPARSENESS SURFACE PLOT FOR E vs. I CELL COUNT  ==============

% These surface plots show sparseness for different numbers of E and I cells
% This is based on the same simulation data used for Figure 4

load experiment3_EICount

% plot sparseness (figures not included in paper)
figure(44);  surf(X, Y, results.sparse_V1e.time);     set(44, 'Name','E-cell Lifetime Sparseness');
figure(45);  surf(X, Y, results.sparse_V1e.pop);      set(45, 'Name','E-cell Population Sparseness');
figure(46);  surf(X, Y, results.sparse_V1i.time);     set(46, 'Name','I-cell Lifetime Sparseness');
figure(47);  surf(X, Y, results.sparse_V1i.pop);      set(47, 'Name','I-cell Population Sparseness');



% ===========  SURFACE PLOT OF I-TRIGGERED E SPIKERATE  ================

% this 3D surface plot is a generalization of Figure 6D which appears in the paper

% load model from thoroughly trained, stabilized E-I Net used in J Neuroscience
load experiment2_EINet

% run model on 2000 samples with disabled learning
model.numSamplesPerBatch     = 2000;
model.learningRate           = 0;
model.stats.keepSpikeHistory = true;
model = RunVisionNetworkSimulation({'numInputSamples',model.numSamplesPerBatch}, model);

% run E-I Net PSTH analysis
psthParams = struct();
psthParams.numStrongDeciles = 2;
results = AnalyzeV1eV1iPSTH(model, psthParams);


% figure: full 3D surface plot of I triggered E spikerate histogram
[X,Y] = meshgrid(results.timeBinVal, results.strengthBinVal);
figure(65);  set(65, 'Name','I spike-triggered E histogram');
  surf(X, Y, results.histMean);
  xlabel('\Delta time (simulation time steps)','FontSize',11);
  ylabel('I->E connection strength','FontSize',11);
  zlabel('E cell spike rate','FontSize',11);
  set(gca,'XLim',[-15 15], 'XTick',-10:5:10, 'YLim',[0 4], 'ZLim',[0 .5]);
  set(gca,'box','off', 'FontSize',14);
