
% ==============   V1 E-I NET CIRCUIT MODEL - J NEUROSCIENCE PLOTS  ==============

% MatLab code to run simulations and generate plots for:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.
%
% The whitened image training data used here is from:
%    Olshausen BA, Field DJ (1997). Sparse Coding with an Overcomplete
%        Basis Set: A Strategy Employed by V1?  Vision Research, 37: 3311-3325. 
%        http://redwood.berkeley.edu/bruno/sparsenet/
%
% Layers of model (high to low):
%    5) demo scripts: demo1, demo2, demo3, ...
%    4) RunVisionNetworkSimulation, NetModel_Figure
%    3) NetModel_InitV1eV1i
%    2) NetModel_Init, NetModel_UpdateFast, NetModel_Stats
%    1) CellGroup_Init, CellGroup_AddInput
%
% The learning rates that appear in the paper result from the combination
% of several different learning rates that appear in the simulator. The weight-change
% learning rate that appears in the paper is calculated as "lrateW" here:
%    lrateW = (model.lrateScale / model.simTimeStep) * model.learningRate
%               * model.cg_xxx.learningRate * model.cg_xxx.ib_yyy.learningRate
%           = .1 * model.learningRate * model.cg_xxx.ib_yyy.learningRate
%
%    lrate_in_E = .1 * .4 * .2  = .008
%    lrate_E_I  = .1 * .4 * .7  = .028
%    lrate_E_I  = .1 * .4 * 1.5 = .06


% run demo networks which draw plots in real-time (these are < 1 minute)
demo1_EINet
demo2_SAILnet

% train EI Net using fast learning (3 minutes)
demo3_TrainEINetFast

% train EI Net using a medium-speed learning (20 minutes)
demo3_TrainEINet

% run all training simulations presented in the paper
% (these take around 10 hours total on a 4-core Intel machine)
train1_SAILNetSlow_Fig3A        % 190 minutes
train2_EINetSlow_Fig3B          % 110 minutes
train3_EICount_Fig4             % 150 minutes
train4_IFastSpiking_Fig7        % 190 minutes

% plot all figures that appear in the paper using the training data
% saved during train1 to train4, above.
plot_figures

