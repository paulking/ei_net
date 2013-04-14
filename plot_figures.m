% draw_figures - Draw all figures that appear in paper after training
%
% Before this script is run, all "train" scripts must be run to prepare the data.
% Result data from the train scripts are saved in MAT files, which are then
% loaded here for plotting.
%
% Reference:
%    King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.


% ========================  FIGURE 2: SPIKE RASTER PLOTS  ============================

load experiment2_EINet

% Figure 2: Example spike raster plot (points on white background)
model.cellGroup{2}.displayColor        = [0 .7 0];            % use darker green for E cells
model.cellGroup{3}.displayColor        = [.9 0 0];
model.stats.keepSpikeHistory           = true;
model = RunVisionNetworkSimulation({'numInputSamples',100}, model);
NetModel_Figure(model, 'spikeRasterAll', 'figureNum',21, 'bgColor','white', 'markerSize',18, 'networkID',1);
set(gca, 'FontSize',16, 'XTick',10:10:50, 'YTick',0:100:400); % fix up plot labels



% ==============  FIGURE 3: RECEPTIVE FIELDS FOR SAILNET, EI-NET  =====================

% Figure 3A: Receptive fields for SAILnet
load experiment1_SAILnet
NetModel_Figure(model, 'STA', 'measureName','STA_V1e');


% Figure 3B: Receptive fields for EINet
load experiment2_EINet
NetModel_Figure(model, 'STA', 'measureName','STA_V1e');
NetModel_Figure(model, 'STA', 'measureName','STA_V1i');



% =================  FIGURE 4: SURFACE PLOT FOR E vs. I CELL COUNT  =================

% Figure 4 (J Neuroscience): surface plot of ecorr & resErr for varying # E and I cells

% load results and generate a 3D plot of residual error and E cell correlation
load experiment3_EICount
[X,Y] = meshgrid(results.propValues{2}, results.propValues{1});
figure(41);  set(41, 'Name', 'EI: E-cell Correlation');
  surf(X, Y, results.finalCorr);   ylim([0 600]);  zlim([.11 .17]);
  set(gca, 'FontSize',16, 'XTick',[0 20 40 60 80 100], 'YTick',[0 200 400 600]);
  set(gca, 'YDir','reverse');
  xlabel('# I cells', 'FontSize',12);  ylabel('# E cells', 'FontSize',12);
figure(42);  set(42, 'Name', 'EI: RMS Residual Error');
  surf(X, Y, results.finalResErr); ylim([0 600]);  zlim([.7 1]);
  set(gca, 'FontSize',16, 'XTick',[0 20 40 60 80 100], 'YTick',[0 200 400 600], 'ZTick',.7:.1:1);
  set(gca, 'YDir','reverse');
  xlabel('# I cells', 'FontSize',12);  ylabel('# E cells', 'FontSize',12);

% plot I correlation
figure(43);  set(43, 'Name','EI: I-cell Correlation'); 
  surf(X, Y, results.finalCorrI);  ylim([0 600]);
  set(gca, 'FontSize',16, 'XTick',[0 20 40 60 80 100], 'YTick',[0 200 400 600]);
  set(gca, 'YDir','reverse');
  xlabel('# I cells', 'FontSize',12);  ylabel('# E cells', 'FontSize',12);



% ===============  FIGURE 5: E-I WEIGHTS SORTED BY ORIENTATION/PHASE  ====================

% Figure 5 (J Neuroscience): Weight matrix sorted by orientation and phase
load experiment2_EINet

% analyze response of E and I cells to sine gratings
resultsE = AnalyzeVisualResponse(model, struct('cellGroupName', 'cg_V1e'));
resultsI = AnalyzeVisualResponse(model, struct('cellGroupName', 'cg_V1i'));

% get W = I->E weights (row = I, col = E)
[i_cg, i_ib] = NetModel_FindElement(model, 'cg_V1e.in_V1i');
W = model.cellGroup{i_cg}.inputBlock(i_ib).weight';
W_sorted = W;

% sort by orientation and plot
sortValE = [resultsE.cellStats.gaborOrientationMean];
sortValI = [resultsI.cellStats.gaborOrientationMean];
[~,idxE] = sort(sortValE);
[~,idxI] = sort(sortValI);
W_sorted = W(idxI, idxE);
figure(51);  set(51, 'Name', 'I->E weights sorted by Gabor ORIENTATION');
  imagesc(W_sorted);  colorbar;  colormap hot;  
  set(gca,'XTick',[100 200 300 400], 'YTick', [10 20 30 40], 'FontSize',14);
  colorbar('YTick', [0 5], 'FontSize', 16);
  xlabel('E cell ID', 'FontSize',12);  ylabel('I cell ID', 'FontSize',12);

% sort by phase and plot
sortValE = [resultsE.cellStats.gaborPhaseMean];
sortValI = [resultsI.cellStats.gaborPhaseMean];
[~,idxE] = sort(sortValE);
[~,idxI] = sort(sortValI);
W_sorted = W(idxI, idxE);
figure(52);  set(52, 'Name', 'I->E weights sorted by Gabor PHASE');
  imagesc(W_sorted);  colorbar;  colormap hot;  
  set(gca,'XTick',[100 200 300 400], 'YTick', [10 20 30 40], 'FontSize',14);
  colorbar('YTick', [0 5], 'FontSize', 16);
  xlabel('E cell ID', 'FontSize',12);  ylabel('I cell ID', 'FontSize',12);



% ===========  FIGURE 6: PSTH AND ACTIVITY PLOTS TO SHOW DYNAMICS  ================

% Figure 6 (J Neuroscience): network dynamics
%    panel A: fig 61 (below)
%    panel B: fig 62 (below)
%    panel C: fig 63 (below)
%    panel D: fig 64 (below)


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


% figure 6A: E spike rate and inhibitory feedback as function of net image input current
NetModel_Figure(model, 'cellGroupVars', 'figureNum',61, 'sourceName','cg_V1e', ...
        'varNames',{'ib_input', 'spikeRate','ib_V1i'}, 'showPDF',2,  ...
        'legendText',{'image input current to E cell (input\rightarrowE)', 'E cell spike rate', ...
                'inhibitory feedback to E cell'}, ...
        'fontSize',16, 'lineWidth',3, 'colormap',[0 .75 .5; 1 0 0], ...
        'title','figure 5: E cell activity as function of net image input');

% figure 6B: E-triggered PSTH showing total I feedback to E cells
% (timeBinVal is shifted 1/2 time step because inhibition isn't received until next time-step)
figure(62);  set(62, 'Name','E spike-triggered I input histogram');
  LW = 3;   LW_mean = 1.5;   LW_axes = 1;    % line widths
  plot(results.timeBinVal+.5, results.histInhibMean, '-r', 'LineWidth',LW);
  set(gca,'YLim',[0 2.7], 'YTick',0:2, 'XLim',[-15 15], 'XTick',-10:5:10);
  set(gca,'box','off', 'FontSize',16);
  hold on; plot(xlim(), [1 1]*mean(results.histInhibMean), '--k', 'LineWidth',LW_mean);   hold off;   % plot Y=avg inhib
  hold on; plot([0 0], ylim(), '-k');  hold off;                                % plot Y axis
  xlabel('\Delta time (simulation time steps)','FontSize',13)
  ylabel('Total feedback inhibition to E cells','FontSize',13)

% figure 6C: E spike rate as function of inhibitory feedback
% spike rate is in spikes-per-time unit
NetModel_Figure(model, 'cellGroupVars', 'figureNum',63, 'sourceName','cg_V1e', ...
        'varNames',{'ib_V1i', 'spikeRate'}, 'showPDF',false, ...
        'legendText',{'total inhibitory feedback to E cell (I\rightarrowE)', 'E cell spike rate'}, ...
        'fontSize',16, 'lineWidth',3, 'colormap',[0 .75 .5; 1 0 0], ...
        'title','figure 6: E spikeRate as function of I feedback');
set(gca,'XLim', [0 1.1], 'XTick', [0 .5 1], 'YLim',[0 .14], 'YTick',[0 .1]);

% figure 6D: I-triggered PSTH showing E spike rate of strong (green) vs. weak (blue) connections
figure(64);  set(64, 'Name','I triggered E spike rate, weak vs. strong');
  LW = 3;   LW_mean = 1.5;   LW_axes = 1;    % line widths
  i_cg_E = NetModel_FindElement(model, 'cg_V1e');
  targetSpikeRateE = model.cellGroup{i_cg_E}.targetSpikeRate;
  plot(results.timeBinVal, results.histMeanStrong, '-', 'LineWidth',LW,  'Color',[0 .75 .5]);
  set(gca,'box','off', 'FontSize',16);
  set(gca,'YLim',[0 .4], 'YTick',0:.1:.4, 'XLim',[-15 15], 'XTick',-10:5:10);
  hold on; plot(results.timeBinVal, results.histMeanWeak,   '-b', 'LineWidth',LW);   hold off;
  legend('strong targets', 'weak targets', 'Location','Northeast');
  hold on; plot(xlim(), [1 1]*targetSpikeRateE, '--k', 'LineWidth',LW_mean); hold off;  % plot Y=targetSR_E
  hold on; plot([0 0], ylim(), '-k');  hold off;                                % plot Y axis
  xlabel('\Delta time (simulation time steps)','FontSize',13);
  ylabel('E cell spike rate','FontSize',13);



% =====  FIGURE 7: SURFACE PLOT FOR I CELL FAST-SPIKING vs. CONNECTION WEIGHT MULTIPLIER  ======

% Figure 7 (J Neuroscience): surface plot of recorr and resErr for FS I cells

% load results and generate a 3D plot of residual error and E correlation
load experiment4_FSBW
[X,Y] = meshgrid(results.blockWeight, results.spikeRate/results.modelInitParams.cg_V1e.targetSpikeRate);
figure(71);  set(71, 'Name','SR vs BW: E-cell Correlation');
  surf(X, Y, results.finalCorr);   
  set(gca, 'FontSize',16, 'XTick',[1 2 4 6], 'YTick',[.5 1 2 3 4]);  ylim([.5 4]);  xlim([1 6]);
  xlabel('E\rightarrowI blockWeight', 'FontSize',12);  ylabel('I spikeRate (norm)', 'FontSize',12);
figure(72); set(72, 'Name','SR vs BW: RMS Residual Error');
  surf(X, Y, results.finalResErr); 
  set(gca, 'FontSize',16, 'XTick',[1 2 4 6], 'YTick',[.5 1 2 3 4]);  ylim([.5 4]);  xlim([1 6]);
  xlabel('E\rightarrowI blockWeight', 'FontSize',12);  ylabel('I spikeRate (norm)', 'FontSize',12);

