function [AnalysisResults] = Fig2_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data.Kernel
data.Kernel.GammaHbT.gamma = [];
data.Kernel.GammaHbT.HbT = [];
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(neuralBands)
        neuralBand = neuralBands{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            % pull HRFs from AnalysisResults.mat structure - dc shift each IR function by offset
            data.Kernel.IR.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.HRFs.(animalID).(neuralBand).cortLH.(behavior).IR_function_short - AnalysisResults.HRFs.(animalID).(neuralBand).cortLH.(behavior).IR_function_short(1);
            data.Kernel.IR.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.HRFs.(animalID).(neuralBand).cortRH.(behavior).IR_function_short - AnalysisResults.HRFs.(animalID).(neuralBand).cortRH.(behavior).IR_function_short(1);
            data.Kernel.IR_gamma.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.HRFs.(animalID).(neuralBand).cortLH.(behavior).IR_gammaFunction;
            data.Kernel.IR_gamma.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.HRFs.(animalID).(neuralBand).cortRH.(behavior).IR_gammaFunction;
        end
    end
end
% concatenate the data.Kernel from left and right into a single data.Kernel set
for dd = 1:length(neuralBands)
    neuralBand = neuralBands{1,dd};
    for ee = 1:length(behaviors)
        behavior = behaviors{1,ee};
        % impulse
        data.Kernel.IR.(neuralBand).(behavior).cat = cat(1,data.Kernel.IR.(neuralBand).(behavior).cortLH,data.Kernel.IR.(neuralBand).(behavior).cortRH);
        data.Kernel.IR_gamma.(neuralBand).(behavior).cat = cat(1,data.Kernel.IR_gamma.(neuralBand).(behavior).cortLH,data.Kernel.IR_gamma.(neuralBand).(behavior).cortRH);
    end
end
% mean and std of each arousal-state
for ii = 1:length(neuralBands)
    neuralBand = neuralBands{1,ii};
    for jj = 1:length(behaviors)
        behavior = behaviors{1,jj};
        % IR mean
        data.Kernel.IR.(neuralBand).(behavior).mean = mean(data.Kernel.IR.(neuralBand).(behavior).cat,1);
    end
end
%% gamma HRF based on impulse deconvolution
for kk = 1:length(neuralBands)
    neuralBand = neuralBands{1,kk};
    for ll = 1:length(behaviors)
        behavior = behaviors{1,ll};
        [peak,peakIndex] = max(data.Kernel.IR.(neuralBand).(behavior).mean);
        peakTime = peakIndex/samplingRate;
        threeQuarterMax = max(data.Kernel.IR.(neuralBand).(behavior).mean)/(4/3);
        index1 = find(data.Kernel.IR.(neuralBand).(behavior).mean >= threeQuarterMax,1,'first');
        % find where the data.Kernel last rises above half the max.
        index2 = find(data.Kernel.IR.(neuralBand).(behavior).mean >= threeQuarterMax,1,'last');
        threeQuarterWidth = (index2 - index1 + 1)/samplingRate; % FWHM in indexes.
        initVals = [peak,peakTime,threeQuarterWidth];
        % create gamma function based on impulse values
        t = 0:1/samplingRate:5;
        IR_a = ((initVals(2)/initVals(3))^2*8*log10(2));
        IR_beta = ((initVals(3)^2)/initVals(2)/8/log10(2));
        data.Kernel.IR_gamma.(neuralBand).(behavior).repFunc = initVals(1)*(t/initVals(2)).^IR_a.*exp((t - initVals(2))/(-1*IR_beta));
    end
end
%% Fig. 2
summaryFigure = figure('Name','Fig2 (A-H)');
sgtitle('Figure 2 - Turner et al. 2020')
%% [2A - top] cortical MUA contra stim
ax1 = subplot(3,4,1);
% p1 = plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortMUA,'color',colors_HRF2020('rich black'),'LineWidth',2);
hold on
% plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortMUA + data.Stim.Contra.std_CortMUA,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
% plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortMUA - data.Stim.Contra.std_CortMUA,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
p2 = plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortGam,'color',colors_HRF2020('dark cyan'),'LineWidth',2);
hold on
p1 = plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortMUA,'color',colors_HRF2020('rich black'),'LineWidth',2);
% plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortGam + data.Stim.Contra.std_CortGam,'color',colors_HRF2020('cyan'),'LineWidth',0.5)
% plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortGam - data.Stim.Contra.std_CortGam,'color',colors_HRF2020('cyan'),'LineWidth',0.5)
title('[2A top] Stim cortical neural')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2],'MUA','Gamma')
axis square
set(gca,'box','off')
xlim([-2,5])
ylim([-100,1200])
ax1.TickLength = [0.03,0.03];
%% [2B - top] HbT contra stim
ax5 = subplot(3,4,5);
plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_HbT,'color',colors_HRF2020('rich black'),'LineWidth',2)
hold on
plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_HbT + data.Stim.Contra.std_HbT,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_HbT - data.Stim.Contra.std_HbT,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
title('[2A bottom] Stim \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
xlim([-2,5])
ylim([-10,25])
ax5.TickLength = [0.03,0.03];
%% [2A - bottom] moderate whisks cortical MUA
ax2 = subplot(3,4,2);
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortMUA,'color',colors_HRF2020('rich black'),'LineWidth',2);
hold on
% plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortMUA + data.Whisk.IntermediateWhisks.stdCortMUA,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
% plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortMUA - data.Whisk.IntermediateWhisks.stdCortMUA,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortGam,'color',colors_HRF2020('dark cyan'),'LineWidth',2);
% plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortGam + data.Whisk.IntermediateWhisks.stdCortGam,'color',colors_HRF2020('cyan'),'LineWidth',0.5)
% plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortGam - data.Whisk.IntermediateWhisks.stdCortGam,'color',colors_HRF2020('cyan'),'LineWidth',0.5)
title('[2B top] Whisk cortical neural')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
xlim([-2,5])
ylim([-10,120])
ax2.TickLength = [0.03,0.03];
%% [2B - bottom] moderate whisks HbT
ax6 = subplot(3,4,6);
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanHbT,'color',colors_HRF2020('rich black'),'LineWidth',2);
hold on
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanHbT + data.Whisk.IntermediateWhisks.stdHbT,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanHbT - data.Whisk.IntermediateWhisks.stdHbT,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
title('[2B bottom] Whisk \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
xlim([-2,5])
ylim([-10,25])
ax6.TickLength = [0.03,0.03];
%% [2C] Gamma-HbT XCorr
freq = 30;
lag = 5;
ax1 = subplot(3,4,[3,7]);
plot(data.XCorr.Rest.meanGAM_lags,data.XCorr.Rest.meanHbTvGAMxcVals,'color',colorRest,'LineWidth',2)
hold on
plot(data.XCorr.Rest.meanGAM_lags,data.XCorr.Rest.meanHbTvGAMxcVals + data.XCorr.Rest.stdHbTvGAMxcVals,'color',colorRest,'LineWidth',0.5)
plot(data.XCorr.Rest.meanGAM_lags,data.XCorr.Rest.meanHbTvGAMxcVals - data.XCorr.Rest.stdHbTvGAMxcVals,'color',colorRest,'LineWidth',0.5)
plot(data.XCorr.NREM.meanGAM_lags,data.XCorr.NREM.meanHbTvGAMxcVals,'color',colorNREM,'LineWidth',2)
plot(data.XCorr.NREM.meanGAM_lags,data.XCorr.NREM.meanHbTvGAMxcVals + data.XCorr.NREM.stdHbTvGAMxcVals,'color',colorNREM,'LineWidth',0.5)
plot(data.XCorr.NREM.meanGAM_lags,data.XCorr.NREM.meanHbTvGAMxcVals - data.XCorr.NREM.stdHbTvGAMxcVals,'color',colorNREM,'LineWidth',0.5)
plot(data.XCorr.REM.meanGAM_lags,data.XCorr.REM.meanHbTvGAMxcVals,'color',colorREM,'LineWidth',2)
plot(data.XCorr.REM.meanGAM_lags,data.XCorr.REM.meanHbTvGAMxcVals + data.XCorr.REM.stdHbTvGAMxcVals,'color',colorREM,'LineWidth',0.5)
plot(data.XCorr.REM.meanGAM_lags,data.XCorr.REM.meanHbTvGAMxcVals - data.XCorr.REM.stdHbTvGAMxcVals,'color',colorREM,'LineWidth',0.5)
title('[2C] Gam-[HbT] XCorr')
xticks([-lag*freq,-lag*freq/2,0,lag*freq/2,lag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-lag*freq,lag*freq])
xlabel('Lags (s)')
ylabel({'Corr. coefficient';'Gamma vs. \Delta[HbT] (\muM)'})
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ylim([-0.1,0.5])
%% [2D] MUA-HbT XCorr
ax1 = subplot(3,4,[4,8]);
plot(data.XCorr.Rest.meanMUA_lags,data.XCorr.Rest.meanHbTvMUAxcVals,'color',colorRest,'LineWidth',2)
hold on
plot(data.XCorr.Rest.meanMUA_lags,data.XCorr.Rest.meanHbTvMUAxcVals + data.XCorr.Rest.stdHbTvMUAxcVals,'color',colorRest,'LineWidth',0.5)
plot(data.XCorr.Rest.meanMUA_lags,data.XCorr.Rest.meanHbTvMUAxcVals - data.XCorr.Rest.stdHbTvMUAxcVals,'color',colorRest,'LineWidth',0.5)
plot(data.XCorr.NREM.meanMUA_lags,data.XCorr.NREM.meanHbTvMUAxcVals,'color',colorNREM,'LineWidth',2)
plot(data.XCorr.NREM.meanMUA_lags,data.XCorr.NREM.meanHbTvMUAxcVals + data.XCorr.NREM.stdHbTvMUAxcVals,'color',colorNREM,'LineWidth',0.5)
plot(data.XCorr.NREM.meanMUA_lags,data.XCorr.NREM.meanHbTvMUAxcVals - data.XCorr.NREM.stdHbTvMUAxcVals,'color',colorNREM,'LineWidth',0.5)
plot(data.XCorr.REM.meanMUA_lags,data.XCorr.REM.meanHbTvMUAxcVals,'color',colorREM,'LineWidth',2)
plot(data.XCorr.REM.meanMUA_lags,data.XCorr.REM.meanHbTvMUAxcVals + data.XCorr.REM.stdHbTvMUAxcVals,'color',colorREM,'LineWidth',0.5)
plot(data.XCorr.REM.meanMUA_lags,data.XCorr.REM.meanHbTvMUAxcVals - data.XCorr.REM.stdHbTvMUAxcVals,'color',colorREM,'LineWidth',0.5)
title('[2D] MUA-[HbT] XCorr')
xticks([-lag*freq,-lag*freq/2,0,lag*freq/2,lag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-lag*freq,lag*freq])
xlabel('Lags (s)')
ylabel({'Corr. coefficient';'MUA vs. \Delta[HbT] (\muM)'})
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ylim([-0.1,0.5])
%% [2E] gamma-band impulse function
subplot(3,4,9);
p1 = plot(t,data.Kernel.IR.gammaBandPower.Rest.mean,'color',colorRest,'LineWidth',2);
hold on
p2 = plot(t,data.Kernel.IR.gammaBandPower.Whisk.mean,'color',colorWhisk,'LineWidth',2);
p3 = plot(t,data.Kernel.IR.gammaBandPower.Contra.mean,'color',colorStim,'LineWidth',2);
p4 = plot(t,data.Kernel.IR.gammaBandPower.NREM.mean,'color',colorNREM,'LineWidth',2);
p5 = plot(t,data.Kernel.IR.gammaBandPower.REM.mean,'color',colorREM,'LineWidth',2);
p6 = plot(t,data.Kernel.IR.gammaBandPower.All.mean,'color',colorAll,'LineWidth',2);
title('[2E] Gamma [30-100 Hz] IRF')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
legend([p1,p2,p3,p4,p5,p6],'Rest','Whisk','Stim','NREM','REM','All','Location','NorthEast')
axis square
xlim([0,5])
ylim([-0.1,0.8])
set(gca,'box','off')
%% [2F] gamma-band gamma function
subplot(3,4,10);
plot(t,data.Kernel.IR_gamma.gammaBandPower.Rest.repFunc,'color',colorRest,'LineWidth',2);
hold on
plot(t,data.Kernel.IR_gamma.gammaBandPower.Whisk.repFunc,'color',colorWhisk,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.gammaBandPower.Contra.repFunc,'color',colorStim,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.gammaBandPower.NREM.repFunc,'color',colorNREM,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.gammaBandPower.REM.repFunc,'color',colorREM,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.gammaBandPower.All.repFunc,'color',colorAll,'LineWidth',2);
title('[2F] Gamma-band [30-100 Hz] GF')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
axis square
xlim([0,5])
ylim([-0.1,0.8])
set(gca,'box','off')
%% [2G] MUA impulse function
subplot(3,4,11);
plot(t,data.Kernel.IR.muaPower.Rest.mean,'color',colorRest,'LineWidth',2);
hold on
plot(t,data.Kernel.IR.muaPower.Whisk.mean,'color',colorWhisk,'LineWidth',2);
plot(t,data.Kernel.IR.muaPower.Contra.mean,'color',colorStim,'LineWidth',2);
plot(t,data.Kernel.IR.muaPower.NREM.mean,'color',colorNREM,'LineWidth',2);
plot(t,data.Kernel.IR.muaPower.REM.mean,'color',colorREM,'LineWidth',2);
plot(t,data.Kernel.IR.muaPower.All.mean,'color',colorAll,'LineWidth',2);
title('[2G] MUA [300-3000 Hz] IRF')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
axis square
xlim([0,5])
ylim([-0.2,1.5])
set(gca,'box','off')
%% [2H] MUA gamma function
subplot(3,4,12);
plot(t,data.Kernel.IR_gamma.muaPower.Rest.repFunc,'color',colorRest,'LineWidth',2);
hold on
plot(t,data.Kernel.IR_gamma.muaPower.Whisk.repFunc,'color',colorWhisk,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.muaPower.Contra.repFunc,'color',colorStim,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.muaPower.NREM.repFunc,'color',colorNREM,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.muaPower.REM.repFunc,'color',colorREM,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.muaPower.All.repFunc,'color',colorAll,'LineWidth',2);
title('[2H] MUA [300-3000 Hz] GF')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
axis square
xlim([0,5])
ylim([-0.2,1.5])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig2']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig2'])
    %     %% text diary
    %     diaryFile = [dirpath 'Fig1-S2_Statistics.txt'];
    %     if exist(diaryFile,'file') == 2
    %         delete(diaryFile)
    %     end
    %     diary(diaryFile)
    %     diary on
    %     % text values
    %     disp('======================================================================================================================')
    %     disp('[1-S2] Text values for gamma/HbT/reflectance changes')
    %     disp('======================================================================================================================')
    %     disp('----------------------------------------------------------------------------------------------------------------------')
    %     % cortical MUA/LFP
    %     [~,index] = max(data.Stim.Contra.mean_CortMUA);
    %     disp(['Contra stim Cort gamma MUA P/P (%): ' num2str(round(data.Stim.Contra.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Stim.Contra.std_CortMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Ipsi.mean_CortMUA);
    %     disp(['Ipsil stim Cort gamma MUA P/P (%): ' num2str(round(data.Stim.Ipsi.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_CortMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Auditory.mean_CortMUA);
    %     disp(['Audit stim Cort gamma MUA P/P (%): ' num2str(round(data.Stim.Auditory.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Stim.Auditory.std_CortMUA(index),1))]); disp(' ')
    %     % cortical LFP
    %     disp(['Contra stim Cort gamma LFP P/P (%): ' num2str(round(data.Stim.Contra.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Stim.Contra.std_CortS_Gam,1))]); disp(' ')
    %     disp(['Ipsil stim Cort gamma LFP P/P (%): ' num2str(round(data.Stim.Ipsi.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_CortS_Gam,1))]); disp(' ')
    %     disp(['Audit stim Cort gamma LFP P/P (%): ' num2str(round(data.Stim.Auditory.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Stim.Auditory.std_CortS_Gam,1))]); disp(' ')
    %     % hippocampal MUA
    %     [~,index] = max(data.Stim.Contra.mean_HipMUA);
    %     disp(['Contra stim Hip gamma MUA P/P (%): ' num2str(round(data.Stim.Contra.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Stim.Contra.std_HipMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Ipsi.mean_HipMUA);
    %     disp(['Ipsil stim Hip gamma MUA P/P (%): ' num2str(round(data.Stim.Ipsi.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_HipMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Auditory.mean_HipMUA);
    %     disp(['Audit stim Hip gamma MUA P/P (%): ' num2str(round(data.Stim.Auditory.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Stim.Auditory.std_HipMUA(index),1))]); disp(' ')
    %     % hippocampal LFP
    %     disp(['Contra stim Hip gamma LFP P/P (%): ' num2str(round(data.Stim.Contra.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Stim.Contra.std_HipS_Gam,1))]); disp(' ')
    %     disp(['Ipsil stim Hip gamma LFP P/P (%): ' num2str(round(data.Stim.Ipsi.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_HipS_Gam,1))]); disp(' ')
    %     disp(['Auditory stim Hip gamma LFP P/P (%): ' num2str(round(data.Stim.Auditory.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Stim.Auditory.std_HipS_Gam,1))]); disp(' ')
    %     % HbT
    %     [~,index] = max(data.Stim.Contra.mean_HbT);
    %     disp(['Contra stim [HbT] (uM): ' num2str(round(data.Stim.Contra.mean_HbT(index),1)) ' +/- ' num2str(round(data.Stim.Contra.std_HbT(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Ipsi.mean_HbT);
    %     disp(['Ipsil stim [HbT] (uM): ' num2str(round(data.Stim.Ipsi.mean_HbT(index),1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_HbT(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Auditory.mean_HbT);
    %     disp(['Audit stim [HbT] (uM): ' num2str(round(data.Stim.Auditory.mean_HbT(index),1)) ' +/- ' num2str(round(data.Stim.Auditory.std_HbT(index),1))]); disp(' ')
    %     % R/R
    %     [~,index] = min(data.Stim.Contra.mean_CBV);
    %     disp(['Contra stim refl R/R (%): ' num2str(round(data.Stim.Contra.mean_CBV(index),1)) ' +/- ' num2str(round(data.Stim.Contra.std_CBV(index),1))]); disp(' ')
    %     [~,index] = min(data.Stim.Ipsi.mean_CBV);
    %     disp(['Ipsil stim refl R/R (%): ' num2str(round(data.Stim.Ipsi.mean_CBV(index),1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_CBV(index),1))]); disp(' ')
    %     [~,index] = min(data.Stim.Auditory.mean_CBV);
    %     disp(['Audit stim refl R/R (%): ' num2str(round(data.Stim.Auditory.mean_CBV(index),1)) ' +/- ' num2str(round(data.Stim.Auditory.std_CBV(index),1))]); disp(' ')
    %     disp('----------------------------------------------------------------------------------------------------------------------')
    %     diary off
    %
    %
    %
    %
    %
    %     %%
    %      % text values
    %     disp('======================================================================================================================')
    %     disp('[1-S3] Text values for gamma/HbT changes')
    %     disp('======================================================================================================================')
    %     disp('----------------------------------------------------------------------------------------------------------------------')
    %      % cortical MUA/LFP
    %     [~,index] = max(data.Whisk.ShortWhisks.meanCortMUA);
    %     disp(['Brief whisk Cort gamma MUA P/P (%): ' num2str(round(data.Whisk.ShortWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.stdCortMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.IntermediateWhisks.meanCortMUA);
    %     disp(['Moderate whisk Cort gamma MUA P/P (%): ' num2str(round(data.Whisk.IntermediateWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.stdCortMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.LongWhisks.meanCortMUA);
    %     disp(['Extended whisk Cort gamma MUA P/P (%): ' num2str(round(data.Whisk.LongWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.stdCortMUA(index),1))]); disp(' ')
    %     % cortical LFP
    %     disp(['Brief whisk Cort gamma LFP P/P (%): ' num2str(round(data.Whisk.ShortWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.std_CortS_Gam,1))]); disp(' ')
    %     disp(['Moderate whisk Cort gamma LFP P/P (%): ' num2str(round(data.Whisk.IntermediateWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.std_CortS_Gam,1))]); disp(' ')
    %     disp(['Extended whisk Cort gamma LFP P/P (%): ' num2str(round(data.Whisk.LongWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.std_CortS_Gam,1))]); disp(' ')
    %     % hippocampal MUA
    %     [~,index] = max(data.Whisk.ShortWhisks.meanHipMUA);
    %     disp(['Brief whisk Hip gamma MUA P/P (%): ' num2str(round(data.Whisk.ShortWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.stdHipMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.IntermediateWhisks.meanHipMUA);
    %     disp(['Moderate whisk Hip gamma MUA P/P (%): ' num2str(round(data.Whisk.IntermediateWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.stdHipMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.LongWhisks.meanHipMUA);
    %     disp(['Extended whisk Hip gamma MUA P/P (%): ' num2str(round(data.Whisk.LongWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.stdHipMUA(index),1))]); disp(' ')
    %     % hippocampal LFP
    %     disp(['Brief whisk Hip gamma LFP P/P (%): ' num2str(round(data.Whisk.ShortWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.std_HipS_Gam,1))]); disp(' ')
    %     disp(['Moderate whisk Hip gamma LFP P/P (%): ' num2str(round(data.Whisk.IntermediateWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.std_HipS_Gam,1))]); disp(' ')
    %     disp(['Extended whisk Hip gamma LFP P/P (%): ' num2str(round(data.Whisk.LongWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.std_HipS_Gam,1))]); disp(' ')
    %     % HbT
    %     [~,index] = max(data.Whisk.ShortWhisks.meanHbT);
    %     disp(['Brief whisk [HbT] (uM): ' num2str(round(data.Whisk.ShortWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.stdHbT(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.IntermediateWhisks.meanHbT);
    %     disp(['Moderate whisk [HbT] (uM): ' num2str(round(data.Whisk.IntermediateWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.stdHbT(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.LongWhisks.meanHbT);
    %     disp(['Extended whisk [HbT] (uM): ' num2str(round(data.Whisk.LongWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.stdHbT(index),1))]); disp(' ')
    %     % R/R
    %     [~,index] = min(data.Whisk.ShortWhisks.meanCBV);
    %     disp(['Brief whisk refl R/R (%): ' num2str(round(data.Whisk.ShortWhisks.meanCBV(index),1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.stdCBV(index),1))]); disp(' ')
    %     [~,index] = min(data.Whisk.IntermediateWhisks.meanCBV);
    %     disp(['Moderate whisk refl R/R (%): ' num2str(round(data.Whisk.IntermediateWhisks.meanCBV(index),1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.stdCBV(index),1))]); disp(' ')
    %     [~,index] = min(data.Whisk.LongWhisks.meanCBV);
    %     disp(['Extended whisk refl R/R (%): ' num2str(round(data.Whisk.LongWhisks.meanCBV(index),1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.stdCBV(index),1))]); disp(' ')
    %     disp('----------------------------------------------------------------------------------------------------------------------')
    %     diary off
end

end
