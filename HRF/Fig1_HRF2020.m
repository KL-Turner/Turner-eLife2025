function [AnalysisResults] = Fig1_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel 1 for Turner_Drew_HRF2020
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
colorRfcREM = [(190/256),(30/256),(45/256)];
% colorRest = [(0/256),(166/256),(81/256)];
% colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
% colorNREM = [(191/256),(0/256),(255/256)];
% colorREM = [(254/256),(139/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% information and data for example
if isfield(AnalysisResults,'ExampleTrials') == false
    AnalysisResults.ExampleTrials = [];
end
if isfield(AnalysisResults.ExampleTrials,'T122AA') == true
    dsFs = AnalysisResults.ExampleTrials.T122A.dsFs;
    filtEMG = AnalysisResults.ExampleTrials.T122A.filtEMG;
    filtForceSensor = AnalysisResults.ExampleTrials.T122A.filtForceSensor;
    filtWhiskerAngle = AnalysisResults.ExampleTrials.T122A.filtWhiskerAngle;
    heartRate = AnalysisResults.ExampleTrials.T122A.heartRate;
    filt_HbT = AnalysisResults.ExampleTrials.T122A.filt_HbT;
    T = AnalysisResults.ExampleTrials.T122A.T;
    F = AnalysisResults.ExampleTrials.T122A.F;
    cortical_LHnormS = AnalysisResults.ExampleTrials.T122A.cortical_LHnormS;
    hippocampusNormS = AnalysisResults.ExampleTrials.T122A.hippocampusNormS;
else
    animalID = 'T122';
    dataLocation = [rootFolder '\' animalID '\Single Hemisphere\'];
    cd(dataLocation)
    exampleProcDataFileID = 'T122_200208_15_35_11_ProcData.mat';
    load(exampleProcDataFileID,'-mat')
    exampleSpecDataFileID = 'T122_200208_15_35_11_SpecDataA.mat';
    load(exampleSpecDataFileID,'-mat')
    exampleBaselineFileID = 'T122_RestingBaselines.mat';
    load(exampleBaselineFileID,'-mat')
    [~,fileDate,~] = GetFileInfo_HRF2020(exampleProcDataFileID);
    strDay = ConvertDate_HRF2020(fileDate);
    dsFs = ProcData.notes.dsFs;
    % setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
    [z1,p1,k1] = butter(4,10/(dsFs/2),'low');
    [sos1,g1] = zp2sos(z1,p1,k1);
    [z2,p2,k2] = butter(4,0.5/(dsFs/2),'low');
    [sos2,g2] = zp2sos(z2,p2,k2);
    % whisker angle
    filtWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
    % force sensor
    filtForceSensor = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
    % EMG
    EMG = ProcData.data.EMG.emg;
    normEMG = EMG - RestingBaselines.manualSelection.EMG.emg.(strDay);
    filtEMG = filtfilt(sos1,g1,normEMG);
    % heart rate
    heartRate = ProcData.data.heartRate;
    % HbT data (LH)
    HbT = ProcData.data.CBV_HbT.adjBarrels;
    filt_HbT = filtfilt(sos2,g2,HbT);
    % cortical and hippocampal spectrograms
    cortical_LHnormS = SpecData.cortical_LH.normS.*100;
    hippocampusNormS = SpecData.hippocampus.normS.*100;
    T = SpecData.cortical_LH.T;
    F = SpecData.cortical_LH.F;
    % update analysis structure
    AnalysisResults.ExampleTrials.T122A.dsFs = dsFs;
    AnalysisResults.ExampleTrials.T122A.filtEMG = filtEMG;
    AnalysisResults.ExampleTrials.T122A.filtForceSensor = filtForceSensor;
    AnalysisResults.ExampleTrials.T122A.filtWhiskerAngle = filtWhiskerAngle;
    AnalysisResults.ExampleTrials.T122A.heartRate = heartRate;
    AnalysisResults.ExampleTrials.T122A.filt_HbT = filt_HbT;
    AnalysisResults.ExampleTrials.T122A.T = T;
    AnalysisResults.ExampleTrials.T122A.F = F;
    AnalysisResults.ExampleTrials.T122A.cortical_LHnormS = cortical_LHnormS;
    AnalysisResults.ExampleTrials.T122A.hippocampusNormS = hippocampusNormS;
    % save results
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end
%% Fig. 1
summaryFigure = figure('Name','Fig1 (e-j)');
sgtitle('Figure 1 - Turner et al. HRF 2020')
%% [1e-j] IOS sleep example
% EMG and force sensor
ax1 = subplot(6,1,1);
p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors_HRF2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','power (a.u.)'})
ylim([-2,2.5])
yyaxis right
p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','Pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([30,90,150,210,270,330,390,450,510,570,630])
xlim([30,630])
ylim([-0.1,1])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_HRF2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% whisker angle and heart rate
ax2 = subplot(6,1,2);
p3 = plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors_HRF2020('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
xlim([30,630])
ylim([-10,50])
yyaxis right
p4 = plot((1:length(heartRate)),heartRate,'color',colors_HRF2020('deep carrot orange'),'LineWidth',0.5);
ylabel({'Heart rate','Freq (Hz)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p3,p4],'Whisker angle','Heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([30,90,150,210,270,330,390,450,510,570,630])
xlim([30,630])
ylim([5,10])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_HRF2020('rich black');
ax2.YAxis(2).Color = colors_HRF2020('deep carrot orange');
% HbT and behavioral indeces
ax34 =subplot(6,1,[3,4]);
p5 = plot((1:length(filt_HbT))/dsFs,filt_HbT,'color',colors_HRF2020('sapphire'),'LineWidth',1);
hold on
x1 = xline(30,'color',colorRfcNREM,'LineWidth',2);
x2 = xline(73,'color',colorRfcREM,'LineWidth',2);
x3 = xline(200,'color',colorRfcAwake,'LineWidth',2);
xline(227,'color',colorRfcNREM,'LineWidth',2);
xline(299,'color',colorRfcAwake,'LineWidth',2);
xline(313,'color',colorRfcNREM,'LineWidth',2);
xline(390,'color',colorRfcAwake,'LineWidth',2);
xline(443,'color',colorRfcNREM,'LineWidth',2);
xline(597,'color',colorRfcAwake,'LineWidth',2);
ylabel('\Delta[HbT] (\muM)')
legend([p5,x3,x1,x2],'Barrels','Awake','NREM','REM')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([30,90,150,210,270,330,390,450,510,570,630])
xlim([30,630])
ylim([-35,135])
ax34.TickLength = [0.01,0.01];
% left cortical electrode spectrogram
ax5 = subplot(6,1,5);
SemilogImageSC_HRF2020(T,F,cortical_LHnormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
ylabel({'LH cort LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([30,90,150,210,270,330,390,450,510,570,630])
xlim([30,630])
ax5.TickLength = [0.01,0.01];
% hippocampal electrode spectrogram
ax6 = subplot(6,1,6);
SemilogImageSC_HRF2020(T,F,hippocampusNormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([30,90,150,210,270,330,390,450,510,570,630])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([30,630])
ax6.TickLength = [0.01,0.01];
% axes properties
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1']);
    % remove surface subplots because they take forever to render
    cla(ax5);
    set(ax5,'YLim',[1,99]);
    cla(ax6);
    set(ax6,'YLim',[1,99]);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1'])
    close(summaryFigure)
    %% subplot figures
    subplotImgs = figure;
    % example 1 LH cortical LFP
    subplot(2,1,1);
    SemilogImageSC_HRF2020(T,F,cortical_LHnormS,'y')
    caxis([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([30,630])
    % example 1 hippocampal LFP
    subplot(2,1,2);
    SemilogImageSC_HRF2020(T,F,hippocampusNormS,'y')
    caxis([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([30,630])
    print('-painters','-dtiffn',[dirpath 'Fig1_SpecImages'])
    close(subplotImgs)
    %% Fig. 1
    figure('Name','Fig1 (e-j)');
    sgtitle('Figure Panel 1 - Turner et al. 2020')
    %% [1e-j] IOS sleep example
    % EMG and force sensor
    ax1 = subplot(6,1,1);
    p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors_HRF2020('rich black'),'LineWidth',0.5);
    ylabel({'EMG','power (a.u.)'})
    ylim([-2,2.5])
    yyaxis right
    p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
    ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
    legend([p1,p2],'EMG','Pressure')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([30,90,150,210,270,330,390,450,510,570,630])
    xlim([30,630])
    ylim([-0.1,1])
    ax1.TickLength = [0.01,0.01];
    ax1.YAxis(1).Color = colors_HRF2020('rich black');
    ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
    % whisker angle and heart rate
    ax2 = subplot(6,1,2);
    p3 = plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors_HRF2020('rich black'),'LineWidth',0.5);
    ylabel({'Whisker','angle (deg)'})
    xlim([30,630])
    ylim([-10,50])
    yyaxis right
    p4 = plot((1:length(heartRate)),heartRate,'color',colors_HRF2020('deep carrot orange'),'LineWidth',0.5);
    ylabel({'Heart rate','Freq (Hz)'},'rotation',-90,'VerticalAlignment','bottom')
    legend([p3,p4],'Whisker angle','Heart rate')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([30,90,150,210,270,330,390,450,510,570,630])
    xlim([30,630])
    ylim([5,10])
    ax2.TickLength = [0.01,0.01];
    ax2.YAxis(1).Color = colors_HRF2020('rich black');
    ax2.YAxis(2).Color = colors_HRF2020('deep carrot orange');
    % HbT and behavioral indeces
    ax34 =subplot(6,1,[3,4]);
    p5 = plot((1:length(filt_HbT))/dsFs,filt_HbT,'color',colors_HRF2020('sapphire'),'LineWidth',1);
    hold on
    x1 = xline(30,'color',colorRfcNREM,'LineWidth',2);
    x2 = xline(73,'color',colorRfcREM,'LineWidth',2);
    x3 = xline(200,'color',colorRfcAwake,'LineWidth',2);
    xline(227,'color',colorRfcNREM,'LineWidth',2);
    xline(299,'color',colorRfcAwake,'LineWidth',2);
    xline(313,'color',colorRfcNREM,'LineWidth',2);
    xline(390,'color',colorRfcAwake,'LineWidth',2);
    xline(443,'color',colorRfcNREM,'LineWidth',2);
    xline(597,'color',colorRfcAwake,'LineWidth',2);
    ylabel('\Delta[HbT] (\muM)')
    legend([p5,x3,x1,x2],'Barrels','Awake','NREM','REM')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([30,90,150,210,270,330,390,450,510,570,630])
    xlim([30,630])
    ylim([-35,135])
    ax34.TickLength = [0.01,0.01];
    % left cortical electrode spectrogram
    ax5 = subplot(6,1,5);
    SemilogImageSC_HRF2020(T,F,cortical_LHnormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,100])
    ylabel({'LH cort LFP','Freq (Hz)'})
    set(gca,'Yticklabel','10^1')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([30,90,150,210,270,330,390,450,510,570,630])
    xlim([30,630])
    ax5.TickLength = [0.01,0.01];
    % hippocampal electrode spectrogram
    ax6 = subplot(6,1,6);
    SemilogImageSC_HRF2020(T,F,hippocampusNormS,'y')
    axis xy
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,100])
    xlabel('Time (min)')
    ylabel({'Hipp LFP','Freq (Hz)'})
    set(gca,'box','off')
    xticks([30,90,150,210,270,330,390,450,510,570,630])
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    xlim([30,630])
    ax6.TickLength = [0.01,0.01];
    % axes properties
    ax1Pos = get(ax1,'position');
    ax5Pos = get(ax5,'position');
    ax6Pos = get(ax6,'position');
    ax5Pos(3:4) = ax1Pos(3:4);
    ax6Pos(3:4) = ax1Pos(3:4);
    set(ax5,'position',ax5Pos);
    set(ax6,'position',ax6Pos);
end

end
