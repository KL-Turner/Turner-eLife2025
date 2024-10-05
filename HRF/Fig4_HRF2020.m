function [AnalysisResults] = Fig4_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
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
% kernels
gammaKernel = AnalysisResults.HRFs.T122.gammaBandPower.cortLH.All.IR_gammaFunction;
muaKernel = AnalysisResults.HRFs.T122.muaPower.cortLH.All.IR_gammaFunction;
if isfield(AnalysisResults.ExampleTrials,'T122B') == true
    % example B
    dsFs_B = AnalysisResults.ExampleTrials.T122B.dsFs;
    binForceSensor_B = AnalysisResults.ExampleTrials.T122B.binForceSensor;
    binWhiskerAngle_B = AnalysisResults.ExampleTrials.T122B.binWhiskerAngle;
    actHbT_B = AnalysisResults.ExampleTrials.T122B.actHbT;
    gammaPredHbT_B = AnalysisResults.ExampleTrials.T122B.gammaPredHbT;
    muaPredHbT_B = AnalysisResults.ExampleTrials.T122B.muaPredHbT;
    gammaR_B = AnalysisResults.ExampleTrials.T122B.gammaR;
    gammaR2_B = AnalysisResults.ExampleTrials.T122B.gammaR2;
    muaR_B = AnalysisResults.ExampleTrials.T122B.muaR;
    muaR2_B = AnalysisResults.ExampleTrials.T122B.muaR2;
    T_B = AnalysisResults.ExampleTrials.T122B.T;
    F_B = AnalysisResults.ExampleTrials.T122B.F;
    cortical_LHnormS_B = AnalysisResults.ExampleTrials.T122B.cortical_LHnormS;
    hippocampusNormS_B = AnalysisResults.ExampleTrials.T122B.hippocampusNormS;
    % example C
    dsFs_C = AnalysisResults.ExampleTrials.T122C.dsFs;
    binForceSensor_C = AnalysisResults.ExampleTrials.T122C.binForceSensor;
    binWhiskerAngle_C = AnalysisResults.ExampleTrials.T122C.binWhiskerAngle;
    actHbT_C = AnalysisResults.ExampleTrials.T122C.actHbT;
    gammaPredHbT_C = AnalysisResults.ExampleTrials.T122C.gammaPredHbT;
    muaPredHbT_C = AnalysisResults.ExampleTrials.T122C.muaPredHbT;
    gammaR_C = AnalysisResults.ExampleTrials.T122C.gammaR;
    gammaR2_C = AnalysisResults.ExampleTrials.T122C.gammaR2;
    muaR_C = AnalysisResults.ExampleTrials.T122C.muaR;
    muaR2_C = AnalysisResults.ExampleTrials.T122C.muaR2;
    T_C = AnalysisResults.ExampleTrials.T122C.T;
    F_C = AnalysisResults.ExampleTrials.T122C.F;
    cortical_LHnormS_C = AnalysisResults.ExampleTrials.T122C.cortical_LHnormS;
    hippocampusNormS_C = AnalysisResults.ExampleTrials.T122C.hippocampusNormS;
else
    animalID = 'T122';
    dataLocation = [rootFolder '\' animalID '\Bilateral Imaging\'];
    cd(dataLocation)
    exampleBaselineFileID = 'T122_RestingBaselines.mat';
    load(exampleBaselineFileID,'-mat')
    % example B
    exampleProcDataFileID_B = 'T122_200218_11_42_04_ProcData.mat';
    load(exampleProcDataFileID_B,'-mat')
    exampleSpecDataFileID_B = 'T122_200218_11_42_04_SpecDataA.mat';
    load(exampleSpecDataFileID_B,'-mat')
    [~,fileDate_B,~] = GetFileInfo_HRF2020(exampleProcDataFileID_B);
    strDay_B = ConvertDate_HRF2020(fileDate_B);
    dsFs_B = ProcData.notes.dsFs;
    % setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
    [z2,p2,k2] = butter(4,0.5/(dsFs_B/2),'low');
    [sos2,g2] = zp2sos(z2,p2,k2);
    % whisker angle
    binWhiskerAngle_B = ProcData.data.binWhiskerAngle;
    % force sensor
    binForceSensor_B = ProcData.data.binForceSensor;
    % HbT data (LH)
    HbT_B = ProcData.data.CBV_HbT.adjLH;
    % neural data
    gammaPower_B = (ProcData.data.cortical_LH.gammaBandPower - RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay_B))./RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay_B);
    muaPower_B = (ProcData.data.cortical_LH.muaPower - RestingBaselines.manualSelection.cortical_LH.muaPower.(strDay_B))./RestingBaselines.manualSelection.cortical_LH.muaPower.(strDay_B);
    % process gamma-band neural data
    gammaTemplate_B = zeros(size(gammaPower_B));
    strt_B = 2*dsFs_B;
    stp_B = size(gammaTemplate_B,2);
    gammaTemplate_B(:,strt_B:stp_B) = gammaPower_B(:,strt_B:stp_B) - mean(gammaPower_B(:,strt_B:stp_B));
    % process MUA neural data
    muaTemplate_B = zeros(size(muaPower_B));
    muaTemplate_B(:,strt_B:stp_B) = muaPower_B(:,strt_B:stp_B) - mean(muaPower_B(:,strt_B:stp_B));
    % process HbT data
    HbTTemplate_B = zeros(size(HbT_B));
    offset_B = mean(HbT_B)*ones(1,stp_B - strt_B + 1);
    HbTTemplate_B(:,strt_B:stp_B) = detrend(HbT_B(:,strt_B:stp_B) - offset_B);
    % gamma kernel predictions
    [gammaAct_B,gammaPred_B] = ConvolveHRF_HRF2020(gammaKernel,detrend(gammaTemplate_B),detrend(HbTTemplate_B),0);
    gamma_mPred_B = filtfilt(sos2,g2,(gammaPred_B(strt_B:stp_B) - mean(gammaPred_B(strt_B:stp_B))));
    gamma_mAct_B = filtfilt(sos2,g2,(gammaAct_B(strt_B:stp_B) - mean(gammaAct_B(strt_B:stp_B))));
    % MUA kernel predictions
    [muaAct_B,muaPred_B] = ConvolveHRF_HRF2020(muaKernel,detrend(muaTemplate_B),detrend(HbTTemplate_B),0);
    mua_mPred_B = filtfilt(sos2,g2,(muaPred_B(strt_B:stp_B) - mean(muaPred_B(strt_B:stp_B))));
    mua_mAct_B = filtfilt(sos2,g2,(muaAct_B(strt_B:stp_B) - mean(muaAct_B(strt_B:stp_B))));
    % gamma R and R2
    gammaR_B = CalculateR_HRF2020(gamma_mPred_B,gamma_mAct_B);
    gammaR2_B = CalculateRsquared_HRF2020(gamma_mPred_B,gamma_mAct_B);
    % MUA R and R2
    muaR_B = CalculateR_HRF2020(mua_mPred_B,mua_mAct_B);
    muaR2_B = CalculateRsquared_HRF2020(mua_mPred_B,mua_mAct_B);
    % pad data
    actHbT_B = [zeros(1,2*dsFs_B - 1),gamma_mAct_B];
    gammaPredHbT_B = [zeros(1,2*dsFs_B - 1),gamma_mPred_B];
    muaPredHbT_B = [zeros(1,2*dsFs_B - 1),mua_mPred_B];
    % cortical and hippocampal spectrograms
    cortical_LHnormS_B = SpecData.cortical_LH.normS.*100;
    hippocampusNormS_B = SpecData.hippocampus.normS.*100;
    T_B = SpecData.cortical_LH.T;
    F_B = SpecData.cortical_LH.F;
    % update analysis structure
    AnalysisResults.ExampleTrials.T122B.dsFs = dsFs_B;
    AnalysisResults.ExampleTrials.T122B.binForceSensor = binForceSensor_B;
    AnalysisResults.ExampleTrials.T122B.binWhiskerAngle = binWhiskerAngle_B;
    AnalysisResults.ExampleTrials.T122B.actHbT = actHbT_B;
    AnalysisResults.ExampleTrials.T122B.gammaPredHbT = gammaPredHbT_B;
    AnalysisResults.ExampleTrials.T122B.muaPredHbT = muaPredHbT_B;
    AnalysisResults.ExampleTrials.T122B.gammaR = gammaR_B;
    AnalysisResults.ExampleTrials.T122B.gammaR2 = gammaR2_B;
    AnalysisResults.ExampleTrials.T122B.muaR = muaR_B;
    AnalysisResults.ExampleTrials.T122B.muaR2 = muaR2_B;
    AnalysisResults.ExampleTrials.T122B.T = T_B;
    AnalysisResults.ExampleTrials.T122B.F = F_B;
    AnalysisResults.ExampleTrials.T122B.cortical_LHnormS = cortical_LHnormS_B;
    AnalysisResults.ExampleTrials.T122B.hippocampusNormS = hippocampusNormS_B;
    % example C
    clear ProcData
    exampleProcDataFileID_C = 'T122_200203_13_51_00_ProcData.mat';
    load(exampleProcDataFileID_C,'-mat')
    exampleSpecDataFileID_C = 'T122_200203_13_51_00_SpecDataA.mat';
    load(exampleSpecDataFileID_C,'-mat')
    [~,fileDate_C,~] = GetFileInfo_HRF2020(exampleProcDataFileID_C);
    strDay_C = ConvertDate_HRF2020(fileDate_C);
    dsFs_C = ProcData.notes.dsFs;
    % setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
    [z2,p2,k2] = butter(4,0.5/(dsFs_C/2),'low');
    [sos2,g2] = zp2sos(z2,p2,k2);
    % whisker angle
    binWhiskerAngle_C = ProcData.data.binWhiskerAngle;
    % force sensor
    binForceSensor_C = ProcData.data.binForceSensor;
    % HbT data (LH)
    HbT_C = ProcData.data.CBV_HbT.adjLH;
    % neural data
    gammaPower_C = (ProcData.data.cortical_LH.gammaBandPower - RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay_C))./RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay_C);
    muaPower_C = (ProcData.data.cortical_LH.muaPower - RestingBaselines.manualSelection.cortical_LH.muaPower.(strDay_C))./RestingBaselines.manualSelection.cortical_LH.muaPower.(strDay_C);
    % process gamma-band neural data
    gammaTemplate_C = zeros(size(gammaPower_C));
    strt_C = 2*dsFs_C;
    stp_C = size(gammaTemplate_C,2);
    gammaTemplate_C(:,strt_C:stp_C) = gammaPower_C(:,strt_C:stp_C) - mean(gammaPower_C(:,strt_C:stp_C));
    % process MUA neural data
    muaTemplate_C = zeros(size(muaPower_C));
    muaTemplate_C(:,strt_C:stp_C) = muaPower_C(:,strt_C:stp_C) - mean(muaPower_C(:,strt_C:stp_C));
    % process HbT data
    HbTTemplate_C = zeros(size(HbT_C));
    offset_C = mean(HbT_C)*ones(1,stp_C - strt_C + 1);
    HbTTemplate_C(:,strt_C:stp_C) = detrend(HbT_C(:,strt_C:stp_C) - offset_C);
    % gamma kernel predictions
    [gammaAct_C,gammaPred_C] = ConvolveHRF_HRF2020(gammaKernel,detrend(gammaTemplate_C),detrend(HbTTemplate_C),0);
    gamma_mPred_C = filtfilt(sos2,g2,(gammaPred_C(strt_C:stp_C) - mean(gammaPred_C(strt_C:stp_C))));
    gamma_mAct_C = filtfilt(sos2,g2,(gammaAct_C(strt_C:stp_C) - mean(gammaAct_C(strt_C:stp_C))));
    % MUA kernel predictions
    [muaAct_C,muaPred_C] = ConvolveHRF_HRF2020(muaKernel,detrend(muaTemplate_C),detrend(HbTTemplate_C),0);
    mua_mPred_C = filtfilt(sos2,g2,(muaPred_C(strt_C:stp_C) - mean(muaPred_C(strt_C:stp_C))));
    mua_mAct_C = filtfilt(sos2,g2,(muaAct_C(strt_C:stp_C) - mean(muaAct_C(strt_C:stp_C))));
    % gamma R and R2
    gammaR_C = CalculateR_HRF2020(gamma_mPred_C,gamma_mAct_C);
    gammaR2_C = CalculateRsquared_HRF2020(gamma_mPred_C,gamma_mAct_C);
    % MUA R and R2
    muaR_C = CalculateR_HRF2020(mua_mPred_C,mua_mAct_C);
    muaR2_C = CalculateRsquared_HRF2020(mua_mPred_C,mua_mAct_C);
    % pad data
    actHbT_C = [zeros(1,2*dsFs_C - 1),gamma_mAct_C];
    gammaPredHbT_C = [zeros(1,2*dsFs_C - 1),gamma_mPred_C];
    muaPredHbT_C = [zeros(1,2*dsFs_C - 1),mua_mPred_C];
    % cortical and hippocampal spectrograms
    cortical_LHnormS_C = SpecData.cortical_LH.normS.*100;
    hippocampusNormS_C = SpecData.hippocampus.normS.*100;
    T_C = SpecData.cortical_LH.T;
    F_C = SpecData.cortical_LH.F;
    % update analysis structure
    AnalysisResults.ExampleTrials.T122C.dsFs = dsFs_C;
    AnalysisResults.ExampleTrials.T122C.binForceSensor = binForceSensor_C;
    AnalysisResults.ExampleTrials.T122C.binWhiskerAngle = binWhiskerAngle_C;
    AnalysisResults.ExampleTrials.T122C.actHbT = actHbT_C;
    AnalysisResults.ExampleTrials.T122C.gammaPredHbT = gammaPredHbT_C;
    AnalysisResults.ExampleTrials.T122C.muaPredHbT = muaPredHbT_C;
    AnalysisResults.ExampleTrials.T122C.gammaR = gammaR_C;
    AnalysisResults.ExampleTrials.T122C.gammaR2 = gammaR2_C;
    AnalysisResults.ExampleTrials.T122C.muaR = muaR_C;
    AnalysisResults.ExampleTrials.T122C.muaR2 = muaR2_C;
    AnalysisResults.ExampleTrials.T122C.T = T_C;
    AnalysisResults.ExampleTrials.T122C.F = F_C;
    AnalysisResults.ExampleTrials.T122C.cortical_LHnormS = cortical_LHnormS_C;
    AnalysisResults.ExampleTrials.T122C.hippocampusNormS = hippocampusNormS_C;
    % save results
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end
%% Fig. 3
summaryFigure = figure('Name','Fig3 (e-j)');
sgtitle('Figure 3 - Turner et al. HRF 2020')
%% [3] IOS awake/sleep prediction examples
% HbT and behavioral indeces
ax1 = subplot(4,2,[1,3]);
whisking_Yvals_B = 35*ones(size(binWhiskerAngle_B));
force_Yvals_B = 33*ones(size(binForceSensor_B));
whiskInds_B = binWhiskerAngle_B.*whisking_Yvals_B;
forceInds_B = binForceSensor_B.*force_Yvals_B;
for x = 1:length(whiskInds_B)
    % set whisk indeces
    if whiskInds_B(1,x) == 0
        whiskInds_B(1,x) = NaN;
    end
    % set force indeces
    if forceInds_B(1,x) == 0
        forceInds_B(1,x) = NaN;
    end
end
p1 = plot((1:length(actHbT_B))/dsFs_B,actHbT_B,'color',colors_HRF2020('dark candy apple red'),'LineWidth',1);
hold on
p2 = plot((1:length(gammaPredHbT_B))/dsFs_B,gammaPredHbT_B,'color',colors_HRF2020('dark cyan'),'LineWidth',1);
p3 = plot((1:length(muaPredHbT_B))/dsFs_B,muaPredHbT_B,'color',colors_HRF2020('rich black'),'LineWidth',1);
s1 = scatter((1:length(binWhiskerAngle_B))/dsFs_B,whiskInds_B,'.','MarkerEdgeColor',colors_HRF2020('deep carrot orange'));
s2 = scatter((1:length(binForceSensor_B))/dsFs_B,forceInds_B,'.','MarkerEdgeColor',colors_HRF2020('sapphire'));
ylabel('\Delta[HbT] (\muM)')
legend([p1,p2,p3,s1,s2],'Actual','Gamma pred','MUA pred','Whisking','Movement')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([250,310,370,430,490,550,610,670,730,790,850])
xlim([250,850])
ylim([-20,35])
ax1.TickLength = [0.01,0.01];
% left cortical electrode spectrogram
ax2 = subplot(4,2,5);
SemilogImageSC_HRF2020(T_B,F_B,cortical_LHnormS_B,'y')
axis xy
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([250,310,370,430,490,550,610,670,730,790,850])
xlim([250,850])
ax2.TickLength = [0.01,0.01];
% hippocampal electrode spectrogram
ax3 = subplot(4,2,7);
SemilogImageSC_HRF2020(T_B,F_B,hippocampusNormS_B,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([250,310,370,430,490,550,610,670,730,790,850])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([250,850])
ax3.TickLength = [0.01,0.01];
% axes properties
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
% HbT and behavioral indeces
ax4 = subplot(4,2,[2,4]);
whisking_Yvals_C = 70*ones(size(binWhiskerAngle_C));
force_Yvals_C = 65*ones(size(binForceSensor_C));
whiskInds_C = binWhiskerAngle_C.*whisking_Yvals_C;
forceInds_C = binForceSensor_C.*force_Yvals_C;
for x = 1:length(whiskInds_C)
    % set whisk indeces
    if whiskInds_C(1,x) == 0
        whiskInds_C(1,x) = NaN;
    end
    % set force indeces
    if forceInds_C(1,x) == 0
        forceInds_C(1,x) = NaN;
    end
end
plot((1:length(actHbT_C))/dsFs_C,actHbT_C,'color',colors_HRF2020('dark candy apple red'),'LineWidth',1);
hold on
plot((1:length(gammaPredHbT_C))/dsFs_C,gammaPredHbT_C,'color',colors_HRF2020('dark cyan'),'LineWidth',1);
plot((1:length(muaPredHbT_C))/dsFs_C,muaPredHbT_C,'color',colors_HRF2020('rich black'),'LineWidth',1);
scatter((1:length(binWhiskerAngle_C))/dsFs_C,whiskInds_C,'.','MarkerEdgeColor',colors_HRF2020('deep carrot orange'));
scatter((1:length(binForceSensor_C))/dsFs_C,forceInds_C,'.','MarkerEdgeColor',colors_HRF2020('sapphire'));
ylabel('\Delta[HbT] (\muM)')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730])
ylim([-50,70])
ax4.TickLength = [0.01,0.01];
% left cortical electrode spectrogram
ax5 = subplot(4,2,6);
SemilogImageSC_HRF2020(T_C,F_C,cortical_LHnormS_C,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730])
ax5.TickLength = [0.01,0.01];
% hippocampal electrode spectrogram
ax6 = subplot(4,2,8);
SemilogImageSC_HRF2020(T_C,F_C,hippocampusNormS_C,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([130,730])
ax6.TickLength = [0.01,0.01];
% axes properties
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax5Pos(3:4) = ax4Pos(3:4);
ax6Pos(3:4) = ax4Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
% %% save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure,[dirpath 'Fig1']);
%     % remove surface subplots because they take forever to render
%     cla(ax5);
%     set(ax5,'YLim',[1,99]);
%     cla(ax3);
%     set(ax3,'YLim',[1,99]);
%     set(summaryFigure,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'Fig1'])
%     close(summaryFigure)
%     %% subplot figures
%     subplotImgs = figure;
%     % example 1 LH cortical LFP
%     subplot(2,1,1);
%     SemilogImageSC_HRF2020(T_B,F_B,cortical_LHnormS_B,'y')
%     caxis([-100,100])
%     set(gca,'box','off')
%     axis xy
%     axis tight
%     axis off
%     xlim([30,630])
%     % example 1 hippocampal LFP
%     subplot(2,1,2);
%     SemilogImageSC_HRF2020(T_B,F_B,hippocampusNormS_B,'y')
%     caxis([-100,100])
%     set(gca,'box','off')
%     axis xy
%     axis tight
%     axis off
%     xlim([30,630])
%     print('-painters','-dtiffn',[dirpath 'Fig1_SpecImages'])
%     close(subplotImgs)
%     %% Fig. 1
%     figure('Name','Fig1 (e-j)');
%     sgtitle('Figure Panel 1 - Turner et al. 2020')
%     %% [1e-j] IOS sleep example
%     % EMG and force sensor
%     ax1 = subplot(6,1,1);
%     p1 = plot((1:length(filtEMG_B))/dsFs_B,filtEMG_B,'color',colors_HRF2020('rich black'),'LineWidth',0.5);
%     ylabel({'EMG','power (a.u.)'})
%     ylim([-2,2.5])
%     yyaxis right
%     p2 = plot((1:length(binForceSensor_B))/dsFs_B,binForceSensor_B,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
%     ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
%     legend([p1,p2],'EMG','Pressure')
%     set(gca,'Xticklabel',[])
%     set(gca,'box','off')
%     xticks([30,90,150,210,270,330,390,450,510,570,630])
%     xlim([30,630])
%     ylim([-0.1,1])
%     ax1.TickLength = [0.01,0.01];
%     ax1.YAxis(1).Color = colors_HRF2020('rich black');
%     ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
%     % whisker angle and heart rate
%     ax2 = subplot(6,1,2);
%     p3 = plot((1:length(binWhiskerAngle_B))/dsFs_B,-binWhiskerAngle_B,'color',colors_HRF2020('rich black'),'LineWidth',0.5);
%     ylabel({'Whisker','angle (deg)'})
%     xlim([30,630])
%     ylim([-10,50])
%     yyaxis right
%     p4 = plot((1:length(heartRate_B)),heartRate_B,'color',colors_HRF2020('deep carrot orange'),'LineWidth',0.5);
%     ylabel({'Heart rate','Freq (Hz)'},'rotation',-90,'VerticalAlignment','bottom')
%     legend([p3,p4],'Whisker angle','Heart rate')
%     set(gca,'Xticklabel',[])
%     set(gca,'box','off')
%     xticks([30,90,150,210,270,330,390,450,510,570,630])
%     xlim([30,630])
%     ylim([5,10])
%     ax2.TickLength = [0.01,0.01];
%     ax2.YAxis(1).Color = colors_HRF2020('rich black');
%     ax2.YAxis(2).Color = colors_HRF2020('deep carrot orange');
%     % HbT and behavioral indeces
%     ax1 =subplot(6,1,[3,4]);
%     p5 = plot((1:length(filt_HbT_B))/dsFs_B,filt_HbT_B,'color',colors_HRF2020('sapphire'),'LineWidth',1);
%     hold on
%     x1 = xline(30,'color',colorRfcNREM,'LineWidth',2);
%     x2 = xline(73,'color',colorRfcREM,'LineWidth',2);
%     x3 = xline(200,'color',colorRfcAwake,'LineWidth',2);
%     xline(227,'color',colorRfcNREM,'LineWidth',2);
%     xline(299,'color',colorRfcAwake,'LineWidth',2);
%     xline(313,'color',colorRfcNREM,'LineWidth',2);
%     xline(390,'color',colorRfcAwake,'LineWidth',2);
%     xline(443,'color',colorRfcNREM,'LineWidth',2);
%     xline(597,'color',colorRfcAwake,'LineWidth',2);
%     ylabel('\Delta[HbT] (\muM)')
%     legend([p5,x3,x1,x2],'Barrels','Awake','NREM','REM')
%     set(gca,'TickLength',[0,0])
%     set(gca,'Xticklabel',[])
%     set(gca,'box','off')
%     xticks([30,90,150,210,270,330,390,450,510,570,630])
%     xlim([30,630])
%     ylim([-35,135])
%     ax1.TickLength = [0.01,0.01];
%     % left cortical electrode spectrogram
%     ax5 = subplot(6,1,5);
%     SemilogImageSC_HRF2020(T_B,F_B,cortical_LHnormS_B,'y')
%     axis xy
%     c2 = colorbar;
%     ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%     caxis([-100,100])
%     ylabel({'LH cort LFP','Freq (Hz)'})
%     set(gca,'Yticklabel','10^1')
%     set(gca,'Xticklabel',[])
%     set(gca,'box','off')
%     xticks([30,90,150,210,270,330,390,450,510,570,630])
%     xlim([30,630])
%     ax5.TickLength = [0.01,0.01];
%     % hippocampal electrode spectrogram
%     ax3 = subplot(6,1,6);
%     SemilogImageSC_HRF2020(T_B,F_B,hippocampusNormS_B,'y')
%     axis xy
%     c3 = colorbar;
%     ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%     caxis([-100,100])
%     xlabel('Time (min)')
%     ylabel({'Hipp LFP','Freq (Hz)'})
%     set(gca,'box','off')
%     xticks([30,90,150,210,270,330,390,450,510,570,630])
%     xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
%     xlim([30,630])
%     ax3.TickLength = [0.01,0.01];
%     % axes properties
%     ax1Pos = get(ax1,'position');
%     ax5Pos = get(ax5,'position');
%     ax6Pos = get(ax3,'position');
%     ax5Pos(3:4) = ax1Pos(3:4);
%     ax6Pos(3:4) = ax1Pos(3:4);
%     set(ax5,'position',ax5Pos);
%     set(ax3,'position',ax6Pos);
% end

end
