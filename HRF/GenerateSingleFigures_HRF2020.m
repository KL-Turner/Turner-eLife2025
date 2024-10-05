function [AnalysisResults] = GenerateSingleFigures_HRF2020(procDataFileID,hemDataType,RestingBaselines,rootFolder,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: 
%________________________________________________________________________________________________________________________

%% information and data for example
load(procDataFileID,'-mat')
[animalID,fileDate,fileID] = GetFileInfo_HRF2020(procDataFileID);
specDataFileID = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFileID,'-mat')
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
% heart rate
HbT = ProcData.data.CBV_HbT.(hemDataType);
% cortical and hippocampal spectrograms
cortical_normS = SpecData.(['cortical_' hemDataType(4:end)]).normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.(['cortical_' hemDataType(4:end)]).T;
F = SpecData.(['cortical_' hemDataType(4:end)]).F;
% solenoids
LPadSol = ProcData.data.solenoids.LPadSol;
RPadSol = ProcData.data.solenoids.RPadSol;
AudSol = ProcData.data.solenoids.AudSol;
% stim indeces
LPad_Yvals = 1.30*max(HbT)*ones(size(LPadSol));
RPad_Yvals = 1.30*max(HbT)*ones(size(RPadSol));
Aud_Yvals = 1.30*max(HbT)*ones(size(AudSol));
%% HRF predictions
% kernels
gammaKernel = AnalysisResults.HRFs.(animalID).gammaBandPower.(['cort' hemDataType(4:end)]).All.FM_gammaFunction;
muaKernel = AnalysisResults.HRFs.(animalID).muaPower.(['cort' hemDataType(4:end)]).All.FM_gammaFunction;
% neural data
gammaPower = (ProcData.data.(['cortical_' hemDataType(4:end)]).gammaBandPower - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:end)]).gammaBandPower.(strDay))./RestingBaselines.manualSelection.(['cortical_' hemDataType(4:end)]).gammaBandPower.(strDay);
muaPower = (ProcData.data.(['cortical_' hemDataType(4:end)]).muaPower - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:end)]).muaPower.(strDay))./RestingBaselines.manualSelection.(['cortical_' hemDataType(4:end)]).muaPower.(strDay);
% process gamma-band neural data
gammaTemplate = zeros(size(gammaPower));
strt = 2*dsFs;
stp = size(gammaTemplate,2);
gammaTemplate(:,strt:stp) = gammaPower(:,strt:stp) - mean(gammaPower(:,strt:stp));
% process MUA neural data
muaTemplate = zeros(size(muaPower));
muaTemplate(:,strt:stp) = muaPower(:,strt:stp) - mean(muaPower(:,strt:stp));
% process HbT data
HbTTemplate = zeros(size(HbT));
offset = mean(HbT)*ones(1,stp - strt + 1);
HbTTemplate(:,strt:stp) = detrend(HbT(:,strt:stp) - offset);
% gamma kernel predictions
[gammaAct,gammaPred] = ConvolveHRF_HRF2020(gammaKernel,detrend(gammaTemplate),detrend(HbTTemplate),0);
gamma_mPred = filtfilt(sos2,g2,(gammaPred(strt:stp) - mean(gammaPred(strt:stp))));
gamma_mAct = filtfilt(sos2,g2,(gammaAct(strt:stp) - mean(gammaAct(strt:stp))));
% MUA kernel predictions
[muaAct,muaPred] = ConvolveHRF_HRF2020(muaKernel,detrend(muaTemplate),detrend(HbTTemplate),0);
mua_mPred = filtfilt(sos2,g2,(muaPred(strt:stp) - mean(muaPred(strt:stp))));
mua_mAct = filtfilt(sos2,g2,(muaAct(strt:stp) - mean(muaAct(strt:stp))));
% gamma R and R2
gammaR = CalculateR_HRF2020(gamma_mPred,gamma_mAct);
gammaR2 = CalculateRsquared_HRF2020(gamma_mPred,gamma_mAct);
% MUA R and R2
muaR = CalculateR_HRF2020(mua_mPred,mua_mAct);
muaR2 = CalculateRsquared_HRF2020(mua_mPred,mua_mAct);
% pad data 
actHbT = [zeros(1,2*dsFs - 1),gamma_mAct];
gammaPredHbT = [zeros(1,2*dsFs - 1),gamma_mPred];
muaPredHbT = [zeros(1,2*dsFs - 1),mua_mPred];
%% Fig.
summaryFigure = figure;
sgtitle([animalID ' ' strrep(fileID,'_',' ') ' ' hemDataType])
% EMG and force sensor
ax1 = subplot(6,1,1);
p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors_HRF2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','power (a.u.)'})
yyaxis right
p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','Pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xlim([0,900])
ax1.YAxis(1).Color = colors_HRF2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% whisker angle and heart rate
ax2 = subplot(6,1,2);
p3 = plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors_HRF2020('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
xlim([0,900])
yyaxis right
p4 = plot((1:length(heartRate)),heartRate,'color',colors_HRF2020('deep carrot orange'),'LineWidth',0.5);
ylabel({'Heart rate','Freq (Hz)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p3,p4],'Whisker angle','Heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xlim([0,900])
ax2.YAxis(1).Color = colors_HRF2020('rich black');
ax2.YAxis(2).Color = colors_HRF2020('deep carrot orange');
% HbT and behavioral indeces
ax34 = subplot(6,1,[3,4]);
p5 = plot((1:length(actHbT))/dsFs,actHbT,'color',colors_HRF2020('dark candy apple red'),'LineWidth',1);
hold on
p6 = plot((1:length(gammaPredHbT))/dsFs,gammaPredHbT,'color',colors_HRF2020('sapphire'),'LineWidth',1);
p7 = plot((1:length(muaPredHbT))/dsFs,muaPredHbT,'color',colors_HRF2020('rich black'),'LineWidth',1);
s1 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s2 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s3 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
ylabel('\Delta[HbT] (\muM)')
title(['Gamma R = ' num2str(round(gammaR,2)) '      Gamma R^2 = '  num2str(round(gammaR2,2)) '      MUA R = '  num2str(round(muaR,2)) '      MUA R^2 = '  num2str(round(muaR2,2))])
legend([p5,p6,p7,s1,s2,s3],'HbT','Gamma-band','MUA kernel','LPad','RPad','Aud')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xlim([0,900])
% cortical electrode spectrogram
ax5 = subplot(6,1,5);
SemilogImageSC_HRF2020(T,F,cortical_normS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'LH cort LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xlim([0,900])
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
xlim([0,900])
% axes properties
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
%% save figure(s)
dirpath = [rootFolder delim animalID delim 'Figures' delim 'Full Trial Predictions' delim];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath animalID '_' fileID '_' hemDataType '_SingleTrialPrediction']);
close(summaryFigure)

end
