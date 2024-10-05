function [figHandle,ax1,ax2,ax3,ax4,ax5] = GenerateSingleFigures_GCaMP_Sleep_IOS_nNOS(procDataFileID)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% load file and gather information
load(procDataFileID)
[animalID,~,fileID] = GetFileInfo_IOS_nNOS(procDataFileID);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,1/(ProcData.notes.CBVCamSamplingRate/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filteredWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
binWhiskers = ProcData.data.binWhiskerAngle;
% force sensor
filtForceSensor = filtfilt(sos1,g1,ProcData.data.forceSensor);
binForce = ProcData.data.binForceSensor;
% emg
EMG = ProcData.data.EMG.emg;
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
% CBV data
LH_HbT = ProcData.data.CBV_HbT.LH;
filtLH_HbT = filtfilt(sos2,g2,LH_HbT);
RH_HbT = ProcData.data.CBV_HbT.RH;
filtRH_HbT = filtfilt(sos2,g2,RH_HbT);
% cortical and hippocampal spectrograms
specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFile,'-mat');
cortical_RHnormS = SpecData.cortical_RH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
% Yvals for behavior Indices
if max(filtLH_HbT) >= max(filtRH_HbT)
    whisking_Yvals = 1.10*max(filtLH_HbT)*ones(size(binWhiskers));
    force_Yvals = 1.20*max(filtLH_HbT)*ones(size(binForce));
    LPad_Yvals = 1.50*max(filtLH_HbT)*ones(size(LPadSol));
    RPad_Yvals = 1.50*max(filtLH_HbT)*ones(size(RPadSol));
    Aud_Yvals = 1.50*max(filtLH_HbT)*ones(size(AudSol));
else
    whisking_Yvals = 1.10*max(filtRH_HbT)*ones(size(binWhiskers));
    force_Yvals = 1.20*max(filtRH_HbT)*ones(size(binForce));
    LPad_Yvals = 1.50*max(filtRH_HbT)*ones(size(LPadSol));
    RPad_Yvals = 1.50*max(filtRH_HbT)*ones(size(RPadSol));
    Aud_Yvals = 1.50*max(filtRH_HbT)*ones(size(AudSol));
end
forceInds = binForce.*force_Yvals;
whiskInds = binWhiskers.*whisking_Yvals;
% set force indeces
for x = 1:length(forceInds)
    if forceInds(1,x) == 0
        forceInds(1,x) = NaN;
    end
end
% set whisk indeces
for x = 1:length(whiskInds)
    if whiskInds(1,x) == 0
        whiskInds(1,x) = NaN;
    end
end
% figure
figHandle = figure;
% force sensor and EMG
ax1 = subplot(5,1,1);
fileID2 = strrep(fileID,'_',' ');
p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('electric purple'),'LineWidth',1);
title([animalID ' ' fileID2])
ylabel('Force Sensor (Volts)')
xlim([0,ProcData.notes.trialDuration_sec])
yyaxis right
p2 = plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color',colors('dark pink'),'LineWidth',1);
ylabel('EMG (Volts^2)')
xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'force sensor','EMG')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% Whisker angle and heart rate
ax2 = subplot(5,1,2);
plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('smoky black'),'LineWidth',1);
ylabel('Angle (deg)')
xlim([0,ProcData.notes.trialDuration_sec])
ylim([-20,60])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% CBV and behavioral indeces
ax3 = subplot(5,1,3);
s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('electric purple'));
hold on
s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('smoky black'));
s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
p5 = plot((1:length(filtLH_HbT))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT,'color',colors('dark candy apple red'),'LineWidth',1);
p6 = plot((1:length(filtRH_HbT))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT,'color',colors('sapphire'),'LineWidth',1);
ylabel('\DeltaHbT')
legend([p5,p6,s1,s2,s3,s4,s5],'LH HbT','RH HbT','movement','whisking',',LPad sol','RPad sol','Aud sol')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
ax4 = subplot(5,1,4);
Semilog_ImageSC(T,F,cortical_RHnormS,'y')
axis xy
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)')
caxis([-100,100])
ylabel('Frequency (Hz)')
set(gca,'Yticklabel','10^1')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
yyaxis right
ylabel('Right cortical LFP')
set(gca,'Yticklabel',[])
% Hippocampal electrode spectrogram
ax5 = subplot(5,1,5);
Semilog_ImageSC(T,F,hippocampusNormS,'y')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)')
caxis([-100,100])
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'box','off')
yyaxis right
ylabel('Hippocampal LFP')
set(gca,'Yticklabel',[])
% Axes properties
linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
ax1Pos = get(ax1,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax1Pos(3:4);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);

end
