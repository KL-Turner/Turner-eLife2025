function [singleTrialFig] = GenerateSingleFigures_2P(mergedDataFileID,baselineType,saveFigs,RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Create a summary figure for a single n minute two photon trial
%________________________________________________________________________________________________________________________

% load file and gather information
load(mergedDataFileID,'-mat')
[animalID,hem,fileDate,fileID,imageID,vesselID] = GetFileInfo2_2P(mergedDataFileID);
strDay = ConvertDate_2P(fileDate);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(MergedData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,1/(MergedData.notes.p2Fs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filteredWhiskerAngle = filtfilt(sos1,g1,MergedData.data.whiskerAngle);
binWhiskers = MergedData.data.binWhiskerAngle;
% force sensor
filtForceSensor = filtfilt(sos1,g1,MergedData.data.forceSensorL);
binForce = MergedData.data.binForceSensorM;
% emg
EMG = MergedData.data.EMG.data;
% solenoids
LPadSol = MergedData.data.solenoids.LPadSol;
RPadSol = MergedData.data.solenoids.RPadSol;
AudSol = MergedData.data.solenoids.AudSol;
% vessel diameter
if strcmp(vesselID(1),'A') == true || strcmp(vesselID(1),'P') == true || strcmp(vesselID(1),'V') == true
    vesselDiameter = MergedData.data.vesselDiameter.data;
    normVesselDiameter = (vesselDiameter - RestingBaselines.(baselineType).vesselDiameter.data.(vesselID).(strDay))./(RestingBaselines.(baselineType).vesselDiameter.data.(vesselID).(strDay));
    filtVesselDiameter = filtfilt(sos2,g2,normVesselDiameter)*100;
elseif strcmp(vesselID(1),'C') == true || strcmp(vesselID(1),'D') == true
    vesselDiameter = MergedData.data.vesselDiameter.data;
    filtVesselDiameter = filtfilt(sos2,g2,vesselDiameter);
end
% cortical and hippocampal spectrograms
specDataFile = [animalID '_' hem '_' fileID '_' imageID '_' vesselID '_SpecData.mat'];
load(specDataFile,'-mat');
cortNormS = SpecData.corticalNeural.fiveSec.normS.*100;
hipNormS = SpecData.hippocampalNeural.fiveSec.normS.*100;
T = SpecData.corticalNeural.fiveSec.T;
F = SpecData.corticalNeural.fiveSec.F;
% Yvals for behavior Indices
whisking_YVals = 1.10*max(filtVesselDiameter)*ones(size(binWhiskers));
force_YVals = 1.20*max(filtVesselDiameter)*ones(size(binForce));
LPad_Yvals = 1.30*max(filtVesselDiameter)*ones(size(LPadSol));
RPad_Yvals = 1.30*max(filtVesselDiameter)*ones(size(RPadSol));
Aud_Yvals = 1.30*max(filtVesselDiameter)*ones(size(AudSol));
whiskInds = binWhiskers.*whisking_YVals;
forceInds = binForce.*force_YVals;
for x = 1:length(whiskInds)
    % set whisk indeces
    if whiskInds(1,x) == 0
        whiskInds(1,x) = NaN;
    end
    % set force indeces
    if forceInds(1,x) == 0
        forceInds(1,x) = NaN;
    end
end
% Figure
singleTrialFig = figure;
fileID2 = strrep(fileID,'_',' ');
sgtitle([animalID ' Two-photon behavioral characterization and vessel ' vesselID ' dynamics for ' fileID2 ' image ' imageID])
% check XCorr
ax1 = subplot(6,1,1);
plot(MergedData.notes.MScan.shiftLags,MergedData.notes.MScan.shiftXCorr,'color',colors('rich black'),'LineWidth',1)
ylabel('XCorr (A.U.)')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% force sensor and EMG
ax2 = subplot(6,1,2);
plot((1:length(filtForceSensor))/MergedData.notes.dsFs,filtForceSensor,'color',colors('sapphire'),'LineWidth',1)
ylabel('Force (V)')
xlim([0,MergedData.notes.trialDuration_Sec])
yyaxis right
plot((1:length(EMG))/MergedData.notes.dsFs,EMG,'color',colors('deep carrot orange'),'LineWidth',1)
ylabel('EMG (Volts^2)')
xlim([0,MergedData.notes.trialDuration_Sec])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% whisker angle
ax3 = subplot(6,1,3);
plot((1:length(filteredWhiskerAngle))/MergedData.notes.dsFs,-filteredWhiskerAngle,'color',colors('rich black'),'LineWidth',1)
ylabel('Whisker angle (deg)')
xlim([0,MergedData.notes.trialDuration_Sec])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% vessel diameter
ax4 = subplot(6,1,4);
plot((1:length(filtVesselDiameter))/MergedData.notes.p2Fs,filtVesselDiameter,'color',colors('dark candy apple red'),'LineWidth',1)
hold on;
s1 = scatter((1:length(binForce))/MergedData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('sapphire'));
s2 = scatter((1:length(binWhiskers))/MergedData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('electric purple'));
s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
if strcmp(vesselID(1),'A') == true || strcmp(vesselID(1),'P') == true || strcmp(vesselID(1),'V') == true
    ylabel('\DeltaD/D (%)')
elseif strcmp(vesselID(1),'C') == true || strcmp(vesselID(1),'D') == true
    ylabel('\DeltaV (\muM/sec)')
end
legend([s1,s2,s3,s4,s5],'Movement','Whisking','LPadSol','RPadSol','AudSol')
xlim([0,MergedData.notes.trialDuration_Sec])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% cortical LFP
ax5 = subplot(6,1,5);
Semilog_ImageSC(T,F,cortNormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)')
caxis([-100,100])
ylabel({'Cortical LFP','Freq (Hz)'})
xlim([0,MergedData.notes.trialDuration_Sec])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% hippocampal LFP
ax6 = subplot(6,1,6);
Semilog_ImageSC(T,F,hipNormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)')
caxis([-100,100])
xlabel('Time (sec)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
xlim([0,MergedData.notes.trialDuration_Sec])
set(gca,'box','off')
axis tight
% Axes properties
linkaxes([ax2,ax3,ax4,ax5,ax6],'x')
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
% Save the file to directory.
if strcmp(saveFigs,'y') == true
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(singleTrialFig,[dirpath animalID '_' hem '_' fileID '_' imageID '_' vesselID '_SingleTrialFig']);
end

end
