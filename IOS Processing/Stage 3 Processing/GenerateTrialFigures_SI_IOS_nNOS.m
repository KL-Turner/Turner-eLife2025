function [figHandle] = GenerateTrialFigures_SI_IOS_nNOS(procDataFileID,RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Create a summary figure for a single n minute IOS trial
%________________________________________________________________________________________________________________________

% load file and gather information
load(procDataFileID)
[animalID,fileDate,fileID] = GetFileInfo_IOS_nNOS(procDataFileID);
strDay = ConvertDate_IOS_nNOS(fileDate);
fileID2 = strrep(fileID,'_',' ');
% imaging wavelengths
imagingWavelengths = ProcData.notes.imagingWavelengths;
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,1/(ProcData.notes.wavelengthSamplingRate/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filteredWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
binWhiskers = ProcData.data.binWhiskerAngle;
% force sensor
filtForceSensor = filtfilt(sos1,g1,ProcData.data.forceSensor);
binForce = ProcData.data.binForceSensor;
% EMG
EMG = ProcData.data.EMG.emg;
% heart rate
heartRate = ProcData.data.heartRate;
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
OptoLED = ProcData.data.stimulations.OptoLED;
% IOS data
if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Lime, Green, & Blue'})) == true
    % HbT
    barrels_HbT = ProcData.data.CBV_HbT.Barrels;
    filtBarrels_HbT = filtfilt(sos2,g2,barrels_HbT);
    % GCaMP
    barrels_GCaMP7s = ProcData.data.GCaMP7s.corBarrels;
    normBarrels_GCaMP7s = (barrels_GCaMP7s - 1)*100;
    filtBarrels_GCaMP7s = filtfilt(sos2,g2,normBarrels_GCaMP7s);
    % Deoxy signal
    barrels_deoxy = ProcData.data.Deoxy.LH;
    normBarrels_Deoxy = ((barrels_deoxy - RestingBaselines.manualSelection.Deoxy.LH.(strDay).mean)./RestingBaselines.manualSelection.Deoxy.LH.(strDay).mean)*100;
    filtBarrels_Deoxy = filtfilt(sos2,g2,normBarrels_Deoxy);
elseif any(strcmp(imagingWavelengths,{'Green & Blue','Lime & Blue'})) == true
    % HbT
    barrels_HbT = ProcData.data.CBV_HbT.Barrels;
    filtBarrels_HbT = filtfilt(sos2,g2,barrels_HbT);
    % GCaMP
    barrels_GCaMP7s = ProcData.data.GCaMP7s.corBarrels;
    normBarrels_GCaMP7s = (barrels_GCaMP7s - 1)*100;
    filtBarrels_GCaMP7s = filtfilt(sos2,g2,normBarrels_GCaMP7s);
elseif any(strcmp(imagingWavelengths,{'Green','Lime'})) == true
    % HbT
    barrels_HbT = ProcData.data.CBV_HbT.Barrels;
    filtBarrels_HbT = filtfilt(sos2,g2,barrels_HbT);
elseif strcmp(imagingWavelengths,'Blue') == true
    % dR/R
    barrels_CBV = ProcData.data.CBV.Barrels;
    normBarrels_CBV = (barrels_CBV - RestingBaselines.manualSelection.CBV.Barrels.(strDay).mean)./(RestingBaselines.manualSelection.CBV.Barrels.(strDay).mean);
    filtBarrels_CBV = (filtfilt(sos2,g2,normBarrels_CBV))*100;
end
% cortical and hippocampal spectrograms
specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFile,'-mat');
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
cortical_RHnormS = SpecData.cortical_RH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
% Yvals for behavior Indices
if strcmp(imagingWavelengths,'Blue') == true
    whisking_Yvals = 1.10*max(filtBarrels_CBV)*ones(size(binWhiskers));
    force_Yvals = 1.20*max(filtBarrels_CBV)*ones(size(binForce));
    LPad_Yvals = 1.30*max(filtBarrels_CBV)*ones(size(LPadSol));
    RPad_Yvals = 1.30*max(filtBarrels_CBV)*ones(size(RPadSol));
    Aud_Yvals = 1.30*max(filtBarrels_CBV)*ones(size(AudSol));
    Opto_Yvals = 1.30*max(filtBarrels_CBV)*ones(size(OptoLED));
else
    whisking_Yvals = 1.10*max(filtBarrels_HbT)*ones(size(binWhiskers));
    force_Yvals = 1.20*max(filtBarrels_HbT)*ones(size(binForce));
    LPad_Yvals = 1.30*max(filtBarrels_HbT)*ones(size(LPadSol));
    RPad_Yvals = 1.30*max(filtBarrels_HbT)*ones(size(RPadSol));
    Aud_Yvals = 1.30*max(filtBarrels_HbT)*ones(size(AudSol));
    Opto_Yvals = 1.30*max(filtBarrels_HbT)*ones(size(OptoLED));
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
% Figure
if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Lime, Green, & Blue'})) == true
    figHandle = figure;
    % force sensor and EMG
    ax1 = subplot(8,1,1);
    p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('sapphire'),'LineWidth',1);
    title([animalID ' ' fileID2])
    ylabel('Force Sensor (Volts)')
    xlim([0,ProcData.notes.trialDuration_sec])
    yyaxis right
    p2 = plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color',colors('deep carrot orange'),'LineWidth',1);
    ylabel('EMG (Volts^2)')
    xlim([0,ProcData.notes.trialDuration_sec])
    legend([p1,p2],'force sensor','EMG')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % whisker angle and heart rate
    ax2 = subplot(8,1,2);
    p3 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('blue-green'),'LineWidth',1);
    ylabel('Angle (deg)')
    xlim([0,ProcData.notes.trialDuration_sec])
    ylim([-20,60])
    yyaxis right
    p4 = plot((1:length(heartRate)),heartRate,'color',colors('dark sea green'),'LineWidth',1);
    ylabel('Heart Rate (Hz)')
    ylim([6,15])
    legend([p3,p4],'whisker angle','heart rate')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % CBV and behavioral indeces
    ax3 = subplot(8,1,3);
    s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('sapphire'));
    hold on
    s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('blue-green'));
    s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
    s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
    s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
    s6 = scatter(OptoLED,Opto_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','b');
    p5 = plot((1:length(filtBarrels_HbT))/ProcData.notes.wavelengthSamplingRate,filtBarrels_HbT,'color',colors('dark candy apple red'),'LineWidth',1);
    legend([p5,s1,s2,s3,s4,s5,s6],'Barrels','movement','whisking',',LPad sol','RPad sol','Aud sol','Opto LED')
    ylabel('\DeltaHbT')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % GCaMP
    ax4 = subplot(8,1,4);
    plot((1:length(filtBarrels_GCaMP7s))/ProcData.notes.wavelengthSamplingRate,filtBarrels_GCaMP7s,'color',colors('dark candy apple red'),'LineWidth',1);
    ylabel('GCaMP7s /DeltaF/F')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % Deoxy
    ax5 = subplot(8,1,5);
    plot((1:length(filtBarrels_Deoxy))/ProcData.notes.wavelengthSamplingRate,filtBarrels_Deoxy,'color',colors('dark candy apple red'),'LineWidth',1);
    ylabel('Deoxy /DeltaR/R')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % left cortical electrode spectrogram
    ax6 = subplot(8,1,6);
    Semilog_ImageSC(T,F,cortical_LHnormS,'y')
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
    ylabel('Left cortical LFP')
    set(gca,'Yticklabel',[])
    % right cortical electrode spectrogram
    ax7 = subplot(8,1,7);
    Semilog_ImageSC(T,F,cortical_RHnormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)')
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
    % hippocampal electrode spectrogram
    ax8 = subplot(8,1,8);
    Semilog_ImageSC(T,F,hippocampusNormS,'y')
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)')
    caxis([-100,100])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'box','off')
    yyaxis right
    ylabel('Hippocampal LFP')
    set(gca,'Yticklabel',[])
    % axes properties
    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],'x')
    ax1Pos = get(ax1,'position');
    ax6Pos = get(ax6,'position');
    ax7Pos = get(ax7,'position');
    ax8Pos = get(ax8,'position');
    ax6Pos(3:4) = ax1Pos(3:4);
    ax7Pos(3:4) = ax1Pos(3:4);
    ax8Pos(3:4) = ax1Pos(3:4);
    set(ax6,'position',ax6Pos);
    set(ax7,'position',ax7Pos);
    set(ax8,'position',ax8Pos);
elseif any(strcmp(imagingWavelengths,{'Green & Blue','Lime & Blue'})) == true
    figHandle = figure;
    % force sensor and EMG
    ax1 = subplot(7,1,1);
    p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('sapphire'),'LineWidth',1);
    title([animalID ' ' fileID2])
    ylabel('Force Sensor (Volts)')
    xlim([0,ProcData.notes.trialDuration_sec])
    yyaxis right
    p2 = plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color',colors('deep carrot orange'),'LineWidth',1);
    ylabel('EMG (Volts^2)')
    xlim([0,ProcData.notes.trialDuration_sec])
    legend([p1,p2],'force sensor','EMG')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % whisker angle and heart rate
    ax2 = subplot(7,1,2);
    p3 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('blue-green'),'LineWidth',1);
    ylabel('Angle (deg)')
    xlim([0,ProcData.notes.trialDuration_sec])
    ylim([-20,60])
    yyaxis right
    p4 = plot((1:length(heartRate)),heartRate,'color',colors('dark sea green'),'LineWidth',1);
    ylabel('Heart Rate (Hz)')
    ylim([6,15])
    legend([p3,p4],'whisker angle','heart rate')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % CBV and behavioral indeces
    ax3 = subplot(7,1,3);
    s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('sapphire'));
    hold on
    s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('blue-green'));
    s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
    s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
    s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
    s6 = scatter(OptoLED,Opto_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','b');
    p5 = plot((1:length(filtBarrels_HbT))/ProcData.notes.wavelengthSamplingRate,filtBarrels_HbT,'color',colors('dark candy apple red'),'LineWidth',1);
    legend([p5,s1,s2,s3,s4,s5,s6],'Barrels','movement','whisking',',LPad sol','RPad sol','Aud sol','Opto LED')
    ylabel('\DeltaHbT')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % GCaMP
    ax4 = subplot(7,1,4);
    plot((1:length(filtBarrels_GCaMP7s))/ProcData.notes.wavelengthSamplingRate,filtBarrels_GCaMP7s,'color',colors('dark candy apple red'),'LineWidth',1);
    ylabel('GCaMP7s /DeltaF/F')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % left cortical electrode spectrogram
    ax5 = subplot(7,1,5);
    Semilog_ImageSC(T,F,cortical_LHnormS,'y')
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
    ylabel('Left cortical LFP')
    set(gca,'Yticklabel',[])
    % right cortical electrode spectrogram
    ax6 = subplot(7,1,6);
    Semilog_ImageSC(T,F,cortical_RHnormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)')
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
    % hippocampal electrode spectrogram
    ax7 = subplot(7,1,7);
    Semilog_ImageSC(T,F,hippocampusNormS,'y')
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)')
    caxis([-100,100])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'box','off')
    yyaxis right
    ylabel('Hippocampal LFP')
    set(gca,'Yticklabel',[])
    % axes properties
    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],'x')
    ax1Pos = get(ax1,'position');
    ax6Pos = get(ax6,'position');
    ax7Pos = get(ax7,'position');
    ax8Pos = get(ax8,'position');
    ax6Pos(3:4) = ax1Pos(3:4);
    ax7Pos(3:4) = ax1Pos(3:4);
    ax8Pos(3:4) = ax1Pos(3:4);
    set(ax6,'position',ax6Pos);
    set(ax7,'position',ax7Pos);
    set(ax8,'position',ax8Pos);
elseif any(strcmp(imagingWavelengths,{'Green','Lime'})) == true
    figHandle = figure;
    % force sensor and EMG
    ax1 = subplot(6,1,1);
    p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('sapphire'),'LineWidth',1);
    title([animalID ' ' fileID2])
    ylabel('Force Sensor (Volts)')
    xlim([0,ProcData.notes.trialDuration_sec])
    yyaxis right
    p2 = plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color',colors('deep carrot orange'),'LineWidth',1);
    ylabel('EMG (Volts^2)')
    xlim([0,ProcData.notes.trialDuration_sec])
    legend([p1,p2],'force sensor','EMG')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % whisker angle and heart rate
    ax2 = subplot(6,1,2);
    p3 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('blue-green'),'LineWidth',1);
    ylabel('Angle (deg)')
    xlim([0,ProcData.notes.trialDuration_sec])
    ylim([-20,60])
    yyaxis right
    p4 = plot((1:length(heartRate)),heartRate,'color',colors('dark sea green'),'LineWidth',1);
    ylabel('Heart Rate (Hz)')
    ylim([6,15])
    legend([p3,p4],'whisker angle','heart rate')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % CBV and behavioral indeces
    ax3 = subplot(6,1,3);
    s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('sapphire'));
    hold on
    s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('blue-green'));
    s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
    s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
    s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
    s6 = scatter(OptoLED,Opto_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','b');
    p5 = plot((1:length(filtBarrels_HbT))/ProcData.notes.wavelengthSamplingRate,filtBarrels_HbT,'color',colors('dark candy apple red'),'LineWidth',1);
    legend([p5,s1,s2,s3,s4,s5,s6],'Barrels','movement','whisking',',LPad sol','RPad sol','Aud sol','Opto LED')
    ylabel('\DeltaHbT')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % left cortical electrode spectrogram
    ax4 = subplot(6,1,4);
    Semilog_ImageSC(T,F,cortical_LHnormS,'y')
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
    ylabel('Left cortical LFP')
    set(gca,'Yticklabel',[])
    % right cortical electrode spectrogram
    ax5 = subplot(6,1,5);
    Semilog_ImageSC(T,F,cortical_RHnormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)')
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
    % hippocampal electrode spectrogram
    ax6 = subplot(6,1,6);
    Semilog_ImageSC(T,F,hippocampusNormS,'y')
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)')
    caxis([-100,100])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'box','off')
    yyaxis right
    ylabel('Hippocampal LFP')
    set(gca,'Yticklabel',[])
    % axes properties
    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x')
    ax1Pos = get(ax1,'position');
    ax4Pos = get(ax4,'position');
    ax5Pos = get(ax5,'position');
    ax6Pos = get(ax6,'position');
    ax4Pos(3:4) = ax1Pos(3:4);
    ax5Pos(3:4) = ax1Pos(3:4);
    ax6Pos(3:4) = ax1Pos(3:4);
    set(ax4,'position',ax4Pos);
    set(ax5,'position',ax5Pos);
    set(ax6,'position',ax6Pos);
elseif strcmp(imagingWavelengths,'Blue') == true
    figHandle = figure;
    % force sensor and EMG
    ax1 = subplot(6,1,1);
    p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('sapphire'),'LineWidth',1);
    title([animalID ' ' fileID2])
    ylabel('Force Sensor (Volts)')
    xlim([0,ProcData.notes.trialDuration_sec])
    yyaxis right
    p2 = plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color',colors('deep carrot orange'),'LineWidth',1);
    ylabel('EMG (Volts^2)')
    xlim([0,ProcData.notes.trialDuration_sec])
    legend([p1,p2],'force sensor','EMG')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % whisker angle and heart rate
    ax2 = subplot(6,1,2);
    p3 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('blue-green'),'LineWidth',1);
    ylabel('Angle (deg)')
    xlim([0,ProcData.notes.trialDuration_sec])
    ylim([-20,60])
    yyaxis right
    p4 = plot((1:length(heartRate)),heartRate,'color',colors('dark sea green'),'LineWidth',1);
    ylabel('Heart Rate (Hz)')
    ylim([6,15])
    legend([p3,p4],'whisker angle','heart rate')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % CBV and behavioral indeces
    ax3 = subplot(6,1,3);
    s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('sapphire'));
    hold on
    s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('blue-green'));
    s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
    s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
    s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
    s6 = scatter(OptoLED,Opto_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','b');
    p5 = plot((1:length(filtBarrels_CBV))/ProcData.notes.wavelengthSamplingRate,filtBarrels_CBV,'color',colors('dark candy apple red'),'LineWidth',1);
    legend([p5,s1,s2,s3,s4,s5,s6],'Barrels','movement','whisking',',LPad sol','RPad sol','Aud sol','Opto LED')
    ylabel('\DeltaF/F')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % left cortical electrode spectrogram
    ax4 = subplot(6,1,4);
    Semilog_ImageSC(T,F,cortical_LHnormS,'y')
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
    ylabel('Left cortical LFP')
    set(gca,'Yticklabel',[])
    % right cortical electrode spectrogram
    ax5 = subplot(6,1,5);
    Semilog_ImageSC(T,F,cortical_RHnormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)')
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
    % hippocampal electrode spectrogram
    ax6 = subplot(6,1,6);
    Semilog_ImageSC(T,F,hippocampusNormS,'y')
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)')
    caxis([-100,100])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'box','off')
    yyaxis right
    ylabel('Hippocampal LFP')
    set(gca,'Yticklabel',[])
    % axes properties
    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x')
    ax1Pos = get(ax1,'position');
    ax4Pos = get(ax4,'position');
    ax5Pos = get(ax5,'position');
    ax6Pos = get(ax6,'position');
    ax4Pos(3:4) = ax1Pos(3:4);
    ax5Pos(3:4) = ax1Pos(3:4);
    ax6Pos(3:4) = ax1Pos(3:4);
    set(ax4,'position',ax4Pos);
    set(ax5,'position',ax5Pos);
    set(ax6,'position',ax6Pos);
end

end
