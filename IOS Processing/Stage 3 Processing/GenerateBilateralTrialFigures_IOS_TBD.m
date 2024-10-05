function [figHandle] = GenerateBilateralTrialFigures_IOS(procDataFileID)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% load file and gather information
load(procDataFileID)
[animalID,~,fileID] = GetFileInfo_IOS(procDataFileID);
fileID2 = strrep(fileID,'_',' ');
% imaging wavelengths
imagingWavelengths = ProcData.notes.imagingWavelengths;
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,1/(ProcData.notes.wavelengthSamplingRate/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filteredWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle.angle);
binWhiskers = ProcData.data.whiskerAngle.binarization;
% force sensor
filtForceSensor = filtfilt(sos1,g1,ProcData.data.forceSensor.force);
binForce = ProcData.data.forceSensor.binarization;
% EMG
EMG = ProcData.data.EMG.power;
% heart rate
heartRate = ProcData.data.heartRate.frequency;
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
OptoLED = ProcData.data.stimulations.OptoLED;
% IOS data
if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Lime, Green, & Blue'})) == true
    % HbT
    LH_HbT = ProcData.data.HbT.LH;
    filtLH_HbT = filtfilt(sos2,g2,LH_HbT);
    RH_HbT = ProcData.data.HbT.RH;
    filtRH_HbT = filtfilt(sos2,g2,RH_HbT);
    % HbO
    LH_HbO = ProcData.data.HbO.LH;
    filtLH_HbO = filtfilt(sos2,g2,LH_HbO);
    RH_HbO = ProcData.data.HbT.RH;
    filtRH_HbO = filtfilt(sos2,g2,RH_HbO);
    % HbR
    LH_HbR = ProcData.data.HbT.LH;
    filtLH_HbR = filtfilt(sos2,g2,LH_HbR);
    RH_HbR = ProcData.data.HbT.RH;
    filtRH_HbR = filtfilt(sos2,g2,RH_HbR);
    % GCaMP
    LH_GCaMP = ProcData.data.GCaMP.LH;
    normLH_GCaMP = (LH_GCaMP - 1)*100;
    filtLH_GCaMP = filtfilt(sos2,g2,normLH_GCaMP);
    RH_GCaMP = ProcData.data.GCaMP.RH;
    normRH_GCaMP = (RH_GCaMP - 1)*100;
    filtRH_GCaMP = filtfilt(sos2,g2,normRH_GCaMP);
elseif any(strcmp(imagingWavelengths,{'Green & Blue','Lime & Blue'})) == true
    % HbT
    LH_HbT = ProcData.data.HbT.LH;
    filtLH_HbT = filtfilt(sos2,g2,LH_HbT);
    RH_HbT = ProcData.data.HbT.RH;
    filtRH_HbT = filtfilt(sos2,g2,RH_HbT);
    % GCaMP
    LH_GCaMP = ProcData.data.GCaMP.LH;
    normLH_GCaMP = (LH_GCaMP - 1)*100;
    filtLH_GCaMP = filtfilt(sos2,g2,normLH_GCaMP);
    RH_GCaMP = ProcData.data.GCaMP.RH;
    normRH_GCaMP = (RH_GCaMP - 1)*100;
    filtRH_GCaMP = filtfilt(sos2,g2,normRH_GCaMP);
elseif any(strcmp(imagingWavelengths,{'Green','Lime'})) == true
    % HbT
    LH_HbT = ProcData.data.HbT.LH;
    filtLH_HbT = filtfilt(sos2,g2,LH_HbT);
    RH_HbT = ProcData.data.HbT.RH;
    filtRH_HbT = filtfilt(sos2,g2,RH_HbT);
end
% cortical and hippocampal spectrograms
specDataFile = [animalID '_' fileID '_SpecData.mat'];
load(specDataFile,'-mat');
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
cortical_RHnormS = SpecData.cortical_RH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
% Yvals for behavior Indices
whisking_Yvals = 1.10*max(filtLH_HbT)*ones(size(binWhiskers));
force_Yvals = 1.20*max(filtLH_HbT)*ones(size(binForce));
LPad_Yvals = 1.30*max(filtLH_HbT)*ones(size(LPadSol));
RPad_Yvals = 1.30*max(filtLH_HbT)*ones(size(RPadSol));
Aud_Yvals = 1.30*max(filtLH_HbT)*ones(size(AudSol));
Opto_Yvals = 1.30*max(filtLH_HbT)*ones(size(OptoLED));
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
    p5 = plot((1:length(filtLH_HbT))/ProcData.notes.wavelengthSamplingRate,filtLH_HbT,'color',colors('red'),'LineWidth',1);
    p6 = plot((1:length(filtRH_HbT))/ProcData.notes.wavelengthSamplingRate,filtRH_HbT,'color',colors('blue'),'LineWidth',1);
    legend([p5,p6,s1,s2,s3,s4,s5,s6],'LH HbT','RH HbT','movement','whisking',',LPad sol','RPad sol','Aud sol','Opto LED')
    ylabel('\DeltaHbT')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % GCaMP
    ax4 = subplot(8,1,4);
    p7 = plot((1:length(filtLH_GCaMP))/ProcData.notes.wavelengthSamplingRate,filtLH_GCaMP,'color',colors('red'),'LineWidth',1);
    hold on
    p8 = plot((1:length(filtRH_GCaMP))/ProcData.notes.wavelengthSamplingRate,filtRH_GCaMP,'color',colors('blue'),'LineWidth',1);
    legend([p7,p8],'LH GCaMP','RH GCaMP')
    ylabel('GCaMP /DeltaF/F')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % HbO,HbR
    ax5 = subplot(8,1,5);
    p9 = plot((1:length(filtLH_HbO))/ProcData.notes.wavelengthSamplingRate,filtLH_HbO,'color',colors('dark candy apple red'),'LineWidth',1);
    hold on
    p10 = plot((1:length(filtRH_HbO))/ProcData.notes.wavelengthSamplingRate,filtRH_HbO,'color',colors('dark candy apple red'),'LineWidth',1);
    p11 = plot((1:length(filtLH_HbR))/ProcData.notes.wavelengthSamplingRate,filtLH_HbR,'color',colors('dark candy apple red'),'LineWidth',1);
    p12 = plot((1:length(filtRH_HbR))/ProcData.notes.wavelengthSamplingRate,filtRH_HbR,'color',colors('dark candy apple red'),'LineWidth',1);
    legend([p9,p10,p11,p12],'LH HbO','RH HbO','LH HbR','RH HbR')
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
    clim([-100,100])
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
    clim([-100,100])
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
    clim([-100,100])
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
    ax3 = subplot(8,1,3:4);
    s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('sapphire'));
    hold on
    s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('blue-green'));
    s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
    s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
    s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
    s6 = scatter(OptoLED,Opto_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','b');
    p5 = plot((1:length(filtLH_HbT))/ProcData.notes.wavelengthSamplingRate,filtLH_HbT,'color',colors('red'),'LineWidth',1);
    p6 = plot((1:length(filtRH_HbT))/ProcData.notes.wavelengthSamplingRate,filtRH_HbT,'color',colors('blue'),'LineWidth',1);
    legend([p5,p6,s1,s2,s3,s4,s5,s6],'LH HbT','RH HbT','movement','whisking',',LPad sol','RPad sol','Aud sol','Opto LED')
    ylabel('\DeltaHbT')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % GCaMP
    ax4 = subplot(8,1,5);
    p7 = plot((1:length(filtLH_GCaMP))/ProcData.notes.wavelengthSamplingRate,filtLH_GCaMP,'color',colors('red'),'LineWidth',1);
    hold on
    p8 = plot((1:length(filtRH_GCaMP))/ProcData.notes.wavelengthSamplingRate,filtRH_GCaMP,'color',colors('blue'),'LineWidth',1);
    legend([p7,p8],'LH GCaMP','RH GCaMP')
    ylabel('GCaMP /DeltaF/F')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % left cortical electrode spectrogram
    ax5 = subplot(8,1,6);
    Semilog_ImageSC(T,F,cortical_LHnormS,'y')
    axis xy
    c4 = colorbar;
    ylabel(c4,'\DeltaP/P (%)')
    clim([-100,100])
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
    ax6 = subplot(8,1,7);
    Semilog_ImageSC(T,F,cortical_RHnormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)')
    clim([-100,100])
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
    ax7 = subplot(8,1,8);
    Semilog_ImageSC(T,F,hippocampusNormS,'y')
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)')
    clim([-100,100])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'box','off')
    yyaxis right
    ylabel('Hippocampal LFP')
    set(gca,'Yticklabel',[])
    % axes properties
    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'x')
    ax1Pos = get(ax1,'position');
    ax6Pos = get(ax6,'position');
    ax7Pos = get(ax7,'position');
    ax6Pos(3:4) = ax1Pos(3:4);
    ax7Pos(3:4) = ax1Pos(3:4);
    set(ax6,'position',ax6Pos);
    set(ax7,'position',ax7Pos);
elseif any(strcmp(imagingWavelengths,{'Green','Lime'})) == true
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
    ax3 = subplot(8,1,3:5);
    s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('sapphire'));
    hold on
    s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('blue-green'));
    s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
    s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
    s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
    s6 = scatter(OptoLED,Opto_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','b');
    p5 = plot((1:length(filtLH_HbT))/ProcData.notes.wavelengthSamplingRate,filtLH_HbT,'color',colors('red'),'LineWidth',1);
    p6 = plot((1:length(filtRH_HbT))/ProcData.notes.wavelengthSamplingRate,filtRH_HbT,'color',colors('blue'),'LineWidth',1);
    legend([p5,p6,s1,s2,s3,s4,s5,s6],'LH HbT','RH HbT','movement','whisking',',LPad sol','RPad sol','Aud sol','Opto LED')
    ylabel('\DeltaHbT')
    xlim([0,ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % left cortical electrode spectrogram
    ax4 = subplot(8,1,6);
    Semilog_ImageSC(T,F,cortical_LHnormS,'y')
    axis xy
    c4 = colorbar;
    ylabel(c4,'\DeltaP/P (%)')
    clim([-100,100])
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
    ax5 = subplot(8,1,7);
    Semilog_ImageSC(T,F,cortical_RHnormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)')
    clim([-100,100])
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
    ax6 = subplot(8,1,8);
    Semilog_ImageSC(T,F,hippocampusNormS,'y')
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)')
    clim([-100,100])
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