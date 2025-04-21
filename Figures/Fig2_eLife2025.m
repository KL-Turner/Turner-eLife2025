function [] = Fig2_eLife2025(rootFolder,saveState,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
Ephys_Fs = 30;
GCaMP_Fs = 10;

%% Awake ephys example
exampleProcFile_AwakeEphys = 'T165_210222_12_45_38_ProcData.mat';
load(exampleProcFile_AwakeEphys)
exampleSpecFile_AwakeEphys = 'T165_210222_12_45_38_SpecDataA.mat';
load(exampleSpecFile_AwakeEphys)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z_AwakeEphys,p_AwakeEphys,k_AwakeEphys] = butter(4,1/(Ephys_Fs/2),'low');
[sos_AwakeEphys,g_AwakeEphys] = zp2sos(z_AwakeEphys,p_AwakeEphys,k_AwakeEphys);

[z,p,k] = butter(1,1/(Ephys_Fs/2),'low');
[sos,g] = zp2sos(z,p,k);

binWhiskers_AwakeEphys = ProcData.data.binWhiskerAngle;
% data
HbT_AwakeEphys = filtfilt(sos_AwakeEphys,g_AwakeEphys,ProcData.data.HbT.RH);
% RH gamma baseline for T165 Jan 30
gammaBaseline_AwakeEphys = 1.016060094087967e-10;
gamma_AwakeEphys = filtfilt(sos,g,((ProcData.data.cortical_RH.gammaBandPower - gammaBaseline_AwakeEphys)./gammaBaseline_AwakeEphys)*100);
% cortical and hippocampal spectrograms
cortNormS_AwakeEphys = SpecData.cortical_RH.normS.*100;
T_AwakeEphys = SpecData.cortical_RH.T;
F_AwakeEphys = SpecData.cortical_RH.F;
% Yvals for behavior Indices
whisking_Yvals_AwakeEphys = 155*ones(size(binWhiskers_AwakeEphys));
whiskInds_AwakeEphys = binWhiskers_AwakeEphys.*whisking_Yvals_AwakeEphys;
% set whisk indeces
for x = 1:length(whiskInds_AwakeEphys)
    if whiskInds_AwakeEphys(1,x) == 0
        whiskInds_AwakeEphys(1,x) = NaN;
    end
end

%% Sleep ephys example
exampleProcFile_AsleepEphys = 'T165_210223_12_37_03_ProcData.mat';
load(exampleProcFile_AsleepEphys)
exampleSpecFile_AsleepEphys = 'T165_210223_12_37_03_SpecDataA.mat';
load(exampleSpecFile_AsleepEphys)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z_AsleepEphys,p_AsleepEphys,k_AsleepEphys] = butter(4,1/(Ephys_Fs/2),'low');
[sos_AsleepEphys,g_AsleepEphys] = zp2sos(z_AsleepEphys,p_AsleepEphys,k_AsleepEphys);
binWhiskers_AsleepEphys = ProcData.data.binWhiskerAngle;
% data
HbT_AsleepEphys = filtfilt(sos_AsleepEphys,g_AsleepEphys,ProcData.data.HbT.RH);
% RH gamma baseline for T151 Jan 18
gammaBaseline_AsleepEphys = 3.513854434241187e-10;
gamma_AsleepEphys = filtfilt(sos,g,((ProcData.data.cortical_RH.gammaBandPower - gammaBaseline_AsleepEphys)./gammaBaseline_AsleepEphys)*100);
% cortical and hippocampal spectrograms
hipNormS_AsleepEphys = SpecData.hippocampus.normS.*100;
T_AsleepEphys = SpecData.cortical_LH.T;
F_AsleepEphys = SpecData.cortical_LH.F;
% Yvals for behavior Indices
whisking_Yvals_AsleepEphys = 155*ones(size(binWhiskers_AsleepEphys));
whiskInds_AsleepEphys = binWhiskers_AsleepEphys.*whisking_Yvals_AsleepEphys;
% set whisk indeces
for x = 1:length(whiskInds_AsleepEphys)
    if whiskInds_AsleepEphys(1,x) == 0
        whiskInds_AsleepEphys(1,x) = NaN;
    end
end

%% Awake GCaMP example
exampleProcFile_AwakeGCaMP = 'T233_220206_13_07_08_ProcData.mat';
load(exampleProcFile_AwakeGCaMP)
exampleSpecFile_AwakeGCaMP = 'T233_220206_13_07_08_SpecDataA.mat';
load(exampleSpecFile_AwakeGCaMP)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z_AwakeGCaMP,p_AwakeGCaMP,k_AwakeGCaMP] = butter(4,1/(GCaMP_Fs/2),'low');
[sos_AwakeGCaMP,g_AwakeGCaMP] = zp2sos(z_AwakeGCaMP,p_AwakeGCaMP,k_AwakeGCaMP);
binWhiskers_AwakeGCaMP = ProcData.data.binWhiskerAngle;
% data
HbT_AwakeGCaMP = filtfilt(sos_AwakeGCaMP,g_AwakeGCaMP,ProcData.data.HbT.RH);
HbO_AwakeGCaMP = filtfilt(sos_AwakeGCaMP,g_AwakeGCaMP,ProcData.data.HbO.RH);
HbR_AwakeGCaMP = filtfilt(sos_AwakeGCaMP,g_AwakeGCaMP,ProcData.data.HbR.RH);
GCaMP_AwakeGCaMP = filtfilt(sos_AwakeGCaMP,g_AwakeGCaMP,ProcData.data.GCaMP.RH);
% cortical and hippocampal spectrograms
cortNormS_AwakeGCaMP = SpecData.cortical_RH.normS.*100;
T_AwakeGCaMP = SpecData.cortical_LH.T;
F_AwakeGCaMP = SpecData.cortical_LH.F;
% Yvals for behavior Indices
whisking_Yvals_AwakeGCaMP = 155*ones(size(binWhiskers_AwakeGCaMP));
whiskInds_AwakeGCaMP = binWhiskers_AwakeGCaMP.*whisking_Yvals_AwakeGCaMP;
% set whisk indeces
for x = 1:length(whiskInds_AwakeGCaMP)
    if whiskInds_AwakeGCaMP(1,x) == 0
        whiskInds_AwakeGCaMP(1,x) = NaN;
    end
end

%% Sleep GCaMP example
exampleProcFile_AsleepGCaMP = 'T233_220206_14_10_44_ProcData.mat';
load(exampleProcFile_AsleepGCaMP)
exampleSpecFile_AsleepGCaMP = 'T233_220206_14_10_44_SpecDataA.mat';
load(exampleSpecFile_AsleepGCaMP)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z_AsleepGCaMP,p_AsleepGCaMP,k_AsleepGCaMP] = butter(4,1/(GCaMP_Fs/2),'low');
[sos_AsleepGCaMP,g_AsleepGCaMP] = zp2sos(z_AsleepGCaMP,p_AsleepGCaMP,k_AsleepGCaMP);
binWhiskers_AsleepGCaMP = ProcData.data.binWhiskerAngle;
% data
HbT_AsleepGCaMP = filtfilt(sos_AsleepGCaMP,g_AsleepGCaMP,ProcData.data.HbT.RH);
HbO_AsleepGCaMP = filtfilt(sos_AsleepGCaMP,g_AsleepGCaMP,ProcData.data.HbO.RH);
HbR_AsleepGCaMP = filtfilt(sos_AsleepGCaMP,g_AsleepGCaMP,ProcData.data.HbR.RH);
GCaMP_AsleepGCaMP = filtfilt(sos_AsleepGCaMP,g_AsleepGCaMP,ProcData.data.GCaMP.RH);
% cortical and hippocampal spectrograms
hipNormS_AsleepGCaMP = SpecData.hippocampus.normS.*100;
T_AsleepGCaMP = SpecData.cortical_LH.T;
F_AsleepGCaMP = SpecData.cortical_LH.F;
% Yvals for behavior Indices
whisking_Yvals_AsleepGCaMP = 155*ones(size(binWhiskers_AsleepGCaMP));
whiskInds_AsleepGCaMP = binWhiskers_AsleepGCaMP.*whisking_Yvals_AsleepGCaMP;
% set whisk indeces
for x = 1:length(whiskInds_AsleepGCaMP)
    if whiskInds_AsleepGCaMP(1,x) == 0
        whiskInds_AsleepGCaMP(1,x) = NaN;
    end
end

%% Figure 2
Fig2 = figure('Name','Fig. 2');

% Hemodynamic and behavioral indeces
ax1 = subplot(4,3,1);
s1 = scatter((1:length(binWhiskers_AwakeEphys))/ProcData.notes.dsFs,whiskInds_AwakeEphys,'.','MarkerEdgeColor',colors('black'));
hold on;
p1 = plot((1:length(HbT_AwakeEphys))/Ephys_Fs,HbT_AwakeEphys,'color',colors('electric purple'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-50,160])
yyaxis right
p2 = plot((1:length(gamma_AwakeEphys))/Ephys_Fs,gamma_AwakeEphys,'color',colors('sapphire'),'LineWidth',1);
ylim([-300,1300])
legend([p1,p2,s1],'HbT','gamma power','whisking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([310,370,430,490])
title('RH SIBF')
xlim([310,490])
ax1.YAxis(1).Color = colors('black');
ax1.YAxis(2).Color = colors('magenta');

% Cortical electrode spectrogram
ax2 = subplot(4,3,4);
Semilog_ImageSC(T_AwakeEphys,F_AwakeEphys,cortNormS_AwakeEphys,'y')
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
title('Cortical LFP')
xlabel('Time (min)')
ylabel('Frequency (Hz)')
xlim([310,490])
set(gca,'box','off')
xticks([310,370,430,490])
xticklabels({'0','1','2','3'});
% axis properties
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax2Pos(3) = ax1Pos(3);
set(ax2,'position',ax2Pos);

% Hemodynamic and behavioral indeces
ax3 = subplot(4,3,[2,3]);
s1 = scatter((1:length(binWhiskers_AsleepEphys))/ProcData.notes.dsFs,whiskInds_AsleepEphys,'.','MarkerEdgeColor',colors('black'));
hold on;
p1 = plot((1:length(HbT_AsleepEphys))/Ephys_Fs,HbT_AsleepEphys,'color',colors('electric purple'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-50,160])
yyaxis right
p2 = plot((1:length(gamma_AsleepEphys))/Ephys_Fs,gamma_AsleepEphys,'color',colors('sapphire'),'LineWidth',1);
ylim([-300,1700])
legend([p1,p2,s1],'HbT','gamma power','whisking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([450,510,570,630,690,750,810])
title('RH SIBF')
xlim([450,810])
ax3.YAxis(1).Color = colors('black');
ax3.YAxis(2).Color = colors('magenta');

% Hippocampal electrode spectrogram
ax4 = subplot(4,3,[5,6]);
Semilog_ImageSC(T_AsleepEphys,F_AsleepEphys,hipNormS_AsleepEphys,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
title('Hippocampal LFP')
xlabel('Time (min)')
ylabel('Frequency (Hz)')
xlim([450,810])
set(gca,'box','off')
xticks([450,510,570,630,690,750,810])
xticklabels({'0','1','2','3','4','5','6'})
% axis properties
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax4Pos(3) = ax3Pos(3);
set(ax4,'position',ax4Pos);

% Hemodynamic and behavioral indeces
ax5 = subplot(4,3,7);
s1 = scatter((1:length(binWhiskers_AwakeGCaMP))/ProcData.notes.dsFs,whiskInds_AwakeGCaMP,'.','MarkerEdgeColor',colors('black'));
hold on;
p1 = plot((1:length(HbO_AwakeGCaMP))/GCaMP_Fs,HbO_AwakeGCaMP,'color',colors('deep carrot orange'),'LineWidth',1);
p2 = plot((1:length(HbR_AwakeGCaMP))/GCaMP_Fs,HbR_AwakeGCaMP,'color',colors('vegas gold'),'LineWidth',1);
p3 = plot((1:length(HbT_AwakeGCaMP))/GCaMP_Fs,HbT_AwakeGCaMP,'color',colors('electric purple'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-50,160])
yyaxis right
p4 = plot((1:length(GCaMP_AwakeGCaMP))/GCaMP_Fs,(GCaMP_AwakeGCaMP - 1)*100,'color',colors('north texas green'),'LineWidth',1);
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
ylim([-7,20])
legend([p1,p2,p3,p4,s1],'HbO','HbR','HbT','GCaMP','whisking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([260,320,380,440])
title('RH SIBF')
xlim([260,440])
ax5.YAxis(1).Color = colors('black');
ax5.YAxis(2).Color = colors('north texas green');

% Cortical electrode spectrogram
ax6 = subplot(4,3,10);
Semilog_ImageSC(T_AwakeGCaMP,F_AwakeGCaMP,cortNormS_AwakeGCaMP,'y')
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
title('Cortical LFP')
xlabel('Time (min)')
ylabel('Frequency (Hz)')
xlim([260,440])
set(gca,'box','off')
xticks([260,320,380,440])
xticklabels({'0','1','2','3'})
% axis properties
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax6Pos(3) = ax5Pos(3);
set(ax6,'position',ax6Pos);

% Hemodynamic and behavioral indeces
ax7 = subplot(4,3,[8,9]);
s1 = scatter((1:length(binWhiskers_AsleepGCaMP))/ProcData.notes.dsFs,whiskInds_AsleepGCaMP,'.','MarkerEdgeColor',colors('black'));
hold on;
p1 = plot((1:length(HbO_AsleepGCaMP))/GCaMP_Fs,HbO_AsleepGCaMP,'color',colors('deep carrot orange'),'LineWidth',1);
p2 = plot((1:length(HbR_AsleepGCaMP))/GCaMP_Fs,HbR_AsleepGCaMP,'color',colors('vegas gold'),'LineWidth',1);
p3 = plot((1:length(HbT_AsleepGCaMP))/GCaMP_Fs,HbT_AsleepGCaMP,'color',colors('electric purple'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-50,160])
yyaxis right
p4 = plot((1:length(GCaMP_AsleepGCaMP))/GCaMP_Fs,(GCaMP_AsleepGCaMP - 1)*100,'color',colors('north texas green'),'LineWidth',1);
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
ylim([-7,20])
legend([p1,p2,p3,p4,s1],'HbO','HbR','HbT','GCaMP','whisking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([215,275,335,395,455,515,575])
title('RH SIBF')
xlim([215,575])
ax7.YAxis(1).Color = colors('black');
ax7.YAxis(2).Color = colors('north texas green');

% Hippocampal electrode spectrogram
ax8 = subplot(4,3,[11,12]);
Semilog_ImageSC(T_AsleepGCaMP,F_AsleepGCaMP,hipNormS_AsleepGCaMP,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
title('Hippocampal LFP')
xlabel('Time (min)')
ylabel('Frequency (Hz)')
xlim([215,575])
set(gca,'box','off')
xticks([215,275,335,395,455,515,575])
xticklabels({'0','1','2','3','4','5','6'})
% axis properties
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');
ax8Pos(3) = ax7Pos(3);
set(ax8,'position',ax8Pos);

%% Save figure and stats
if saveState == true
    dirpath = [rootFolder delim 'MATLAB Figs/Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig2,[dirpath 'Fig2']);
    diaryFile = [dirpath 'Fig2_Stats.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    disp('======================================================================================================================')
    diary off
end
