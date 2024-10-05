function [AnalysisResults] = Fig5_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1 for HRF manuscript
%________________________________________________________________________________________________________________________

%% arousal-state colors
colorRest = [(0/256),(166/256),(81/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAll = [(183/256),(115/256),(51/256)];
colorAlert = [(255/256),(191/256),(0/256)];
colorAsleep = [(0/256),(128/256),(255/256)];
%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
neuralBands = {'gammaBandPower','muaPower'};
behaviors = {'Contra','Whisk','Rest','NREM','REM','All','Alert','Asleep'};
samplingRate = 30;   % Hz
%% extract HRF kernels, R, and R2 for each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(neuralBands)
        neuralBand = neuralBands{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            try
                % coherence
                data.f.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.Coherence.(animalID).(neuralBand).cortLH.(behavior).FM_f;
                data.f.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.Coherence.(animalID).(neuralBand).cortRH.(behavior).FM_f;
                data.C.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.Coherence.(animalID).(neuralBand).cortLH.(behavior).FM_C;
                data.C.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.Coherence.(animalID).(neuralBand).cortRH.(behavior).FM_C;
                % power
                data.predf.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.Power.(animalID).(neuralBand).cortLH.(behavior).FM_pred_f;
                data.predf.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.Power.(animalID).(neuralBand).cortRH.(behavior).FM_pred_f;
                data.predS.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.Power.(animalID).(neuralBand).cortLH.(behavior).FM_pred_S;
                data.predS.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.Power.(animalID).(neuralBand).cortRH.(behavior).FM_pred_S;
                data.actf.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.Power.(animalID).(neuralBand).cortLH.(behavior).FM_act_f;
                data.actf.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.Power.(animalID).(neuralBand).cortRH.(behavior).FM_act_f;
                data.actS.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.Power.(animalID).(neuralBand).cortLH.(behavior).FM_act_S;
                data.actS.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.Power.(animalID).(neuralBand).cortRH.(behavior).FM_act_S;
            catch
                % coherence
                data.f.(neuralBand).(behavior).cortLH(aa,:) = NaN(1,size(data.f.(neuralBand).(behavior).cortLH,2));
                data.f.(neuralBand).(behavior).cortRH(aa,:) = NaN(1,size(data.f.(neuralBand).(behavior).cortRH,2));
                data.C.(neuralBand).(behavior).cortLH(aa,:) = NaN(1,size(data.C.(neuralBand).(behavior).cortLH,2));
                data.C.(neuralBand).(behavior).cortRH(aa,:) = NaN(1,size(data.C.(neuralBand).(behavior).cortRH,2));
                % power
                data.predf.(neuralBand).(behavior).cortLH(aa,:) = NaN(1,size(data.predf.(neuralBand).(behavior).cortLH,2));
                data.predf.(neuralBand).(behavior).cortRH(aa,:) = NaN(1,size(data.predf.(neuralBand).(behavior).cortRH,2));
                data.predS.(neuralBand).(behavior).cortLH(aa,:) = NaN(1,size(data.predS.(neuralBand).(behavior).cortLH,2));
                data.predS.(neuralBand).(behavior).cortRH(aa,:) = NaN(1,size(data.predS.(neuralBand).(behavior).cortRH,2));
                data.actf.(neuralBand).(behavior).cortLH(aa,:) = NaN(1,size(data.actf.(neuralBand).(behavior).cortLH,2));
                data.actf.(neuralBand).(behavior).cortRH(aa,:) = NaN(1,size(data.actf.(neuralBand).(behavior).cortRH,2));
                data.actS.(neuralBand).(behavior).cortLH(aa,:) = NaN(1,size(data.actS.(neuralBand).(behavior).cortLH,2));
                data.actS.(neuralBand).(behavior).cortRH(aa,:) = NaN(1,size(data.actS.(neuralBand).(behavior).cortRH,2));
            end
        end
    end
end
% concatenate the data from left and right into a single data set
for dd = 1:length(neuralBands)
    neuralBand = neuralBands{1,dd};
    for ee = 1:length(behaviors)
        behavior = behaviors{1,ee};
        % coherence
        data.f.(neuralBand).(behavior).cat = cat(1,data.f.(neuralBand).(behavior).cortLH,data.f.(neuralBand).(behavior).cortRH);
        data.C.(neuralBand).(behavior).cat = cat(1,data.C.(neuralBand).(behavior).cortLH,data.C.(neuralBand).(behavior).cortRH);
        % power
        data.predf.(neuralBand).(behavior).cat = cat(1,data.predf.(neuralBand).(behavior).cortLH,data.predf.(neuralBand).(behavior).cortRH);
        data.predS.(neuralBand).(behavior).cat = cat(1,data.predS.(neuralBand).(behavior).cortLH,data.predS.(neuralBand).(behavior).cortRH);
        data.actf.(neuralBand).(behavior).cat = cat(1,data.actf.(neuralBand).(behavior).cortLH,data.actf.(neuralBand).(behavior).cortRH);
        data.actS.(neuralBand).(behavior).cat = cat(1,data.actS.(neuralBand).(behavior).cortLH,data.actS.(neuralBand).(behavior).cortRH);
    end
end
% mean and std of each arousal-state
for ii = 1:length(neuralBands)
    neuralBand = neuralBands{1,ii};
    for jj = 1:length(behaviors)
        behavior = behaviors{1,jj};
        % coherence mean/std
        data.f.(neuralBand).(behavior).meanf = nanmean(data.f.(neuralBand).(behavior).cat,1);
        data.f.(neuralBand).(behavior).StDf = nanstd(data.f.(neuralBand).(behavior).cat,0,1);
        data.C.(neuralBand).(behavior).meanC = nanmean(data.C.(neuralBand).(behavior).cat,1);
        data.C.(neuralBand).(behavior).StDC = nanstd(data.C.(neuralBand).(behavior).cat,0,1);
        % power mean/std
        data.predf.(neuralBand).(behavior).meanf = nanmean(data.predf.(neuralBand).(behavior).cat,1);
        data.predf.(neuralBand).(behavior).StDf = nanstd(data.predf.(neuralBand).(behavior).cat,0,1);
        data.predS.(neuralBand).(behavior).meanS = nanmean(data.predS.(neuralBand).(behavior).cat,1);
        data.predS.(neuralBand).(behavior).StDS = nanstd(data.predS.(neuralBand).(behavior).cat,0,1);
        data.actf.(neuralBand).(behavior).meanf = nanmean(data.actf.(neuralBand).(behavior).cat,1);
        data.actf.(neuralBand).(behavior).StDf = nanstd(data.actf.(neuralBand).(behavior).cat,0,1);
        data.actS.(neuralBand).(behavior).meanS = nanmean(data.actS.(neuralBand).(behavior).cat,1);
        data.actS.(neuralBand).(behavior).StDS = nanstd(data.actS.(neuralBand).(behavior).cat,0,1);
    end
end
%% Fig.5
figure('Name','Fig5 (-)');
sgtitle('Figure 5')
%% [5A] actual HbT power
subplot(2,3,[1,4])
loglog(data.actf.gammaBandPower.Rest.meanf,data.actS.gammaBandPower.Rest.meanS,'color',colorRest,'LineWidth',2)
hold on
loglog(data.actf.gammaBandPower.Whisk.meanf,data.actS.gammaBandPower.Whisk.meanS,'color',colorWhisk,'LineWidth',2)
loglog(data.actf.gammaBandPower.Contra.meanf,data.actS.gammaBandPower.Contra.meanS,'color',colorStim,'LineWidth',2)
loglog(data.actf.gammaBandPower.NREM.meanf,data.actS.gammaBandPower.NREM.meanS,'color',colorNREM,'LineWidth',2)
loglog(data.actf.gammaBandPower.REM.meanf,data.actS.gammaBandPower.REM.meanS,'color',colorREM,'LineWidth',2)
loglog(data.actf.gammaBandPower.All.meanf,data.actS.gammaBandPower.All.meanS,'color',colorAll,'LineWidth',2)
loglog(data.actf.gammaBandPower.Alert.meanf,data.actS.gammaBandPower.Alert.meanS,'color',colorAlert,'LineWidth',2)
loglog(data.actf.gammaBandPower.Asleep.meanf,data.actS.gammaBandPower.Asleep.meanS,'color',colorAsleep,'LineWidth',2)
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
title('[5A] Actual \Delta[HbT] (\muM) power')
xlim([0.003,2])
% ylim([0,1])
set(gca,'box','off')
axis square
%% [5B] gamma-band kernel coherence
subplot(2,3,2)
loglog(data.predf.gammaBandPower.Rest.meanf,data.predS.gammaBandPower.Rest.meanS,'color',colorRest,'LineWidth',2)
hold on
loglog(data.predf.gammaBandPower.Whisk.meanf,data.predS.gammaBandPower.Whisk.meanS,'color',colorWhisk,'LineWidth',2)
loglog(data.predf.gammaBandPower.Contra.meanf,data.predS.gammaBandPower.Contra.meanS,'color',colorStim,'LineWidth',2)
loglog(data.predf.gammaBandPower.NREM.meanf,data.predS.gammaBandPower.NREM.meanS,'color',colorNREM,'LineWidth',2)
loglog(data.predf.gammaBandPower.REM.meanf,data.predS.gammaBandPower.REM.meanS,'color',colorREM,'LineWidth',2)
loglog(data.predf.gammaBandPower.All.meanf,data.predS.gammaBandPower.All.meanS,'color',colorAll,'LineWidth',2)
loglog(data.predf.gammaBandPower.Alert.meanf,data.predS.gammaBandPower.Alert.meanS,'color',colorAlert,'LineWidth',2)
loglog(data.predf.gammaBandPower.Asleep.meanf,data.predS.gammaBandPower.Asleep.meanS,'color',colorAsleep,'LineWidth',2)
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
title('[5B] Gamma-band Pred \Delta[HbT] (\muM) power')
xlim([0.003,2])
% ylim([0,1])
set(gca,'box','off')
axis square
%% [5C] gamma-band kernel coherence
subplot(2,3,3)
semilogx(data.f.gammaBandPower.Rest.meanf,data.C.gammaBandPower.Rest.meanC,'color',colorRest,'LineWidth',2)
hold on
semilogx(data.f.gammaBandPower.Whisk.meanf,data.C.gammaBandPower.Whisk.meanC,'color',colorWhisk,'LineWidth',2)
semilogx(data.f.gammaBandPower.Contra.meanf,data.C.gammaBandPower.Contra.meanC,'color',colorStim,'LineWidth',2)
semilogx(data.f.gammaBandPower.NREM.meanf,data.C.gammaBandPower.NREM.meanC,'color',colorNREM,'LineWidth',2)
semilogx(data.f.gammaBandPower.REM.meanf,data.C.gammaBandPower.REM.meanC,'color',colorREM,'LineWidth',2)
semilogx(data.f.gammaBandPower.All.meanf,data.C.gammaBandPower.All.meanC,'color',colorAll,'LineWidth',2)
semilogx(data.f.gammaBandPower.Alert.meanf,data.C.gammaBandPower.Alert.meanC,'color',colorAlert,'LineWidth',2)
semilogx(data.f.gammaBandPower.Asleep.meanf,data.C.gammaBandPower.Asleep.meanC,'color',colorAsleep,'LineWidth',2)
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence')
xlabel('Freq (Hz)')
title('[5C] Gamma-band Pred vs. Actual \Delta[HbT] (\muM) Coherence')
xlim([0.003,2])
ylim([0,1])
set(gca,'box','off')
axis square
%% [5A] gamma-band kernel coherence
subplot(2,3,5)
loglog(data.predf.muaPower.Rest.meanf,data.predS.muaPower.Rest.meanS,'color',colorRest,'LineWidth',2)
hold on
loglog(data.predf.muaPower.Whisk.meanf,data.predS.muaPower.Whisk.meanS,'color',colorWhisk,'LineWidth',2)
loglog(data.predf.muaPower.Contra.meanf,data.predS.muaPower.Contra.meanS,'color',colorStim,'LineWidth',2)
loglog(data.predf.muaPower.NREM.meanf,data.predS.muaPower.NREM.meanS,'color',colorNREM,'LineWidth',2)
loglog(data.predf.muaPower.REM.meanf,data.predS.muaPower.REM.meanS,'color',colorREM,'LineWidth',2)
loglog(data.predf.muaPower.All.meanf,data.predS.muaPower.All.meanS,'color',colorAll,'LineWidth',2)
loglog(data.predf.muaPower.Alert.meanf,data.predS.muaPower.Alert.meanS,'color',colorAlert,'LineWidth',2)
loglog(data.predf.muaPower.Asleep.meanf,data.predS.muaPower.Asleep.meanS,'color',colorAsleep,'LineWidth',2)
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
title('[5D] MUA Pred \Delta[HbT] (\muM) power')
xlim([0.003,2])
% ylim([0,1])
set(gca,'box','off')
axis square
%% [5E] gamma-band kernel coherence
subplot(2,3,6)
semilogx(data.f.muaPower.Rest.meanf,data.C.muaPower.Rest.meanC,'color',colorRest,'LineWidth',2)
hold on
semilogx(data.f.muaPower.Whisk.meanf,data.C.muaPower.Whisk.meanC,'color',colorWhisk,'LineWidth',2)
semilogx(data.f.muaPower.Contra.meanf,data.C.muaPower.Contra.meanC,'color',colorStim,'LineWidth',2)
semilogx(data.f.muaPower.NREM.meanf,data.C.muaPower.NREM.meanC,'color',colorNREM,'LineWidth',2)
semilogx(data.f.muaPower.REM.meanf,data.C.muaPower.REM.meanC,'color',colorREM,'LineWidth',2)
semilogx(data.f.muaPower.All.meanf,data.C.muaPower.All.meanC,'color',colorAll,'LineWidth',2)
semilogx(data.f.muaPower.Alert.meanf,data.C.muaPower.Alert.meanC,'color',colorAlert,'LineWidth',2)
semilogx(data.f.muaPower.Asleep.meanf,data.C.muaPower.Asleep.meanC,'color',colorAsleep,'LineWidth',2)
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence')
xlabel('Freq (Hz)')
title('[5E] MUA Pred vs. Actual \Delta[HbT] (\muM) Coherence')
xlim([0.003,2])
% ylim([0,1])
set(gca,'box','off')
axis square
%% save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
%     if ~exist(dirpath, 'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure,[dirpath 'Fig1']);
%     set(summaryFigure,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'Fig1'])
% end

end
