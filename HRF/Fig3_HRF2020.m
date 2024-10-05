function [AnalysisResults] = Fig3_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
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
% excel sheet kernel success
msExcelFile = 'HRF_Kernel_Success.xlsx';
[~,~,alldata] = xlsread(msExcelFile);
kernelSheet.animalID_behavior = alldata(:,1);
kernelSheet.IR_gammaBandPower_cortLH = alldata(:,2);
kernelSheet.IR_gammaBandPower_cortRH = alldata(:,3);
kernelSheet.IR_gammaBandPower_hippLH = alldata(:,4);
kernelSheet.IR_muaPower_cortLH = alldata(:,5);
kernelSheet.IR_muaPower_cortRH = alldata(:,6);
kernelSheet.IR_muaPower_hippLH = alldata(:,7);
kernelSheet.FM_gammaBandPower_cortLH = alldata(:,8);
kernelSheet.FM_gammaBandPower_cortRH = alldata(:,9);
kernelSheet.FM_gammaBandPower_hippLH = alldata(:,10);
kernelSheet.FM_muaPower_cortLH = alldata(:,11);
kernelSheet.FM_muaPower_cortRH = alldata(:,12);
kernelSheet.FM_muaPower_hippLH = alldata(:,13);
%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
neuralBands = {'gammaBandPower','muaPower'};
behaviors = {'Contra','Whisk','Rest','NREM','REM','All','Alert','Asleep'};
dataTypes = {'IR_gammaBandPower_cortLH','IR_gammaBandPower_cortRH','FM_gammaBandPower_cortLH','FM_gammaBandPower_cortRH',...
    'IR_muaPower_cortLH','IR_muaPower_cortRH','FM_muaPower_cortLH','FM_muaPower_cortRH'};
%% extract HRF kernels, R, and R2 for each arousal-state
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        % corr coeff and IR_R^2
        animalID_behavior = [animalID '_' behavior];
        for cc = 1:length(kernelSheet.animalID_behavior)
            if strcmp(animalID_behavior,kernelSheet.animalID_behavior{cc,1}) == true
                for dd = 1:length(dataTypes)
                    dataTypeLogical(1,dd) = cell2mat(kernelSheet.(dataTypes{1,dd})(cc,1)); %#ok<AGROW>
                end
            end
        end
        % pull data from valid kernels
        for ee = 1:length(dataTypes)
            dataType = dataTypes{1,ee};
            strBreaks = strfind(dataType(1,:),'_');
            kernelType = dataType(1:strBreaks(1) - 1);
            neuralBand = dataType(strBreaks(1) + 1:strBreaks(2) - 1);
            hemisphere = dataType(strBreaks(2) + 1:end);
            if dataTypeLogical(1,ee) == true
                data.Predictions.([kernelType '_R']).(neuralBand).(behavior).(hemisphere){aa,1} = AnalysisResults.Predictions.(animalID).(neuralBand).(hemisphere).(behavior).([kernelType '_R']);
                data.Predictions.([kernelType '_R2']).(neuralBand).(behavior).(hemisphere){aa,1} = AnalysisResults.Predictions.(animalID).(neuralBand).(hemisphere).(behavior).([kernelType '_R2']);
            else
                data.Predictions.([kernelType '_R']).(neuralBand).(behavior).(hemisphere){aa,1} = NaN;
                data.Predictions.([kernelType '_R2']).(neuralBand).(behavior).(hemisphere){aa,1} = NaN;
            end
        end
    end
end
% concatenate the data.Predictions from left and right into a single data.Predictions set
for dd = 1:length(neuralBands)
    neuralBand = neuralBands{1,dd};
    for ee = 1:length(behaviors)
        behavior = behaviors{1,ee};
        % corr coef and IR_R^2
        data.Predictions.IR_R.(neuralBand).(behavior).cat = cat(1,data.Predictions.IR_R.(neuralBand).(behavior).cortLH,data.Predictions.IR_R.(neuralBand).(behavior).cortRH);
        data.Predictions.IR_R2.(neuralBand).(behavior).cat = cat(1,data.Predictions.IR_R2.(neuralBand).(behavior).cortLH,data.Predictions.IR_R2.(neuralBand).(behavior).cortRH);
        data.Predictions.FM_R.(neuralBand).(behavior).cat = cat(1,data.Predictions.FM_R.(neuralBand).(behavior).cortLH,data.Predictions.FM_R.(neuralBand).(behavior).cortRH);
        data.Predictions.FM_R2.(neuralBand).(behavior).cat = cat(1,data.Predictions.FM_R2.(neuralBand).(behavior).cortLH,data.Predictions.FM_R2.(neuralBand).(behavior).cortRH);
    end
end
% mean of individual hemispheres and cat ind points for histograms
for ff = 1:length(neuralBands)
    neuralBand = neuralBands{1,ff};
    for gg = 1:length(behaviors)
        behavior = behaviors{1,gg};
        % corr coef and IR_R^2
        for hh = 1:length(data.Predictions.IR_R.(neuralBand).(behavior).cat)
            data.Predictions.IR_R.(neuralBand).(behavior).indMeds(hh,1) = nanmedian(data.Predictions.IR_R.(neuralBand).(behavior).cat{hh,1});
            data.Predictions.IR_R2.(neuralBand).(behavior).indMeds(hh,1) = nanmedian(data.Predictions.IR_R2.(neuralBand).(behavior).cat{hh,1});
            data.Predictions.FM_R.(neuralBand).(behavior).indMeds(hh,1) = nanmedian(data.Predictions.FM_R.(neuralBand).(behavior).cat{hh,1});
            data.Predictions.FM_R2.(neuralBand).(behavior).indMeds(hh,1) = nanmedian(data.Predictions.FM_R2.(neuralBand).(behavior).cat{hh,1});
        end
    end
end
% mean and std of each arousal-state
for ii = 1:length(neuralBands)
    neuralBand = neuralBands{1,ii};
    for jj = 1:length(behaviors)
        behavior = behaviors{1,jj};
        % IR_R med/std
        data.Predictions.IR_R.(neuralBand).(behavior).median = nanmean(data.Predictions.IR_R.(neuralBand).(behavior).indMeds,1);
        data.Predictions.IR_R.(neuralBand).(behavior).medStD = nanstd(data.Predictions.IR_R.(neuralBand).(behavior).indMeds,0,1);
        % IR_R^2 med/std
        data.Predictions.IR_R2.(neuralBand).(behavior).median = nanmean(data.Predictions.IR_R2.(neuralBand).(behavior).indMeds,1);
        data.Predictions.IR_R2.(neuralBand).(behavior).medStD = nanstd(data.Predictions.IR_R2.(neuralBand).(behavior).indMeds,0,1);
        % FM_R med/std
        data.Predictions.FM_R.(neuralBand).(behavior).median = nanmean(data.Predictions.FM_R.(neuralBand).(behavior).indMeds,1);
        data.Predictions.FM_R.(neuralBand).(behavior).medStD = nanstd(data.Predictions.FM_R.(neuralBand).(behavior).indMeds,0,1);
        % FM_R^2 med/std
        data.Predictions.FM_R2.(neuralBand).(behavior).median = nanmean(data.Predictions.FM_R2.(neuralBand).(behavior).indMeds,1);
        data.Predictions.FM_R2.(neuralBand).(behavior).medStD = nanstd(data.Predictions.FM_R2.(neuralBand).(behavior).indMeds,0,1);
    end
end
%% Fig.3
summaryFigure = figure('Name','Fig3 (A-D)');
sgtitle('Figure 3')
%% [3A] gamma-band gamma HRF predictions - dedicated kernels
subplot(2,3,1);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.Predictions.IR_R.gammaBandPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Predictions.IR_R.gammaBandPower.Contra.median,data.Predictions.IR_R.gammaBandPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.Predictions.FM_R.gammaBandPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Predictions.FM_R.gammaBandPower.Contra.median,data.Predictions.FM_R.gammaBandPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Predictions.IR_R.gammaBandPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Predictions.IR_R.gammaBandPower.Whisk.median,data.Predictions.IR_R.gammaBandPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.Predictions.FM_R.gammaBandPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.Predictions.FM_R.gammaBandPower.Whisk.median,data.Predictions.FM_R.gammaBandPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Predictions.IR_R.gammaBandPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.Predictions.IR_R.gammaBandPower.Rest.median,data.Predictions.IR_R.gammaBandPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.Predictions.FM_R.gammaBandPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.Predictions.FM_R.gammaBandPower.Rest.median,data.Predictions.FM_R.gammaBandPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.Predictions.IR_R.gammaBandPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.Predictions.IR_R.gammaBandPower.NREM.median,data.Predictions.IR_R.gammaBandPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.Predictions.FM_R.gammaBandPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.Predictions.FM_R.gammaBandPower.NREM.median,data.Predictions.FM_R.gammaBandPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.Predictions.IR_R.gammaBandPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.Predictions.IR_R.gammaBandPower.REM.median,data.Predictions.IR_R.gammaBandPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.Predictions.FM_R.gammaBandPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.Predictions.FM_R.gammaBandPower.REM.median,data.Predictions.FM_R.gammaBandPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.Predictions.IR_R.gammaBandPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.Predictions.IR_R.gammaBandPower.All.median,data.Predictions.IR_R.gammaBandPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.Predictions.FM_R.gammaBandPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.Predictions.FM_R.gammaBandPower.All.median,data.Predictions.FM_R.gammaBandPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title('[3A] Gamma-band [30-100 Hz] IR_R')
ylabel('Corr. coef. (IR_R)')
xlim([0,13])
ylim([-0.3,1])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% [3B] gamma-band gamma HRF predictions
subplot(2,3,2);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.Predictions.IR_R2.gammaBandPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Predictions.IR_R2.gammaBandPower.Contra.median,data.Predictions.IR_R2.gammaBandPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.Predictions.FM_R2.gammaBandPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Predictions.FM_R2.gammaBandPower.Contra.median,data.Predictions.FM_R2.gammaBandPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Predictions.IR_R2.gammaBandPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Predictions.IR_R2.gammaBandPower.Whisk.median,data.Predictions.IR_R2.gammaBandPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.Predictions.FM_R2.gammaBandPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.Predictions.FM_R2.gammaBandPower.Whisk.median,data.Predictions.FM_R2.gammaBandPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Predictions.IR_R2.gammaBandPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.Predictions.IR_R2.gammaBandPower.Rest.median,data.Predictions.IR_R2.gammaBandPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.Predictions.FM_R2.gammaBandPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.Predictions.FM_R2.gammaBandPower.Rest.median,data.Predictions.FM_R2.gammaBandPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.Predictions.IR_R2.gammaBandPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.Predictions.IR_R2.gammaBandPower.NREM.median,data.Predictions.IR_R2.gammaBandPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.Predictions.FM_R2.gammaBandPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.Predictions.FM_R2.gammaBandPower.NREM.median,data.Predictions.FM_R2.gammaBandPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.Predictions.IR_R2.gammaBandPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.Predictions.IR_R2.gammaBandPower.REM.median,data.Predictions.IR_R2.gammaBandPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.Predictions.FM_R2.gammaBandPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.Predictions.FM_R2.gammaBandPower.REM.median,data.Predictions.FM_R2.gammaBandPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.Predictions.IR_R2.gammaBandPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.Predictions.IR_R2.gammaBandPower.All.median,data.Predictions.IR_R2.gammaBandPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.Predictions.FM_R2.gammaBandPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.Predictions.FM_R2.gammaBandPower.All.median,data.Predictions.FM_R2.gammaBandPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title('[3B] Gamma-band [30-100 Hz] IR_R^2')
ylabel('Coeff of det. (IR_R^2)')
xlim([0,13])
ylim([-0.3,1])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% [3C] gamma-band gamma HRF predictions
subplot(2,3,3);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.Predictions.IR_R2.gammaBandPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Predictions.IR_R2.gammaBandPower.All.median,data.Predictions.IR_R2.gammaBandPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.Predictions.FM_R2.gammaBandPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Predictions.FM_R2.gammaBandPower.All.median,data.Predictions.FM_R2.gammaBandPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Predictions.IR_R2.gammaBandPower.Alert.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Predictions.IR_R2.gammaBandPower.Alert.median,data.Predictions.IR_R2.gammaBandPower.Alert.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.Predictions.FM_R2.gammaBandPower.Alert.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAlert)*0.5 + colorAlert,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.Predictions.FM_R2.gammaBandPower.Alert.median,data.Predictions.FM_R2.gammaBandPower.Alert.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Predictions.IR_R2.gammaBandPower.Asleep.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.Predictions.IR_R2.gammaBandPower.Asleep.median,data.Predictions.IR_R2.gammaBandPower.Asleep.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.Predictions.FM_R2.gammaBandPower.Asleep.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAsleep)*0.5 + colorAsleep,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.Predictions.FM_R2.gammaBandPower.Asleep.median,data.Predictions.FM_R2.gammaBandPower.Asleep.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title('[3C] Gamma-band [30-100 Hz] IR_R^2')
ylabel('Coeff of det. (IR_R^2)')
xlim([0,7])
ylim([-0.3,1])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
axis square
%% [3D] MUA impulse HRF predictions
subplot(2,3,4);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.Predictions.IR_R.muaPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Predictions.IR_R.muaPower.Contra.median,data.Predictions.IR_R.muaPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.Predictions.FM_R.muaPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Predictions.FM_R.muaPower.Contra.median,data.Predictions.FM_R.muaPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Predictions.IR_R.muaPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Predictions.IR_R.muaPower.Whisk.median,data.Predictions.IR_R.muaPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.Predictions.FM_R.muaPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.Predictions.FM_R.muaPower.Whisk.median,data.Predictions.FM_R.muaPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Predictions.IR_R.muaPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.Predictions.IR_R.muaPower.Rest.median,data.Predictions.IR_R.muaPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.Predictions.FM_R.muaPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.Predictions.FM_R.muaPower.Rest.median,data.Predictions.FM_R.muaPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.Predictions.IR_R.muaPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.Predictions.IR_R.muaPower.NREM.median,data.Predictions.IR_R.muaPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.Predictions.FM_R.muaPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.Predictions.FM_R.muaPower.NREM.median,data.Predictions.FM_R.muaPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.Predictions.IR_R.muaPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.Predictions.IR_R.muaPower.REM.median,data.Predictions.IR_R.muaPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.Predictions.FM_R.muaPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.Predictions.FM_R.muaPower.REM.median,data.Predictions.FM_R.muaPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.Predictions.IR_R.muaPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.Predictions.IR_R.muaPower.All.median,data.Predictions.IR_R.muaPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.Predictions.FM_R.muaPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.Predictions.FM_R.muaPower.All.median,data.Predictions.FM_R.muaPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title('[3D] MUA [300-3000 Hz] IR_R')
ylabel('Corr. coef. (IR_R)')
xlim([0,13])
ylim([-0.3,1])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% [3E] MUA HRF predictions (IR_R2)
subplot(2,3,5);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.Predictions.IR_R2.muaPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Predictions.IR_R2.muaPower.Contra.median,data.Predictions.IR_R2.muaPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.Predictions.FM_R2.muaPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Predictions.FM_R2.muaPower.Contra.median,data.Predictions.FM_R2.muaPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Predictions.IR_R2.muaPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Predictions.IR_R2.muaPower.Whisk.median,data.Predictions.IR_R2.muaPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.Predictions.FM_R2.muaPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.Predictions.FM_R2.muaPower.Whisk.median,data.Predictions.FM_R2.muaPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Predictions.IR_R2.muaPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.Predictions.IR_R2.muaPower.Rest.median,data.Predictions.IR_R2.muaPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.Predictions.FM_R2.muaPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.Predictions.FM_R2.muaPower.Rest.median,data.Predictions.FM_R2.muaPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.Predictions.IR_R2.muaPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.Predictions.IR_R2.muaPower.NREM.median,data.Predictions.IR_R2.muaPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.Predictions.FM_R2.muaPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.Predictions.FM_R2.muaPower.NREM.median,data.Predictions.FM_R2.muaPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.Predictions.IR_R2.muaPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.Predictions.IR_R2.muaPower.REM.median,data.Predictions.IR_R2.muaPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.Predictions.FM_R2.muaPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.Predictions.FM_R2.muaPower.REM.median,data.Predictions.FM_R2.muaPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.Predictions.IR_R2.muaPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.Predictions.IR_R2.muaPower.All.median,data.Predictions.IR_R2.muaPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.Predictions.FM_R2.muaPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.Predictions.FM_R2.muaPower.All.median,data.Predictions.FM_R2.muaPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title('[3E] MUA [300-3000 Hz] IR_R^2')
ylabel('Coeff of det. (IR_R^2)')
xlim([0,13])
ylim([-0.3,1])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% [3F] MUA gamma HRF predictions
subplot(2,3,6);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.Predictions.IR_R2.muaPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Predictions.IR_R2.muaPower.All.median,data.Predictions.IR_R2.muaPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.Predictions.FM_R2.muaPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Predictions.FM_R2.muaPower.All.median,data.Predictions.FM_R2.muaPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.Predictions.IR_R2.muaPower.Alert.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Predictions.IR_R2.muaPower.Alert.median,data.Predictions.IR_R2.muaPower.Alert.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.Predictions.FM_R2.muaPower.Alert.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAlert)*0.5 + colorAlert,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.Predictions.FM_R2.muaPower.Alert.median,data.Predictions.FM_R2.muaPower.Alert.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.Predictions.IR_R2.muaPower.Asleep.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.Predictions.IR_R2.muaPower.Asleep.median,data.Predictions.IR_R2.muaPower.Asleep.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.Predictions.FM_R2.muaPower.Asleep.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAsleep)*0.5 + colorAsleep,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.Predictions.FM_R2.muaPower.Asleep.median,data.Predictions.FM_R2.muaPower.Asleep.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title('[3F] MUA [300-3000 Hz] IR_R^2')
ylabel('Coeff of det. (IR_R^2)')
xlim([0,7])
ylim([-0.3,1])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
axis square
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig3']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig3'])
end

end
