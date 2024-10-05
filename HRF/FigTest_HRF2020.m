function [AnalysisResults] = Fig1_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
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
%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
hemispheres = {'LH','LH','LH','LH','LH','RH','RH','LH','LH','RH','RH','RH','LH','LH'};
neuralBands = {'gammaBandPower','muaPower'};
behaviors = {'Contra','Whisk','Rest','NREM','REM','All'};
samplingRate = 30;   % Hz
%% extract HRF kernels, R, and R2 for each arousal-state
data.GammaHbT.gamma = [];
data.GammaHbT.HbT = [];
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(neuralBands)
        neuralBand = neuralBands{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            % pull HRFs from AnalysisResults.mat structure - dc shift each IR function by offset
            data.IR.(neuralBand).(behavior).adjLH(aa,:) = AnalysisResults.HRFs.(animalID).HbT.(neuralBand).adjLH.(behavior).IR_function_B - AnalysisResults.HRFs.(animalID).HbT.(neuralBand).adjLH.(behavior).IR_function_B(1);
            data.IR.(neuralBand).(behavior).adjRH(aa,:) = AnalysisResults.HRFs.(animalID).HbT.(neuralBand).adjRH.(behavior).IR_function_B - AnalysisResults.HRFs.(animalID).HbT.(neuralBand).adjRH.(behavior).IR_function_B(1);
            data.IR_gamma.(neuralBand).(behavior).adjLH(aa,:) = AnalysisResults.HRFs.(animalID).HbT.(neuralBand).adjLH.(behavior).IR_gammaFunction;
            data.IR_gamma.(neuralBand).(behavior).adjRH(aa,:) = AnalysisResults.HRFs.(animalID).HbT.(neuralBand).adjRH.(behavior).IR_gammaFunction;
            % pull HRF predictions (R and R2) from AnalysisResults.mat structure
            if strcmp(behavior,'All') == true
                % corr coef and R^2
                data.R.(neuralBand).(behavior).adjLH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjLH.All.(behavior).R;
                data.R.(neuralBand).(behavior).adjRH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjRH.All.(behavior).R;
                data.R2.(neuralBand).(behavior).adjLH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjLH.All.(behavior).R2;
                data.R2.(neuralBand).(behavior).adjRH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjRH.All.(behavior).R2;
                % all data kernel coff coef and R^2
                data.R.(neuralBand).AllData.(behavior).adjLH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjLH.All.(behavior).R;
                data.R.(neuralBand).AllData.(behavior).adjRH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjRH.All.(behavior).R;
                data.R2.(neuralBand).AllData.(behavior).adjLH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjLH.All.(behavior).R2;
                data.R2.(neuralBand).AllData.(behavior).adjRH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjRH.All.(behavior).R2;
                % coherence
                data.f.(neuralBand).(behavior).adjLH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjLH.All.(behavior).f;
                data.f.(neuralBand).(behavior).adjRH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjRH.All.(behavior).f;
                data.C.(neuralBand).(behavior).adjLH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjLH.All.(behavior).C;
                data.C.(neuralBand).(behavior).adjRH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjRH.All.(behavior).C;
                data.f.(neuralBand).AllData.(behavior).adjLH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjLH.All.(behavior).f;
                data.f.(neuralBand).AllData.(behavior).adjRH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjRH.All.(behavior).f;
                data.C.(neuralBand).AllData.(behavior).adjLH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjLH.All.(behavior).C;
                data.C.(neuralBand).AllData.(behavior).adjRH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjRH.All.(behavior).C;
            else
                % corr coeff and R^2
                data.R.(neuralBand).(behavior).adjLH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjLH.(behavior).R;
                data.R.(neuralBand).(behavior).adjRH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjRH.(behavior).R;
                data.R2.(neuralBand).(behavior).adjLH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjLH.(behavior).R2;
                data.R2.(neuralBand).(behavior).adjRH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjRH.(behavior).R2;
                % all data kernel coff coef and R^2
                data.R.(neuralBand).AllData.(behavior).adjLH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjLH.All.(behavior).R;
                data.R.(neuralBand).AllData.(behavior).adjRH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjRH.All.(behavior).R;
                data.R2.(neuralBand).AllData.(behavior).adjLH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjLH.All.(behavior).R2;
                data.R2.(neuralBand).AllData.(behavior).adjRH{aa,1} = AnalysisResults.Predictions.(animalID).HbT.(neuralBand).adjRH.All.(behavior).R2;
                % coherence
                data.f.(neuralBand).(behavior).adjLH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjLH.(behavior).f;
                data.f.(neuralBand).(behavior).adjRH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjRH.(behavior).f;
                data.C.(neuralBand).(behavior).adjLH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjLH.(behavior).C;
                data.C.(neuralBand).(behavior).adjRH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjRH.(behavior).C;
                data.f.(neuralBand).AllData.(behavior).adjLH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjLH.All.(behavior).f;
                data.f.(neuralBand).AllData.(behavior).adjRH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjRH.All.(behavior).f;
                data.C.(neuralBand).AllData.(behavior).adjLH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjLH.All.(behavior).C;
                data.C.(neuralBand).AllData.(behavior).adjRH(aa,:) = AnalysisResults.Coherence.(animalID).HbT.(neuralBand).adjRH.All.(behavior).C;
            end
        end
    end
    data.GammaHbT.gamma = cat(1,data.GammaHbT.gamma,AnalysisResults.GammaHbT.(animalID).LH_gamma,AnalysisResults.GammaHbT.(animalID).RH_gamma);
    data.GammaHbT.HbT = cat(1,data.GammaHbT.HbT,AnalysisResults.GammaHbT.(animalID).LH_HbT,AnalysisResults.GammaHbT.(animalID).RH_HbT);
end
% concatenate the data from left and right into a single data set
for dd = 1:length(neuralBands)
    neuralBand = neuralBands{1,dd};
    for ee = 1:length(behaviors)
        behavior = behaviors{1,ee};
        % impulse
        data.IR.(neuralBand).(behavior).cat = cat(1,data.IR.(neuralBand).(behavior).adjLH,data.IR.(neuralBand).(behavior).adjRH);
        data.IR_gamma.(neuralBand).(behavior).cat = cat(1,data.IR_gamma.(neuralBand).(behavior).adjLH,data.IR_gamma.(neuralBand).(behavior).adjRH);
        % corr coef and R^2
        data.R.(neuralBand).(behavior).cat = cat(1,data.R.(neuralBand).(behavior).adjLH,data.R.(neuralBand).(behavior).adjRH);
        data.R2.(neuralBand).(behavior).cat = cat(1,data.R2.(neuralBand).(behavior).adjLH,data.R2.(neuralBand).(behavior).adjRH);
        data.R.(neuralBand).AllData.(behavior).cat = cat(1,data.R.(neuralBand).AllData.(behavior).adjLH,data.R.(neuralBand).AllData.(behavior).adjRH);
        data.R2.(neuralBand).AllData.(behavior).cat = cat(1,data.R2.(neuralBand).AllData.(behavior).adjLH,data.R2.(neuralBand).AllData.(behavior).adjRH);
        % coherence
        data.f.(neuralBand).(behavior).cat = cat(1,data.f.(neuralBand).(behavior).adjLH,data.f.(neuralBand).(behavior).adjRH);
        data.C.(neuralBand).(behavior).cat = cat(1,data.C.(neuralBand).(behavior).adjLH,data.C.(neuralBand).(behavior).adjRH);
        data.f.(neuralBand).AllData.(behavior).cat = cat(1,data.f.(neuralBand).AllData.(behavior).adjLH,data.f.(neuralBand).AllData.(behavior).adjRH);
        data.C.(neuralBand).AllData.(behavior).cat = cat(1,data.C.(neuralBand).AllData.(behavior).adjLH,data.C.(neuralBand).AllData.(behavior).adjRH);
    end
end
% mean of individual hemispheres and cat ind points for histograms
for ff = 1:length(neuralBands)
    neuralBand = neuralBands{1,ff};
    for gg = 1:length(behaviors)
        behavior = behaviors{1,gg};
        % corr coef and R^2
        data.R.(neuralBand).(behavior).allData = [];
        data.R2.(neuralBand).(behavior).allData = [];
        % all data kernel corr coef and R^2
        data.R.(neuralBand).AllData.(behavior).allData = [];
        data.R2.(neuralBand).AllData.(behavior).allData = [];
        for hh = 1:length(data.R.(neuralBand).(behavior).cat)
            data.R.(neuralBand).(behavior).indMeans(hh,1) = mean(data.R.(neuralBand).(behavior).cat{hh,1});
            data.R.(neuralBand).(behavior).indMeds(hh,1) = median(data.R.(neuralBand).(behavior).cat{hh,1});
            data.R.(neuralBand).(behavior).allData = cat(2,data.R.(neuralBand).(behavior).allData,data.R.(neuralBand).(behavior).cat{hh,1});
            data.R2.(neuralBand).(behavior).indMeans(hh,1) = mean(data.R2.(neuralBand).(behavior).cat{hh,1});
            data.R2.(neuralBand).(behavior).indMeds(hh,1) = median(data.R2.(neuralBand).(behavior).cat{hh,1});
            data.R2.(neuralBand).(behavior).allData = cat(2,data.R2.(neuralBand).(behavior).allData,data.R2.(neuralBand).(behavior).cat{hh,1});
            data.R.(neuralBand).AllData.(behavior).indMeans(hh,1) = mean(data.R.(neuralBand).AllData.(behavior).cat{hh,1});
            data.R.(neuralBand).AllData.(behavior).indMeds(hh,1) = median(data.R.(neuralBand).AllData.(behavior).cat{hh,1});
            data.R.(neuralBand).AllData.(behavior).allData = cat(2,data.R.(neuralBand).AllData.(behavior).allData,data.R.(neuralBand).AllData.(behavior).cat{hh,1});
            data.R2.(neuralBand).AllData.(behavior).indMeans(hh,1) = mean(data.R2.(neuralBand).AllData.(behavior).cat{hh,1});
            data.R2.(neuralBand).AllData.(behavior).indMeds(hh,1) = median(data.R2.(neuralBand).AllData.(behavior).cat{hh,1});
            data.R2.(neuralBand).AllData.(behavior).allData = cat(2,data.R2.(neuralBand).AllData.(behavior).allData,data.R2.(neuralBand).AllData.(behavior).cat{hh,1});
        end
    end
end
% mean and std of each arousal-state
for ii = 1:length(neuralBands)
    neuralBand = neuralBands{1,ii};
    for jj = 1:length(behaviors)
        behavior = behaviors{1,jj};
        % IR mean
        data.IR.(neuralBand).(behavior).mean = mean(data.IR.(neuralBand).(behavior).cat,1);
        % R mean/std
        data.R.(neuralBand).(behavior).mean = mean(data.R.(neuralBand).(behavior).indMeans,1);
        data.R.(neuralBand).(behavior).meanStD = std(data.R.(neuralBand).(behavior).indMeans,0,1);
        % R med/std
        data.R.(neuralBand).(behavior).median = mean(data.R.(neuralBand).(behavior).indMeds,1);
        data.R.(neuralBand).(behavior).medStD = std(data.R.(neuralBand).(behavior).indMeds,0,1);
        % R^2 mean/std
        data.R2.(neuralBand).(behavior).mean = mean(data.R2.(neuralBand).(behavior).indMeans,1);
        data.R2.(neuralBand).(behavior).meanStD = std(data.R2.(neuralBand).(behavior).indMeans,0,1);
        % R^2 med/std
        data.R2.(neuralBand).(behavior).median = mean(data.R2.(neuralBand).(behavior).indMeds,1);
        data.R2.(neuralBand).(behavior).medStD = std(data.R2.(neuralBand).(behavior).indMeds,0,1);
        % all data kernel R mean/std
        data.R.(neuralBand).AllData.(behavior).mean = mean(data.R.(neuralBand).AllData.(behavior).indMeans,1);
        data.R.(neuralBand).AllData.(behavior).meanStD = std(data.R.(neuralBand).AllData.(behavior).indMeans,0,1);
        % all data kernel R med/std
        data.R.(neuralBand).AllData.(behavior).median = mean(data.R.(neuralBand).AllData.(behavior).indMeds,1);
        data.R.(neuralBand).AllData.(behavior).medStD = std(data.R.(neuralBand).AllData.(behavior).indMeds,0,1);
        % all data kernel R^2 mean/std
        data.R2.(neuralBand).AllData.(behavior).mean = mean(data.R2.(neuralBand).AllData.(behavior).indMeans,1);
        data.R2.(neuralBand).AllData.(behavior).meanStD = std(data.R2.(neuralBand).AllData.(behavior).indMeans,0,1);
        % all data kernel R^2 med/std
        data.R2.(neuralBand).AllData.(behavior).median = mean(data.R2.(neuralBand).AllData.(behavior).indMeds,1);
        data.R2.(neuralBand).AllData.(behavior).medStD = std(data.R2.(neuralBand).AllData.(behavior).indMeds,0,1);
        % coherence mean/std
        data.f.(neuralBand).(behavior).meanf = mean(data.f.(neuralBand).(behavior).cat,1);
        data.f.(neuralBand).(behavior).StDf = std(data.f.(neuralBand).(behavior).cat,0,1);
        data.C.(neuralBand).(behavior).meanC = mean(data.C.(neuralBand).(behavior).cat,1);
        data.C.(neuralBand).(behavior).StDC = std(data.C.(neuralBand).(behavior).cat,0,1);
        % all data kernel coherence mean/std
        data.f.(neuralBand).AllData.(behavior).meanf = mean(data.f.(neuralBand).AllData.(behavior).cat,1);
        data.f.(neuralBand).AllData.(behavior).StDf = std(data.f.(neuralBand).AllData.(behavior).cat,0,1);
        data.C.(neuralBand).AllData.(behavior).meanC = mean(data.C.(neuralBand).AllData.(behavior).cat,1);
        data.C.(neuralBand).AllData.(behavior).StDC = std(data.C.(neuralBand).AllData.(behavior).cat,0,1);
    end
end
%% gamma HRF based on impulse deconvolution
for kk = 1:length(neuralBands)
    neuralBand = neuralBands{1,kk};
    for ll = 1:length(behaviors)
        behavior = behaviors{1,ll};
        [peak,peakIndex] = max(data.IR.(neuralBand).(behavior).mean);
        peakTime = peakIndex/samplingRate;
        threeQuarterMax = max(data.IR.(neuralBand).(behavior).mean)/(4/3);
        index1 = find(data.IR.(neuralBand).(behavior).mean >= threeQuarterMax,1,'first');
        % find where the data last rises above half the max.
        index2 = find(data.IR.(neuralBand).(behavior).mean >= threeQuarterMax,1,'last');
        threeQuarterWidth = (index2 - index1 + 1)/samplingRate; % FWHM in indexes.
        initVals = [peak,peakTime,threeQuarterWidth];
        % create gamma function based on impulse values
        t = 0:1/samplingRate:5;
        IR_a = ((initVals(2)/initVals(3))^2*8*log10(2));
        IR_beta = ((initVals(3)^2)/initVals(2)/8/log10(2));
        data.IR_gamma.(neuralBand).(behavior).repFunc = initVals(1)*(t/initVals(2)).^IR_a.*exp((t - initVals(2))/(-1*IR_beta));
    end
end
%% Fig.1
figure('Name','Fig1 (-)');
sgtitle('Figure 1')
%% [1a] gamma-band impulse HRFs
subplot(2,2,1);
p1 = plot(t,data.IR.gammaBandPower.Rest.mean,'color',colorRest,'LineWidth',2);
hold on
p2 = plot(t,data.IR.gammaBandPower.Whisk.mean,'color',colorWhisk,'LineWidth',2);
p3 = plot(t,data.IR.gammaBandPower.Contra.mean,'color',colorStim,'LineWidth',2);
p4 = plot(t,data.IR.gammaBandPower.NREM.mean,'color',colorNREM,'LineWidth',2);
p5 = plot(t,data.IR.gammaBandPower.REM.mean,'color',colorREM,'LineWidth',2);
p6 = plot(t,data.IR.gammaBandPower.All.mean,'color',colorAll,'LineWidth',2);
title({'[1a] Impulse response response function','gamma-band [30-100 Hz] derived'})
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
legend([p1,p2,p3,p4,p5,p6],'Rest','Whisk','Stim','NREM','REM','All','Location','NorthEast')
axis tight
xlim([0,5])
set(gca,'box','off')
%% [1b] gamma-band impulse HRFs
subplot(2,2,2);
plot(t,data.IR_gamma.gammaBandPower.Rest.repFunc,'color',colorRest,'LineWidth',2);
hold on
plot(t,data.IR_gamma.gammaBandPower.Whisk.repFunc,'color',colorWhisk,'LineWidth',2);
plot(t,data.IR_gamma.gammaBandPower.Contra.repFunc,'color',colorStim,'LineWidth',2);
plot(t,data.IR_gamma.gammaBandPower.NREM.repFunc,'color',colorNREM,'LineWidth',2);
plot(t,data.IR_gamma.gammaBandPower.REM.repFunc,'color',colorREM,'LineWidth',2);
plot(t,data.IR_gamma.gammaBandPower.All.repFunc,'color',colorAll,'LineWidth',2);
title({'[1b] Representative gamma function','gamma-band [30-100 Hz] derived'})
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
axis tight
xlim([0,5])
set(gca,'box','off')
%% [1c] MUA impulse HRFs
subplot(2,2,3);
plot(t,data.IR.muaPower.Rest.mean,'color',colorRest,'LineWidth',2);
hold on
plot(t,data.IR.muaPower.Whisk.mean,'color',colorWhisk,'LineWidth',2);
plot(t,data.IR.muaPower.Contra.mean,'color',colorStim,'LineWidth',2);
plot(t,data.IR.muaPower.NREM.mean,'color',colorNREM,'LineWidth',2);
plot(t,data.IR.muaPower.REM.mean,'color',colorREM,'LineWidth',2);
plot(t,data.IR.muaPower.All.mean,'color',colorAll,'LineWidth',2);
title({'[1c] Impulse response response function','MUA [300-3000 Hz] derived'})
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
axis tight
xlim([0,5])
set(gca,'box','off')
%% [1d] MUA impulse HRFs
subplot(2,2,4);
plot(t,data.IR_gamma.muaPower.Rest.repFunc,'color',colorRest,'LineWidth',2);
hold on
plot(t,data.IR_gamma.muaPower.Whisk.repFunc,'color',colorWhisk,'LineWidth',2);
plot(t,data.IR_gamma.muaPower.Contra.repFunc,'color',colorStim,'LineWidth',2);
plot(t,data.IR_gamma.muaPower.NREM.repFunc,'color',colorNREM,'LineWidth',2);
plot(t,data.IR_gamma.muaPower.REM.repFunc,'color',colorREM,'LineWidth',2);
plot(t,data.IR_gamma.muaPower.All.repFunc,'color',colorAll,'LineWidth',2);
title({'[1d] Representative gamma function','MUA [300-3000 Hz] derived'})
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
axis tight
xlim([0,5])
set(gca,'box','off')
%% Fig.2
figure('Name','Fig2 (-)');
sgtitle('Figure 2')
%% [2a] gamma-band gamma HRF predictions - dedicated kernels
subplot(2,1,1);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.R.gammaBandPower.Rest.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.R.gammaBandPower.Rest.mean,data.R.gammaBandPower.Rest.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.R.gammaBandPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.R.gammaBandPower.Rest.median,data.R.gammaBandPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.R.gammaBandPower.Whisk.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.R.gammaBandPower.Whisk.mean,data.R.gammaBandPower.Whisk.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.R.gammaBandPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.R.gammaBandPower.Whisk.median,data.R.gammaBandPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.R.gammaBandPower.Contra.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.R.gammaBandPower.Contra.mean,data.R.gammaBandPower.Contra.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.R.gammaBandPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.R.gammaBandPower.Contra.median,data.R.gammaBandPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.R.gammaBandPower.NREM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.R.gammaBandPower.NREM.mean,data.R.gammaBandPower.NREM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.R.gammaBandPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.R.gammaBandPower.NREM.median,data.R.gammaBandPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.R.gammaBandPower.REM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.R.gammaBandPower.REM.mean,data.R.gammaBandPower.REM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.R.gammaBandPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.R.gammaBandPower.REM.median,data.R.gammaBandPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.R.gammaBandPower.All.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.R.gammaBandPower.All.mean,data.R.gammaBandPower.All.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.R.gammaBandPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.R.gammaBandPower.All.median,data.R.gammaBandPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title({'[2a] gamma function predictions (R)','gamma-band [30-100 Hz] derived'})
ylabel('Corr coef (R)')
xlabel('shaded = mean, lighter = median')
xlim([0,13])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% [2b] MUA impulse HRF predictions
subplot(2,1,2);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.R.muaPower.Rest.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.R.muaPower.Rest.mean,data.R.muaPower.Rest.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.R.muaPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.R.muaPower.Rest.median,data.R.muaPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.R.muaPower.Whisk.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.R.muaPower.Whisk.mean,data.R.muaPower.Whisk.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.R.muaPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.R.muaPower.Whisk.median,data.R.muaPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.R.muaPower.Contra.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.R.muaPower.Contra.mean,data.R.muaPower.Contra.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.R.muaPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.R.muaPower.Contra.median,data.R.muaPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.R.muaPower.NREM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.R.muaPower.NREM.mean,data.R.muaPower.NREM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.R.muaPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.R.muaPower.NREM.median,data.R.muaPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.R.muaPower.REM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.R.muaPower.REM.mean,data.R.muaPower.REM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.R.muaPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.R.muaPower.REM.median,data.R.muaPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.R.muaPower.All.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.R.muaPower.All.mean,data.R.muaPower.All.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.R.muaPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.R.muaPower.All.median,data.R.muaPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title({'[2b] gamma function predictions (R)','MUA [300-3000 Hz] derived'})
ylabel('Corr coef (R)')
xlabel('shaded = mean, lighter = median')
xlim([0,13])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% Fig.3
figure('Name','Fig3 (-)');
sgtitle('Figure 3')
%% [3a] gamma-band gamma HRF predictions
subplot(2,1,1);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.R2.gammaBandPower.Rest.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.R2.gammaBandPower.Rest.mean,data.R2.gammaBandPower.Rest.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.R2.gammaBandPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.R2.gammaBandPower.Rest.median,data.R2.gammaBandPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.R2.gammaBandPower.Whisk.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.R2.gammaBandPower.Whisk.mean,data.R2.gammaBandPower.Whisk.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.R2.gammaBandPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.R2.gammaBandPower.Whisk.median,data.R2.gammaBandPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.R2.gammaBandPower.Contra.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.R2.gammaBandPower.Contra.mean,data.R2.gammaBandPower.Contra.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.R2.gammaBandPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.R2.gammaBandPower.Contra.median,data.R2.gammaBandPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.R2.gammaBandPower.NREM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.R2.gammaBandPower.NREM.mean,data.R2.gammaBandPower.NREM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.R2.gammaBandPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.R2.gammaBandPower.NREM.median,data.R2.gammaBandPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.R2.gammaBandPower.REM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.R2.gammaBandPower.REM.mean,data.R2.gammaBandPower.REM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.R2.gammaBandPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.R2.gammaBandPower.REM.median,data.R2.gammaBandPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.R2.gammaBandPower.All.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.R2.gammaBandPower.All.mean,data.R2.gammaBandPower.All.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.R2.gammaBandPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.R2.gammaBandPower.All.median,data.R2.gammaBandPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title({'[3a] gamma function predictions (R^2)','gamma-band [30-100 Hz] derived'})
ylabel('Coeff of det. (R^2)')
xlabel('shaded = mean, lighter = median')
xlim([0,13])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% [3b] MUA impulse HRF predictions
subplot(2,1,2);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.R2.muaPower.Rest.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.R2.muaPower.Rest.mean,data.R2.muaPower.Rest.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.R2.muaPower.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.R2.muaPower.Rest.median,data.R2.muaPower.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.R2.muaPower.Whisk.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.R2.muaPower.Whisk.mean,data.R2.muaPower.Whisk.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.R2.muaPower.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.R2.muaPower.Whisk.median,data.R2.muaPower.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.R2.muaPower.Contra.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.R2.muaPower.Contra.mean,data.R2.muaPower.Contra.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.R2.muaPower.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.R2.muaPower.Contra.median,data.R2.muaPower.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.R2.muaPower.NREM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.R2.muaPower.NREM.mean,data.R2.muaPower.NREM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.R2.muaPower.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.R2.muaPower.NREM.median,data.R2.muaPower.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.R2.muaPower.REM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.R2.muaPower.REM.mean,data.R2.muaPower.REM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.R2.muaPower.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.R2.muaPower.REM.median,data.R2.muaPower.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.R2.muaPower.All.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.R2.muaPower.All.mean,data.R2.muaPower.All.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.R2.muaPower.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.R2.muaPower.All.median,data.R2.muaPower.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title({'[3b] gamma function predictions (R^2)','MUA [300-3000 Hz] derived'})
ylabel('Coeff of det. (R^2)')
xlabel('shaded = mean, lighter = median')
xlim([0,13])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% Fig.4
figure('Name','Fig4 (-)');
sgtitle('Figure 4')
%% [4a] gamma-band gamma HRF predictions (R) distribution
subplot(2,2,1);
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.Rest.allData,21,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.Whisk.allData,21,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.Contra.allData,21,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.NREM.allData,21,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.REM.allData,21,'kernel');
[xCurve6,yCurve6] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.All.allData,21,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
plot(xCurve6,yCurve6,'color',colorAll,'LineWidth',2)
title({'[4a] gamma function predictions (R)','gamma-band [30-100 Hz] derived'})
xlabel('Corr coeff (R)')
ylabel('Probability')
set(gca,'box','off')
xlim([-1,1.15])
%% [4b] gamma-band gamma HRF predictions (R^2) distribution
subplot(2,2,2);
% set vals = 0 for below 0
for aa = 1:length(data.R2.gammaBandPower.Rest.allData)
    if data.R2.gammaBandPower.Rest.allData(1,aa) < 0
        data.R2.gammaBandPower.Rest.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.Whisk.allData)
    if data.R2.gammaBandPower.Whisk.allData(1,aa) < 0
        data.R2.gammaBandPower.Whisk.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.Contra.allData)
    if data.R2.gammaBandPower.Contra.allData(1,aa) < 0
        data.R2.gammaBandPower.Contra.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.NREM.allData)
    if data.R2.gammaBandPower.NREM.allData(1,aa) < 0
        data.R2.gammaBandPower.NREM.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.REM.allData)
    if data.R2.gammaBandPower.REM.allData(1,aa) < 0
        data.R2.gammaBandPower.REM.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.All.allData)
    if data.R2.gammaBandPower.All.allData(1,aa) < 0
        data.R2.gammaBandPower.All.allData(1,aa) = 0;
    end
end
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.Rest.allData,21,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.Whisk.allData,21,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.Contra.allData,21,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.NREM.allData,21,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.REM.allData,21,'kernel');
[xCurve6,yCurve6] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.All.allData,21,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
plot(xCurve6,yCurve6,'color',colorAll,'LineWidth',2)
title({'[4b] gamma function predictions (R^2)','gamma-band [30-100 Hz] derived'})
xlabel('Coeff of det. (R^2)')
ylabel('Probability')
set(gca,'box','off')
xlim([-0.15,1.15])
%% [4c] MUA gamma HRF predictions (R) distribution
subplot(2,2,3);
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.Rest.allData,21,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.Whisk.allData,21,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.Contra.allData,21,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.NREM.allData,21,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.REM.allData,21,'kernel');
[xCurve6,yCurve6] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.All.allData,21,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
plot(xCurve6,yCurve6,'color',colorAll,'LineWidth',2)
title({'[4c] gamma function predictions (R)','MUA [300-3000 Hz] derived'})
xlabel('Corr coeff (R)')
ylabel('Probability')
set(gca,'box','off')
xlim([-1,1.15])
%% [4d] MUA gamma HRF predictions (R^2) distribution
subplot(2,2,4);
% set vals = 0 for below 0
for aa = 1:length(data.R2.muaPower.Rest.allData)
    if data.R2.muaPower.Rest.allData(1,aa) < 0
        data.R2.muaPower.Rest.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.Whisk.allData)
    if data.R2.muaPower.Whisk.allData(1,aa) < 0
        data.R2.muaPower.Whisk.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.Contra.allData)
    if data.R2.muaPower.Contra.allData(1,aa) < 0
        data.R2.muaPower.Contra.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.NREM.allData)
    if data.R2.muaPower.NREM.allData(1,aa) < 0
        data.R2.muaPower.NREM.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.REM.allData)
    if data.R2.muaPower.REM.allData(1,aa) < 0
        data.R2.muaPower.REM.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.All.allData)
    if data.R2.muaPower.All.allData(1,aa) < 0
        data.R2.muaPower.All.allData(1,aa) = 0;
    end
end
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.Rest.allData,21,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.Whisk.allData,21,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.Contra.allData,21,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.NREM.allData,21,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.REM.allData,21,'kernel');
[xCurve6,yCurve6] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.All.allData,21,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
plot(xCurve6,yCurve6,'color',colorAll,'LineWidth',2)
title({'[4d] gamma function predictions (R^2)','MUA [300-3000 Hz] derived'})
xlabel('Coeff of det. (R^2)')
ylabel('Probability')
set(gca,'box','off')
xlim([-0.15,1.15])
%% Fig.5
figure('Name','Fig5 (-)');
sgtitle('Figure 5')
%% [5a] gamma-band gamma HRF predictions - dedicated kernels
subplot(2,1,1);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.R.gammaBandPower.AllData.Rest.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.R.gammaBandPower.AllData.Rest.mean,data.R.gammaBandPower.AllData.Rest.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.R.gammaBandPower.AllData.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.R.gammaBandPower.AllData.Rest.median,data.R.gammaBandPower.AllData.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.R.gammaBandPower.AllData.Whisk.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.R.gammaBandPower.AllData.Whisk.mean,data.R.gammaBandPower.AllData.Whisk.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.R.gammaBandPower.AllData.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.R.gammaBandPower.AllData.Whisk.median,data.R.gammaBandPower.AllData.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.R.gammaBandPower.AllData.Contra.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.R.gammaBandPower.AllData.Contra.mean,data.R.gammaBandPower.AllData.Contra.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.R.gammaBandPower.AllData.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.R.gammaBandPower.AllData.Contra.median,data.R.gammaBandPower.AllData.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.R.gammaBandPower.AllData.NREM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.R.gammaBandPower.AllData.NREM.mean,data.R.gammaBandPower.AllData.NREM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.R.gammaBandPower.AllData.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.R.gammaBandPower.AllData.NREM.median,data.R.gammaBandPower.AllData.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.R.gammaBandPower.AllData.REM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.R.gammaBandPower.AllData.REM.mean,data.R.gammaBandPower.AllData.REM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.R.gammaBandPower.AllData.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.R.gammaBandPower.AllData.REM.median,data.R.gammaBandPower.AllData.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.R.gammaBandPower.AllData.All.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.R.gammaBandPower.AllData.All.mean,data.R.gammaBandPower.AllData.All.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.R.gammaBandPower.AllData.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.R.gammaBandPower.AllData.All.median,data.R.gammaBandPower.AllData.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title({'[5a] gamma function predictions (R)','gamma-band [30-100 Hz] derived'})
ylabel('Corr coef (R)')
xlabel('shaded = mean, lighter = median')
xlim([0,13])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% [5b] MUA impulse HRF predictions
subplot(2,1,2);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.R.muaPower.AllData.Rest.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.R.muaPower.AllData.Rest.mean,data.R.muaPower.AllData.Rest.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.R.muaPower.AllData.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.R.muaPower.AllData.Rest.median,data.R.muaPower.AllData.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.R.muaPower.AllData.Whisk.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.R.muaPower.AllData.Whisk.mean,data.R.muaPower.AllData.Whisk.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.R.muaPower.AllData.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.R.muaPower.AllData.Whisk.median,data.R.muaPower.AllData.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.R.muaPower.AllData.Contra.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.R.muaPower.AllData.Contra.mean,data.R.muaPower.AllData.Contra.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.R.muaPower.AllData.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.R.muaPower.AllData.Contra.median,data.R.muaPower.AllData.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.R.muaPower.AllData.NREM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.R.muaPower.AllData.NREM.mean,data.R.muaPower.AllData.NREM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.R.muaPower.AllData.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.R.muaPower.AllData.NREM.median,data.R.muaPower.AllData.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.R.muaPower.AllData.REM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.R.muaPower.AllData.REM.mean,data.R.muaPower.AllData.REM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.R.muaPower.AllData.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.R.muaPower.AllData.REM.median,data.R.muaPower.AllData.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.R.muaPower.AllData.All.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.R.muaPower.AllData.All.mean,data.R.muaPower.AllData.All.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.R.muaPower.AllData.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.R.muaPower.AllData.All.median,data.R.muaPower.AllData.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title({'[5b] gamma function predictions (R)','MUA [300-3000 Hz] derived'})
ylabel('Corr coef (R)')
xlabel('shaded = mean, lighter = median')
xlim([0,13])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% Fig.6
figure('Name','Fig6 (-)');
sgtitle('Figure 6')
%% [6a] gamma-band gamma HRF predictions
subplot(2,1,1);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.R2.gammaBandPower.AllData.Rest.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.R2.gammaBandPower.AllData.Rest.mean,data.R2.gammaBandPower.AllData.Rest.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.R2.gammaBandPower.AllData.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.R2.gammaBandPower.AllData.Rest.median,data.R2.gammaBandPower.AllData.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.R2.gammaBandPower.AllData.Whisk.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.R2.gammaBandPower.AllData.Whisk.mean,data.R2.gammaBandPower.AllData.Whisk.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.R2.gammaBandPower.AllData.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.R2.gammaBandPower.AllData.Whisk.median,data.R2.gammaBandPower.AllData.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.R2.gammaBandPower.AllData.Contra.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.R2.gammaBandPower.AllData.Contra.mean,data.R2.gammaBandPower.AllData.Contra.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.R2.gammaBandPower.AllData.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.R2.gammaBandPower.AllData.Contra.median,data.R2.gammaBandPower.AllData.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.R2.gammaBandPower.AllData.NREM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.R2.gammaBandPower.AllData.NREM.mean,data.R2.gammaBandPower.AllData.NREM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.R2.gammaBandPower.AllData.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.R2.gammaBandPower.AllData.NREM.median,data.R2.gammaBandPower.AllData.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.R2.gammaBandPower.AllData.REM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.R2.gammaBandPower.AllData.REM.mean,data.R2.gammaBandPower.AllData.REM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.R2.gammaBandPower.AllData.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.R2.gammaBandPower.AllData.REM.median,data.R2.gammaBandPower.AllData.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.R2.gammaBandPower.AllData.All.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.R2.gammaBandPower.AllData.All.mean,data.R2.gammaBandPower.AllData.All.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.R2.gammaBandPower.AllData.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.R2.gammaBandPower.AllData.All.median,data.R2.gammaBandPower.AllData.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title({'[6a] gamma function predictions (R^2)','gamma-band [30-100 Hz] derived'})
ylabel('Coeff of det. (R^2)')
xlabel('shaded = mean, lighter = median')
xlim([0,13])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% [6b] MUA impulse HRF predictions
subplot(2,1,2);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.R2.muaPower.AllData.Rest.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.R2.muaPower.AllData.Rest.mean,data.R2.muaPower.AllData.Rest.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.R2.muaPower.AllData.Rest.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorRest)*0.5 + colorRest,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.R2.muaPower.AllData.Rest.median,data.R2.muaPower.AllData.Rest.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.R2.muaPower.AllData.Whisk.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.R2.muaPower.AllData.Whisk.mean,data.R2.muaPower.AllData.Whisk.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.R2.muaPower.AllData.Whisk.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorWhisk)*0.5 + colorWhisk,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.R2.muaPower.AllData.Whisk.median,data.R2.muaPower.AllData.Whisk.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.R2.muaPower.AllData.Contra.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.R2.muaPower.AllData.Contra.mean,data.R2.muaPower.AllData.Contra.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(xInds*6,data.R2.muaPower.AllData.Contra.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorStim)*0.5 + colorStim,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.R2.muaPower.AllData.Contra.median,data.R2.muaPower.AllData.Contra.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(xInds*7,data.R2.muaPower.AllData.NREM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,data.R2.muaPower.AllData.NREM.mean,data.R2.muaPower.AllData.NREM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(xInds*8,data.R2.muaPower.AllData.NREM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorNREM)*0.5 + colorNREM,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,data.R2.muaPower.AllData.NREM.median,data.R2.muaPower.AllData.NREM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(xInds*9,data.R2.muaPower.AllData.REM.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,data.R2.muaPower.AllData.REM.mean,data.R2.muaPower.AllData.REM.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(xInds*10,data.R2.muaPower.AllData.REM.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorREM)*0.5 + colorREM,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,data.R2.muaPower.AllData.REM.median,data.R2.muaPower.AllData.REM.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
scatter(xInds*11,data.R2.muaPower.AllData.All.indMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,data.R2.muaPower.AllData.All.mean,data.R2.muaPower.AllData.All.meanStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(xInds*12,data.R2.muaPower.AllData.All.indMeds,75,'MarkerEdgeColor','k','MarkerFaceColor',(1 - colorAll)*0.5 + colorAll,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,data.R2.muaPower.AllData.All.median,data.R2.muaPower.AllData.All.medStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
title({'[6b] gamma function predictions (R^2)','MUA [300-3000 Hz] derived'})
ylabel('Coeff of det. (R^2)')
xlabel('shaded = mean, lighter = median')
xlim([0,13])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'box','off')
%% Fig.7
figure('Name','Fig7 (-)');
sgtitle('Figure 7')
%% [7a] gamma-band gamma HRF predictions (R) distribution
subplot(2,2,1);
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.AllData.Rest.allData,21,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.AllData.Whisk.allData,21,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.AllData.Contra.allData,21,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.AllData.NREM.allData,21,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.AllData.REM.allData,21,'kernel');
[xCurve6,yCurve6] = SmoothHistogramBinsFit_eLife2020(data.R.gammaBandPower.AllData.All.allData,21,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
plot(xCurve6,yCurve6,'color',colorAll,'LineWidth',2)
title({'[7a] gamma function predictions (R)','gamma-band [30-100 Hz] derived'})
xlabel('Corr coeff (R)')
ylabel('Probability')
set(gca,'box','off')
xlim([-1,1.15])
%% [7b] gamma-band gamma HRF predictions (R^2) distribution
subplot(2,2,2);
% set vals = 0 for below 0
for aa = 1:length(data.R2.gammaBandPower.AllData.Rest.allData)
    if data.R2.gammaBandPower.AllData.Rest.allData(1,aa) < 0
        data.R2.gammaBandPower.AllData.Rest.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.AllData.Whisk.allData)
    if data.R2.gammaBandPower.AllData.Whisk.allData(1,aa) < 0
        data.R2.gammaBandPower.AllData.Whisk.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.AllData.Contra.allData)
    if data.R2.gammaBandPower.AllData.Contra.allData(1,aa) < 0
        data.R2.gammaBandPower.AllData.Contra.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.AllData.NREM.allData)
    if data.R2.gammaBandPower.AllData.NREM.allData(1,aa) < 0
        data.R2.gammaBandPower.AllData.NREM.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.AllData.REM.allData)
    if data.R2.gammaBandPower.AllData.REM.allData(1,aa) < 0
        data.R2.gammaBandPower.AllData.REM.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.gammaBandPower.AllData.All.allData)
    if data.R2.gammaBandPower.AllData.All.allData(1,aa) < 0
        data.R2.gammaBandPower.AllData.All.allData(1,aa) = 0;
    end
end
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.AllData.Rest.allData,21,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.AllData.Whisk.allData,21,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.AllData.Contra.allData,21,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.AllData.NREM.allData,21,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.AllData.REM.allData,21,'kernel');
[xCurve6,yCurve6] = SmoothHistogramBinsFit_eLife2020(data.R2.gammaBandPower.AllData.All.allData,21,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
plot(xCurve6,yCurve6,'color',colorAll,'LineWidth',2)
title({'[7b] gamma function predictions (R^2)','gamma-band [30-100 Hz] derived'})
xlabel('Coeff of det. (R^2)')
ylabel('Probability')
set(gca,'box','off')
xlim([-0.15,1.15])
%% [7c] MUA gamma HRF predictions (R) distribution
subplot(2,2,3);
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.AllData.Rest.allData,21,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.AllData.Whisk.allData,21,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.AllData.Contra.allData,21,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.AllData.NREM.allData,21,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.AllData.REM.allData,21,'kernel');
[xCurve6,yCurve6] = SmoothHistogramBinsFit_eLife2020(data.R.muaPower.AllData.All.allData,21,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
plot(xCurve6,yCurve6,'color',colorAll,'LineWidth',2)
title({'[7c] gamma function predictions (R)','MUA [300-3000 Hz] derived'})
xlabel('Corr coeff (R)')
ylabel('Probability')
set(gca,'box','off')
xlim([-1,1.15])
%% [7d] MUA gamma HRF predictions (R^2) distribution
subplot(2,2,4);
% set vals = 0 for below 0
for aa = 1:length(data.R2.muaPower.AllData.Rest.allData)
    if data.R2.muaPower.AllData.Rest.allData(1,aa) < 0
        data.R2.muaPower.AllData.Rest.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.AllData.Whisk.allData)
    if data.R2.muaPower.AllData.Whisk.allData(1,aa) < 0
        data.R2.muaPower.AllData.Whisk.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.AllData.Contra.allData)
    if data.R2.muaPower.AllData.Contra.allData(1,aa) < 0
        data.R2.muaPower.AllData.Contra.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.AllData.NREM.allData)
    if data.R2.muaPower.AllData.NREM.allData(1,aa) < 0
        data.R2.muaPower.AllData.NREM.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.AllData.REM.allData)
    if data.R2.muaPower.AllData.REM.allData(1,aa) < 0
        data.R2.muaPower.AllData.REM.allData(1,aa) = 0;
    end
end
for aa = 1:length(data.R2.muaPower.AllData.All.allData)
    if data.R2.muaPower.AllData.All.allData(1,aa) < 0
        data.R2.muaPower.AllData.All.allData(1,aa) = 0;
    end
end
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.AllData.Rest.allData,21,'kernel');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.AllData.Whisk.allData,21,'kernel');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.AllData.Contra.allData,21,'kernel');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.AllData.NREM.allData,21,'kernel');
[xCurve5,yCurve5] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.AllData.REM.allData,21,'kernel');
[xCurve6,yCurve6] = SmoothHistogramBinsFit_eLife2020(data.R2.muaPower.AllData.All.allData,21,'kernel');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
plot(xCurve6,yCurve6,'color',colorAll,'LineWidth',2)
title({'[7d] gamma function predictions (R^2)','MUA [300-3000 Hz] derived'})
xlabel('Coeff of det. (R^2)')
ylabel('Probability')
set(gca,'box','off')
xlim([-0.15,1.15])
%% Fig.8
figure('Name','Fig8 (-)');
sgtitle('Figure 8')
%% [8a] gamma-band kernel coherence 
subplot(2,2,1)
semilogx(data.f.gammaBandPower.Rest.meanf,data.C.gammaBandPower.Rest.meanC,'color',colorRest,'LineWidth',2)
hold on
semilogx(data.f.gammaBandPower.Whisk.meanf,data.C.gammaBandPower.Whisk.meanC,'color',colorWhisk,'LineWidth',2)
semilogx(data.f.gammaBandPower.Contra.meanf,data.C.gammaBandPower.Contra.meanC,'color',colorStim,'LineWidth',2)
semilogx(data.f.gammaBandPower.NREM.meanf,data.C.gammaBandPower.NREM.meanC,'color',colorNREM,'LineWidth',2)
semilogx(data.f.gammaBandPower.REM.meanf,data.C.gammaBandPower.REM.meanC,'color',colorREM,'LineWidth',2)
semilogx(data.f.gammaBandPower.All.meanf,data.C.gammaBandPower.All.meanC,'color',colorAll,'LineWidth',2)
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[8a] Gamma-band coherence','Actual vs. predicted \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% [8b] gamma-band all data kernel coherence
subplot(2,2,2)
semilogx(data.f.gammaBandPower.AllData.Rest.meanf,data.C.gammaBandPower.AllData.Rest.meanC,'color',colorRest,'LineWidth',2)
hold on
semilogx(data.f.gammaBandPower.AllData.Whisk.meanf,data.C.gammaBandPower.AllData.Whisk.meanC,'color',colorWhisk,'LineWidth',2)
semilogx(data.f.gammaBandPower.AllData.Contra.meanf,data.C.gammaBandPower.AllData.Contra.meanC,'color',colorStim,'LineWidth',2)
semilogx(data.f.gammaBandPower.AllData.NREM.meanf,data.C.gammaBandPower.AllData.NREM.meanC,'color',colorNREM,'LineWidth',2)
semilogx(data.f.gammaBandPower.AllData.REM.meanf,data.C.gammaBandPower.AllData.REM.meanC,'color',colorREM,'LineWidth',2)
semilogx(data.f.gammaBandPower.AllData.All.meanf,data.C.gammaBandPower.AllData.All.meanC,'color',colorAll,'LineWidth',2)
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[8a] Gamma-band all data kernel coherence','Actual vs. predicted \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% [8c] MUA kernel coherence
subplot(2,2,3)
semilogx(data.f.muaPower.Rest.meanf,data.C.muaPower.Rest.meanC,'color',colorRest,'LineWidth',2)
hold on
semilogx(data.f.muaPower.Whisk.meanf,data.C.muaPower.Whisk.meanC,'color',colorWhisk,'LineWidth',2)
semilogx(data.f.muaPower.Contra.meanf,data.C.muaPower.Contra.meanC,'color',colorStim,'LineWidth',2)
semilogx(data.f.muaPower.NREM.meanf,data.C.muaPower.NREM.meanC,'color',colorNREM,'LineWidth',2)
semilogx(data.f.muaPower.REM.meanf,data.C.muaPower.REM.meanC,'color',colorREM,'LineWidth',2)
semilogx(data.f.muaPower.All.meanf,data.C.muaPower.All.meanC,'color',colorAll,'LineWidth',2)
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[8c] MUA coherence','Actual vs. predicted \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% [8d] MUA all data kernel coherence
subplot(2,2,4)
semilogx(data.f.muaPower.AllData.Rest.meanf,data.C.muaPower.AllData.Rest.meanC,'color',colorRest,'LineWidth',2)
hold on
semilogx(data.f.muaPower.AllData.Whisk.meanf,data.C.muaPower.AllData.Whisk.meanC,'color',colorWhisk,'LineWidth',2)
semilogx(data.f.muaPower.AllData.Contra.meanf,data.C.muaPower.AllData.Contra.meanC,'color',colorStim,'LineWidth',2)
semilogx(data.f.muaPower.AllData.NREM.meanf,data.C.muaPower.AllData.NREM.meanC,'color',colorNREM,'LineWidth',2)
semilogx(data.f.muaPower.AllData.REM.meanf,data.C.muaPower.AllData.REM.meanC,'color',colorREM,'LineWidth',2)
semilogx(data.f.muaPower.AllData.All.meanf,data.C.muaPower.AllData.All.meanC,'color',colorAll,'LineWidth',2)
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[8d] MUA all data kernel coherence','Actual vs. predicted \Delta[HbT] (\muM)'})
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
%% Fig.9
figure('Name','Fig9 (-)');
%% [9a] 
histogram2(data.GammaHbT.gamma,data.GammaHbT.HbT,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',0:0.2:15,'YBinedges',-30:2.5:70,'Normalization','probability');
title({'[9a] Gamma-band vs. \Delta[HbT] (\muM)','max neural vs. max hemo after stimulus'})
xlabel('Gamma-band power [30-100 Hz]')
ylabel('\Delta[HbT] (\muM)')
set(gca,'box','off')
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
