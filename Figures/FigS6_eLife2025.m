function [] = FigS6_eLife2025(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% Ephys Pearsons correlations
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PearsonCorr_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','SSP_SAP','Blank_SAP'};
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PearsonCorr_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            ephysCorrData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(ephysCorrData.(group).(dataType),behavior) == false
                    ephysCorrData.(group).(dataType).(behavior).data = [];
                    ephysCorrData.(group).(dataType).(behavior).group = {};
                    ephysCorrData.(group).(dataType).(behavior).animalID = {};
                end
                if isempty(Results_PearsonCorr_Ephys.(group).(animalID).(dataType).(behavior).R) == false
                    ephysCorrData.(group).(dataType).(behavior).data = cat(1,ephysCorrData.(group).(dataType).(behavior).data,mean(Results_PearsonCorr_Ephys.(group).(animalID).(dataType).(behavior).R));
                    ephysCorrData.(group).(dataType).(behavior).group = cat(1,ephysCorrData.(group).(dataType).(behavior).group,group);
                    ephysCorrData.(group).(dataType).(behavior).animalID = cat(1,ephysCorrData.(group).(dataType).(behavior).animalID,animalID);
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            ephysCorrData.(group).(dataType).(behavior).meanData = mean(ephysCorrData.(group).(dataType).(behavior).data,1);
            ephysCorrData.(group).(dataType).(behavior).stdData = std(ephysCorrData.(group).(dataType).(behavior).data,0,1);
        end
    end
end
% statistics - ttest
for bb = 1:length(dataTypes)
    dataType = dataTypes{1,bb};
    for cc = 1:length(behaviors)
        behavior = behaviors{1,cc};
        % statistics - unpaired ttest
        [EphysPearsonStats.(dataType).(behavior).h,EphysPearsonStats.(dataType).(behavior).p] = ttest2(ephysCorrData.Blank_SAP.(dataType).(behavior).data,ephysCorrData.SSP_SAP.(dataType).(behavior).data);
    end
end

%% GCaMP Pearsons correlations
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PearsonCorr_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'SSP_SAP','Blank_SAP'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PearsonCorr_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            gcampCorrData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(gcampCorrData.(group).(dataType),behavior) == false
                    gcampCorrData.(group).(dataType).(behavior).data = [];
                    gcampCorrData.(group).(dataType).(behavior).group = {};
                    gcampCorrData.(group).(dataType).(behavior).animalID = {};
                end
                if isempty(Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).(behavior).R) == false
                    gcampCorrData.(group).(dataType).(behavior).data = cat(1,gcampCorrData.(group).(dataType).(behavior).data,mean(Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).(behavior).R));
                    gcampCorrData.(group).(dataType).(behavior).group = cat(1,gcampCorrData.(group).(dataType).(behavior).group,group);
                    gcampCorrData.(group).(dataType).(behavior).animalID = cat(1,gcampCorrData.(group).(dataType).(behavior).animalID,animalID);
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            gcampCorrData.(group).(dataType).(behavior).meanData = mean(gcampCorrData.(group).(dataType).(behavior).data,1);
            gcampCorrData.(group).(dataType).(behavior).stdData = std(gcampCorrData.(group).(dataType).(behavior).data,0,1);
        end
    end
end
% statistics - ttest
for bb = 1:length(dataTypes)
    dataType = dataTypes{1,bb};
    for cc = 1:length(behaviors)
        behavior = behaviors{1,cc};
        % statistics - unpaired ttest
        [GCaMPPearsonStats.(dataType).(behavior).h,GCaMPPearsonStats.(dataType).(behavior).p] = ttest2(gcampCorrData.Blank_SAP.(dataType).(behavior).data,gcampCorrData.SSP_SAP.(dataType).(behavior).data);
    end
end

%% Figure S6
FigS6 = figure('Name','Fig. S6');

subplot(1,4,1)
xInds = ones(1,length(ephysCorrData.Blank_SAP.HbT.Rest.data));
scatter(xInds*1,ephysCorrData.Blank_SAP.HbT.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(1,ephysCorrData.Blank_SAP.HbT.Rest.meanData,ephysCorrData.Blank_SAP.HbT.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.HbT.Rest.data));
scatter(xInds*2,ephysCorrData.SSP_SAP.HbT.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,ephysCorrData.SSP_SAP.HbT.Rest.meanData,ephysCorrData.SSP_SAP.HbT.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysCorrData.Blank_SAP.HbT.Alert.data));
scatter(xInds*4,ephysCorrData.Blank_SAP.HbT.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(4,ephysCorrData.Blank_SAP.HbT.Alert.meanData,ephysCorrData.Blank_SAP.HbT.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.HbT.Alert.data));
scatter(xInds*5,ephysCorrData.SSP_SAP.HbT.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,ephysCorrData.SSP_SAP.HbT.Alert.meanData,ephysCorrData.SSP_SAP.HbT.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysCorrData.Blank_SAP.HbT.Asleep.data));
scatter(xInds*6,ephysCorrData.Blank_SAP.HbT.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(6,ephysCorrData.Blank_SAP.HbT.Asleep.meanData,ephysCorrData.Blank_SAP.HbT.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.HbT.Asleep.data));
scatter(xInds*7,ephysCorrData.SSP_SAP.HbT.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(7,ephysCorrData.SSP_SAP.HbT.Asleep.meanData,ephysCorrData.SSP_SAP.HbT.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('Corr Coef')
xlim([0,8])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(1,4,2)
xInds = ones(1,length(ephysCorrData.Blank_SAP.gammaBandPower.Rest.data));
scatter(xInds*1,ephysCorrData.Blank_SAP.gammaBandPower.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(1,ephysCorrData.Blank_SAP.gammaBandPower.Rest.meanData,ephysCorrData.Blank_SAP.gammaBandPower.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.gammaBandPower.Rest.data));
scatter(xInds*2,ephysCorrData.SSP_SAP.gammaBandPower.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,ephysCorrData.SSP_SAP.gammaBandPower.Rest.meanData,ephysCorrData.SSP_SAP.gammaBandPower.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysCorrData.Blank_SAP.gammaBandPower.Alert.data));
scatter(xInds*4,ephysCorrData.Blank_SAP.gammaBandPower.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(4,ephysCorrData.Blank_SAP.gammaBandPower.Alert.meanData,ephysCorrData.Blank_SAP.gammaBandPower.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.gammaBandPower.Alert.data));
scatter(xInds*5,ephysCorrData.SSP_SAP.gammaBandPower.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,ephysCorrData.SSP_SAP.gammaBandPower.Alert.meanData,ephysCorrData.SSP_SAP.gammaBandPower.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysCorrData.Blank_SAP.gammaBandPower.Asleep.data));
scatter(xInds*6,ephysCorrData.Blank_SAP.gammaBandPower.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(6,ephysCorrData.Blank_SAP.gammaBandPower.Asleep.meanData,ephysCorrData.Blank_SAP.gammaBandPower.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysCorrData.SSP_SAP.gammaBandPower.Asleep.data));
scatter(xInds*7,ephysCorrData.SSP_SAP.gammaBandPower.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(7,ephysCorrData.SSP_SAP.gammaBandPower.Asleep.meanData,ephysCorrData.SSP_SAP.gammaBandPower.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('Corr Coef')
xlim([0,8])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(1,4,3)
xInds = ones(1,length(gcampCorrData.Blank_SAP.HbT.Rest.data));
scatter(xInds*1,gcampCorrData.Blank_SAP.HbT.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(1,gcampCorrData.Blank_SAP.HbT.Rest.meanData,gcampCorrData.Blank_SAP.HbT.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.HbT.Rest.data));
scatter(xInds*2,gcampCorrData.SSP_SAP.HbT.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,gcampCorrData.SSP_SAP.HbT.Rest.meanData,gcampCorrData.SSP_SAP.HbT.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampCorrData.Blank_SAP.HbT.Alert.data));
scatter(xInds*4,gcampCorrData.Blank_SAP.HbT.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(4,gcampCorrData.Blank_SAP.HbT.Alert.meanData,gcampCorrData.Blank_SAP.HbT.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.HbT.Alert.data));
scatter(xInds*5,gcampCorrData.SSP_SAP.HbT.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampCorrData.SSP_SAP.HbT.Alert.meanData,gcampCorrData.SSP_SAP.HbT.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampCorrData.Blank_SAP.HbT.Asleep.data));
scatter(xInds*6,gcampCorrData.Blank_SAP.HbT.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(6,gcampCorrData.Blank_SAP.HbT.Asleep.meanData,gcampCorrData.Blank_SAP.HbT.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.HbT.Asleep.data));
scatter(xInds*7,gcampCorrData.SSP_SAP.HbT.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(7,gcampCorrData.SSP_SAP.HbT.Asleep.meanData,gcampCorrData.SSP_SAP.HbT.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('Corr Coef')
xlim([0,8])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(1,4,4)
xInds = ones(1,length(gcampCorrData.Blank_SAP.GCaMP.Rest.data));
scatter(xInds*1,gcampCorrData.Blank_SAP.GCaMP.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(1,gcampCorrData.Blank_SAP.GCaMP.Rest.meanData,gcampCorrData.Blank_SAP.GCaMP.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.GCaMP.Rest.data));
scatter(xInds*2,gcampCorrData.SSP_SAP.GCaMP.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,gcampCorrData.SSP_SAP.GCaMP.Rest.meanData,gcampCorrData.SSP_SAP.GCaMP.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampCorrData.Blank_SAP.GCaMP.Alert.data));
scatter(xInds*4,gcampCorrData.Blank_SAP.GCaMP.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(4,gcampCorrData.Blank_SAP.GCaMP.Alert.meanData,gcampCorrData.Blank_SAP.GCaMP.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.GCaMP.Alert.data));
scatter(xInds*5,gcampCorrData.SSP_SAP.GCaMP.Alert.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampCorrData.SSP_SAP.GCaMP.Alert.meanData,gcampCorrData.SSP_SAP.GCaMP.Alert.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampCorrData.Blank_SAP.GCaMP.Asleep.data));
scatter(xInds*6,gcampCorrData.Blank_SAP.GCaMP.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e2 = errorbar(6,gcampCorrData.Blank_SAP.GCaMP.Asleep.meanData,gcampCorrData.Blank_SAP.GCaMP.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampCorrData.SSP_SAP.GCaMP.Asleep.data));
scatter(xInds*7,gcampCorrData.SSP_SAP.GCaMP.Asleep.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(7,gcampCorrData.SSP_SAP.GCaMP.Asleep.meanData,gcampCorrData.SSP_SAP.GCaMP.Asleep.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('Corr Coef')
xlim([0,8])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

%% save figure and stats
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figs & Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS6,[dirpath 'FigS6']);
    % statistical diary
    diaryFile = [dirpath 'FigS6_Stats.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on

    % [HbT] Pearson's correlation during each arousal-state
    disp('======================================================================================================================')
    disp('[HbT] Pearson''s correlation during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(ephysCorrData.Blank_SAP.HbT.Rest.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.HbT.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(ephysCorrData.SSP_SAP.HbT.Rest.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.HbT.Rest.stdData)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(EphysPearsonStats.HbT.Rest.p)]); disp(' ')
    disp(['Blank-SAP Alert: ' num2str(ephysCorrData.Blank_SAP.HbT.Alert.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.HbT.Alert.stdData)]); disp(' ')
    disp(['SSP-SAP Alert: ' num2str(ephysCorrData.SSP_SAP.HbT.Alert.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.HbT.Alert.stdData)]); disp(' ')
    disp(['Blank vs. SAP Alert ttest p = ' num2str(EphysPearsonStats.HbT.Alert.p)]); disp(' ')
    disp(['Blank-SAP Asleep: ' num2str(ephysCorrData.Blank_SAP.HbT.Asleep.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.HbT.Asleep.stdData)]); disp(' ')
    disp(['SSP-SAP Asleep: ' num2str(ephysCorrData.SSP_SAP.HbT.Asleep.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.HbT.Asleep.stdData)]); disp(' ')
    disp(['Blank vs. SAP Asleep ttest p = ' num2str(EphysPearsonStats.HbT.Asleep.p)]); disp(' ')

    % [gammaBandPower] Pearson's correlation during each arousal-state
    disp('======================================================================================================================')
    disp('[gammaBandPower] Pearson''s correlation during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Rest.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Rest.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Rest.stdData)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(EphysPearsonStats.gammaBandPower.Rest.p)]); disp(' ')
    disp(['Blank-SAP Alert: ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Alert.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Alert.stdData)]); disp(' ')
    disp(['SSP-SAP Alert: ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Alert.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Alert.stdData)]); disp(' ')
    disp(['Blank vs. SAP Alert ttest p = ' num2str(EphysPearsonStats.gammaBandPower.Alert.p)]); disp(' ')
    disp(['Blank-SAP Asleep: ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Asleep.meanData) ' +/- ' num2str(ephysCorrData.Blank_SAP.gammaBandPower.Asleep.stdData)]); disp(' ')
    disp(['SSP-SAP Asleep: ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Asleep.meanData) ' +/- ' num2str(ephysCorrData.SSP_SAP.gammaBandPower.Asleep.stdData)]); disp(' ')
    disp(['Blank vs. SAP Asleep ttest p = ' num2str(EphysPearsonStats.gammaBandPower.Asleep.p)]); disp(' ')

    % [HbT] Pearson's correlation during each arousal-state
    disp('======================================================================================================================')
    disp('[HbT] Pearson''s correlation during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampCorrData.Blank_SAP.HbT.Rest.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.HbT.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampCorrData.SSP_SAP.HbT.Rest.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.HbT.Rest.stdData)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPPearsonStats.HbT.Rest.p)]); disp(' ')
    disp(['Blank-SAP Alert: ' num2str(gcampCorrData.Blank_SAP.HbT.Alert.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.HbT.Alert.stdData)]); disp(' ')
    disp(['SSP-SAP Alert: ' num2str(gcampCorrData.SSP_SAP.HbT.Alert.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.HbT.Alert.stdData)]); disp(' ')
    disp(['Blank vs. SAP Alert ttest p = ' num2str(GCaMPPearsonStats.HbT.Alert.p)]); disp(' ')
    disp(['Blank-SAP Asleep: ' num2str(gcampCorrData.Blank_SAP.HbT.Asleep.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.HbT.Asleep.stdData)]); disp(' ')
    disp(['SSP-SAP Asleep: ' num2str(gcampCorrData.SSP_SAP.HbT.Asleep.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.HbT.Asleep.stdData)]); disp(' ')
    disp(['Blank vs. SAP Asleep ttest p = ' num2str(GCaMPPearsonStats.HbT.Asleep.p)]); disp(' ')

    % [GCaMP] Pearson's correlation during each arousal-state
    disp('======================================================================================================================')
    disp('[GCaMP] Pearson''s correlation during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampCorrData.Blank_SAP.GCaMP.Rest.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.GCaMP.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampCorrData.SSP_SAP.GCaMP.Rest.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.GCaMP.Rest.stdData)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPPearsonStats.GCaMP.Rest.p)]); disp(' ')
    disp(['Blank-SAP Alert: ' num2str(gcampCorrData.Blank_SAP.GCaMP.Alert.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.GCaMP.Alert.stdData)]); disp(' ')
    disp(['SSP-SAP Alert: ' num2str(gcampCorrData.SSP_SAP.GCaMP.Alert.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.GCaMP.Alert.stdData)]); disp(' ')
    disp(['Blank vs. SAP Alert ttest p = ' num2str(GCaMPPearsonStats.GCaMP.Alert.p)]); disp(' ')
    disp(['Blank-SAP Asleep: ' num2str(gcampCorrData.Blank_SAP.GCaMP.Asleep.meanData) ' +/- ' num2str(gcampCorrData.Blank_SAP.GCaMP.Asleep.stdData)]); disp(' ')
    disp(['SSP-SAP Asleep: ' num2str(gcampCorrData.SSP_SAP.GCaMP.Asleep.meanData) ' +/- ' num2str(gcampCorrData.SSP_SAP.GCaMP.Asleep.stdData)]); disp(' ')
    disp(['Blank vs. SAP Asleep ttest p = ' num2str(GCaMPPearsonStats.GCaMP.Asleep.p)]); disp(' ')

    diary off
end