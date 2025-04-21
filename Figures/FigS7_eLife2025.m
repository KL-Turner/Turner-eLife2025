function [] = FigS7_eLife2025(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% arousal-state hemodnyamics [GCaMP]
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','Whisk','Stim','NREM','REM'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
variables = {'avg','p2p','vari'};
fs = 10;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                gcampData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(gcampData.(group).(hemisphere).(dataType),behavior) == false
                        gcampData.(group).(hemisphere).(dataType).(behavior).avg = [];
                        gcampData.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        gcampData.(group).(hemisphere).(dataType).(behavior).vari = [];
                        gcampData.(group).(hemisphere).(dataType).(behavior).group = {};
                        gcampData.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    animalVar = [];
                    animalP2P = [];
                    for ff = 1:length(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData)
                        if strcmp(behavior,'Rest') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1}(2*fs:end);
                        elseif strcmp(behavior,'Stim') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1} - mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData),1);
                        else
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1};
                        end
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    try
                        gcampData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean));
                    catch
                        gcampData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).avg,mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean)));
                    end
                    gcampData.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    gcampData.(group).(hemisphere).(dataType).(behavior).vari = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    gcampData.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).group,group);
                    gcampData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    gcampData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(gcampData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    gcampData.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(gcampData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                end
            end
        end
    end
end
% statistics - ttest
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                % statistics - unpaired ttest
                [GCaMPStats.(hemisphere).(dataType).(behavior).(variable).h,GCaMPStats.(hemisphere).(dataType).(behavior).(variable).p] = ttest2(gcampData.Blank_SAP.(hemisphere).(dataType).(behavior).(variable),gcampData.SSP_SAP.(hemisphere).(dataType).(behavior).(variable));
            end
        end
    end
end


cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Transitions_Ephys';
load(resultsStruct);
cd(rootFolder)
groups = {'Blank_SAP','SSP_SAP'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};

%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Transitions_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(transitions)
            transition = transitions{1,cc};
            % pre-allocate necessary variable fields
            transitionData.(group).(transition).dummCheck = 1;
            if isfield(transitionData.(group).(transition),'HbT') == false
                transitionData.(group).(transition).HbT = [];
                transitionData.(group).(transition).difference = [];
                transitionData.(group).(transition).animalID = {};
                transitionData.(group).(transition).group = {};
            end
            transitionData.(group).(transition).HbT = cat(1,transitionData.(group).(transition).HbT,Results_Transitions_Ephys.(group).(animalID).RH.(transition).HbT);
            leadHbT = mean(Results_Transitions_Ephys.(group).(animalID).RH.(transition).HbT(1:30*20 + 1));
            lagHbT = mean(Results_Transitions_Ephys.(group).(animalID).RH.(transition).HbT(end - 30*20:end));
            transitionData.(group).(transition).difference = cat(1,transitionData.(group).(transition).difference,leadHbT - lagHbT);
            transitionData.(group).(transition).animalID = cat(1,transitionData.(group).(transition).animalID,animalID);
            transitionData.(group).(transition).group = cat(1,transitionData.(group).(transition).group,group);
        end
    end
end
% take average for each behavioral transition
for qq = 1:length(groups)
    group = groups{1,qq};
    for cc = 1:length(transitions)
        transition = transitions{1,cc};
        transitionData.(group).(transition).meanHbT = mean(transitionData.(group).(transition).HbT,1);
        transitionData.(group).(transition).stdHbT = std(transitionData.(group).(transition).HbT,0,1);
        transitionData.(group).(transition).meanDifference = mean(transitionData.(group).(transition).difference,1);
        transitionData.(group).(transition).stdDifference = std(transitionData.(group).(transition).difference,0,1);
    end
end

% GLME comparing peak correlation
for bb = 1:length(transitions)
    transition = transitions{1,bb};
    transitionStats.(transition).tableSize = cat(1,transitionData.Blank_SAP.(transition).difference,transitionData.SSP_SAP.(transition).difference);
    transitionStats.(transition).Table = table('Size',[size(transitionStats.(transition).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'AnimalID','Treatment','Difference'});
    transitionStats.(transition).Table.AnimalID = cat(1,transitionData.Blank_SAP.(transition).animalID,transitionData.SSP_SAP.(transition).animalID);
    transitionStats.(transition).Table.Treatment = cat(1,transitionData.Blank_SAP.(transition).group,transitionData.SSP_SAP.(transition).group);
    transitionStats.(transition).Table.Difference = cat(1,transitionData.Blank_SAP.(transition).difference,transitionData.SSP_SAP.(transition).difference);
    transitionStats.(transition).FitFormula = 'Difference ~ 1 + Treatment + (1|AnimalID)';
    transitionStats.(transition).Stats = fitglme(transitionStats.(transition).Table,transitionStats.(transition).FitFormula);
end
%% Figure S7
FigS7 = figure('Name','Fig. S7');

subplot(3,3,1);
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbT.NREM.avg));
hold on
scatter(xInds*4,gcampData.Blank_SAP.RH.HbT.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(4,gcampData.Blank_SAP.RH.HbT.NREM.mean_avg,gcampData.Blank_SAP.RH.HbT.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbT.NREM.avg));
scatter(xInds*5,gcampData.SSP_SAP.RH.HbT.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampData.SSP_SAP.RH.HbT.NREM.mean_avg,gcampData.SSP_SAP.RH.HbT.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbT.REM.avg));
scatter(xInds*7,gcampData.Blank_SAP.RH.HbT.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(7,gcampData.Blank_SAP.RH.HbT.REM.mean_avg,gcampData.Blank_SAP.RH.HbT.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbT.REM.avg));
scatter(xInds*8,gcampData.SSP_SAP.RH.HbT.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(8,gcampData.SSP_SAP.RH.HbT.REM.mean_avg,gcampData.SSP_SAP.RH.HbT.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('\Delta[HbT] (\muM)')
xlim([0,9])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(3,3,2);
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbO.NREM.avg));
hold on
scatter(xInds*4,gcampData.Blank_SAP.RH.HbO.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(4,gcampData.Blank_SAP.RH.HbO.NREM.mean_avg,gcampData.Blank_SAP.RH.HbO.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbO.NREM.avg));
scatter(xInds*5,gcampData.SSP_SAP.RH.HbO.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampData.SSP_SAP.RH.HbO.NREM.mean_avg,gcampData.SSP_SAP.RH.HbO.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbO.REM.avg));
scatter(xInds*7,gcampData.Blank_SAP.RH.HbO.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(7,gcampData.Blank_SAP.RH.HbO.REM.mean_avg,gcampData.Blank_SAP.RH.HbO.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbO.REM.avg));
scatter(xInds*8,gcampData.SSP_SAP.RH.HbO.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(8,gcampData.SSP_SAP.RH.HbO.REM.mean_avg,gcampData.SSP_SAP.RH.HbO.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('\Delta[HbO] (\muM)')
xlim([0,9])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(3,3,3);
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbR.NREM.avg));
hold on
scatter(xInds*4,gcampData.Blank_SAP.RH.HbR.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(4,gcampData.Blank_SAP.RH.HbR.NREM.mean_avg,gcampData.Blank_SAP.RH.HbR.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbR.NREM.avg));
scatter(xInds*5,gcampData.SSP_SAP.RH.HbR.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,gcampData.SSP_SAP.RH.HbR.NREM.mean_avg,gcampData.SSP_SAP.RH.HbR.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(gcampData.Blank_SAP.RH.HbR.REM.avg));
scatter(xInds*7,gcampData.Blank_SAP.RH.HbR.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(7,gcampData.Blank_SAP.RH.HbR.REM.mean_avg,gcampData.Blank_SAP.RH.HbR.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(gcampData.SSP_SAP.RH.HbR.REM.avg));
scatter(xInds*8,gcampData.SSP_SAP.RH.HbR.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(8,gcampData.SSP_SAP.RH.HbR.REM.mean_avg,gcampData.SSP_SAP.RH.HbR.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('\Delta[HbR] (\muM)')
xlim([0,9])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

T1 = -30 + (1/30):(1/30):30;
for aa = 1:length(transitions)
    transition = transitions{1,aa};
    subplot(3,3,3 + aa);
    p1 = plot(T1,transitionData.Blank_SAP.(transition).meanHbT,'-','color',colors('black'),'LineWidth',2);
    hold on
    plot(T1,transitionData.Blank_SAP.(transition).meanHbT + transitionData.Blank_SAP.(transition).stdHbT,'-','color',colors('black'),'LineWidth',0.5);
    plot(T1,transitionData.Blank_SAP.(transition).meanHbT - transitionData.Blank_SAP.(transition).stdHbT,'-','color',colors('black'),'LineWidth',0.5);
    p2 = plot(T1,transitionData.SSP_SAP.(transition).meanHbT,'-','color',colors('dark candy apple red'),'LineWidth',2);
    plot(T1,transitionData.SSP_SAP.(transition).meanHbT + transitionData.SSP_SAP.(transition).stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,transitionData.SSP_SAP.(transition).meanHbT - transitionData.SSP_SAP.(transition).stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    ylabel('\Delta[HbT] (\muM)')
    title(transition)
    legend([p1,p2],'Blank','SP')
    xlim([-30,30])
end

%% save figure and stats
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figs & Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS7,[dirpath 'FigS7']);
    % statistical diary
    diaryFile = [dirpath 'FigS7_Stats.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on

    % [HbT] during each arousal-state
    disp('======================================================================================================================')
    disp('[HbT] during each arousal-state (GCaMP), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampData.Blank_SAP.RH.HbT.Rest.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbT.Rest.std_avg)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampData.SSP_SAP.RH.HbT.Rest.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbT.Rest.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPStats.RH.HbT.Rest.avg.p)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(gcampData.Blank_SAP.RH.HbT.NREM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbT.NREM.std_avg)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(gcampData.SSP_SAP.RH.HbT.NREM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbT.NREM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(GCaMPStats.RH.HbT.NREM.avg.p)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(gcampData.Blank_SAP.RH.HbT.REM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbT.REM.std_avg)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(gcampData.SSP_SAP.RH.HbT.REM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbT.REM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(GCaMPStats.RH.HbT.REM.avg.p)]); disp(' ')

    % [HbO] during each arousal-state
    disp('======================================================================================================================')
    disp('[HbO] during each arousal-state (GCaMP), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampData.Blank_SAP.RH.HbO.Rest.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbO.Rest.std_avg)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampData.SSP_SAP.RH.HbO.Rest.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbO.Rest.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPStats.RH.HbO.Rest.avg.p)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(gcampData.Blank_SAP.RH.HbO.NREM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbO.NREM.std_avg)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(gcampData.SSP_SAP.RH.HbO.NREM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbO.NREM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(GCaMPStats.RH.HbO.NREM.avg.p)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(gcampData.Blank_SAP.RH.HbO.REM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbO.REM.std_avg)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(gcampData.SSP_SAP.RH.HbO.REM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbO.REM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(GCaMPStats.RH.HbO.REM.avg.p)]); disp(' ')

    % [HbR] during each arousal-state
    disp('======================================================================================================================')
    disp('[HbR] during each arousal-state (GCaMP), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(gcampData.Blank_SAP.RH.HbR.Rest.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbR.Rest.std_avg)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(gcampData.SSP_SAP.RH.HbR.Rest.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbR.Rest.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(GCaMPStats.RH.HbR.Rest.avg.p)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(gcampData.Blank_SAP.RH.HbR.NREM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbR.NREM.std_avg)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(gcampData.SSP_SAP.RH.HbR.NREM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbR.NREM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(GCaMPStats.RH.HbR.NREM.avg.p)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(gcampData.Blank_SAP.RH.HbR.REM.mean_avg) ' +/- ' num2str(gcampData.Blank_SAP.RH.HbR.REM.std_avg)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(gcampData.SSP_SAP.RH.HbR.REM.mean_avg) ' +/- ' num2str(gcampData.SSP_SAP.RH.HbR.REM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(GCaMPStats.RH.HbR.REM.avg.p)]); disp(' ')

    % Awake to NREM
    disp('======================================================================================================================')
    disp('Awake to NREM transition, n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(transitionData.Blank_SAP.AWAKEtoNREM.meanDifference)) ' +/- ' num2str(transitionData.Blank_SAP.AWAKEtoNREM.stdDifference)]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(transitionData.SSP_SAP.AWAKEtoNREM.meanDifference)) ' +/- ' num2str(transitionData.SSP_SAP.AWAKEtoNREM.stdDifference)]); disp(' ')
    disp(transitionStats.AWAKEtoNREM.Stats)

    % NREM to Awake
    disp('======================================================================================================================')
    disp('NREM to Awake transition, n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(transitionData.Blank_SAP.NREMtoAWAKE.meanDifference)) ' +/- ' num2str(transitionData.Blank_SAP.NREMtoAWAKE.stdDifference)]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(transitionData.SSP_SAP.NREMtoAWAKE.meanDifference)) ' +/- ' num2str(transitionData.SSP_SAP.NREMtoAWAKE.stdDifference)]); disp(' ')
    disp(transitionStats.NREMtoAWAKE.Stats)

    % NREM to REM
    disp('======================================================================================================================')
    disp('NREM to REM transition, n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(transitionData.Blank_SAP.NREMtoREM.meanDifference)) ' +/- ' num2str(transitionData.Blank_SAP.NREMtoREM.stdDifference)]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(transitionData.SSP_SAP.NREMtoREM.meanDifference)) ' +/- ' num2str(transitionData.SSP_SAP.NREMtoREM.stdDifference)]); disp(' ')
    disp(transitionStats.NREMtoREM.Stats)

    % REM to Awake
    disp('======================================================================================================================')
    disp('REM to Awake transition, n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(transitionData.Blank_SAP.REMtoAWAKE.meanDifference)) ' +/- ' num2str(transitionData.Blank_SAP.REMtoAWAKE.stdDifference)]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(transitionData.SSP_SAP.REMtoAWAKE.meanDifference)) ' +/- ' num2str(transitionData.SSP_SAP.REMtoAWAKE.stdDifference)]); disp(' ')
    disp(transitionStats.REMtoAWAKE.Stats)

end
