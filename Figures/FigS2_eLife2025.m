function [] = FigS2_eLife2025(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% open field tracking example
path = [rootFolder delim 'Results_Hossain'];
cd(path)
load('T267_Run1.mat');
load('OFTTable.mat')
animalID = OFTData.MouseID;
frameX = 1;
% mouse body parts
mid = OFTData.TrackData(:,7:8);
tailTrunk = OFTData.TrackData(:,11:12);
% outer
OTR = OFTData.TrackData(:,17:18);
OBR = OFTData.TrackData(:,19:20);
OTL = OFTData.TrackData(:,21:22);
OBL = OFTData.TrackData(:,23:24);
Zodd = 1:2:9;
Zeven = 2:2:10;
outerpoints = [OTR,OTL,OBL,OBR,OTR];

%% open field distance traveled
[G,~] = grp2idx(AnimalTable_OFT.DrugGroup);
Cidx = find(G==2);
Bidx = find(G==1);
Sidx = find(G==3);
Distance_Control = AnimalTable_OFT.total_distance_5min(Cidx,:);
Distance_Blank = AnimalTable_OFT.total_distance_5min(Bidx,:);
Distance_SSP = AnimalTable_OFT.total_distance_5min(Sidx,:);
meanDistance_Control = mean(Distance_Control);
meanDistance_Blank = mean(Distance_Blank);
meanDistance_SSP = mean(Distance_SSP);
stdDistance_Control = std(Distance_Control);
stdDistance_Blank = std(Distance_Blank);
stdDistance_SSP = std(Distance_SSP);
naiveGroup = cell(length(Distance_Control),1);
naiveGroup(:) = {'Naive'};
blankGroup = cell(length(Distance_Blank),1);
blankGroup(:) = {'Blank_SAP'};
sspGroup = cell(length(Distance_SSP),1);
sspGroup(:) = {'SSP_SAP'};
% Naive vs blank RH
DistanceStats.NaiveBlank.tableSize = cat(1,Distance_Control,Distance_Blank);
DistanceStats.NaiveBlank.Table = table('Size',[size(DistanceStats.NaiveBlank.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Distance'});
DistanceStats.NaiveBlank.Table.Group = cat(1,naiveGroup,blankGroup);
DistanceStats.NaiveBlank.Table.Distance = cat(1,Distance_Control,Distance_Blank);
DistanceStats.NaiveBlank.FitFormula = 'Distance ~ 1 + Group';
DistanceStats.NaiveBlank.Stats = fitglme(DistanceStats.NaiveBlank.Table,DistanceStats.NaiveBlank.FitFormula);
% Blank vs SSP RH
DistanceStats.BlankSSP.tableSize = cat(1,Distance_Blank,Distance_SSP);
DistanceStats.BlankSSP.Table = table('Size',[size(DistanceStats.BlankSSP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Distance'});
DistanceStats.BlankSSP.Table.Group = cat(1,blankGroup,sspGroup);
DistanceStats.BlankSSP.Table.Distance = cat(1,Distance_Blank,Distance_SSP);
DistanceStats.BlankSSP.FitFormula = 'Distance ~ 1 + Group';
DistanceStats.BlankSSP.Stats = fitglme(DistanceStats.BlankSSP.Table,DistanceStats.BlankSSP.FitFormula);

%% open field time spent in center
[G,~] = grp2idx(AnimalTable_OFT.DrugGroup);
Cidx = find(G==2);
Bidx = find(G==1);
Sidx = find(G==3);
CenterTime_Control = AnimalTable_OFT.center_percentage_5min(Cidx,:);
CenterTime_Blank = AnimalTable_OFT.center_percentage_5min(Bidx,:);
CenterTime_SSP = AnimalTable_OFT.center_percentage_5min(Sidx,:);
meanCenterTime_Control = mean(CenterTime_Control);
meanCenterTime_Blank = mean(CenterTime_Blank);
meanCenterTime_SSP = mean(CenterTime_SSP);
stdCenterTime_Control = std(CenterTime_Control);
stdCenterTime_Blank = std(CenterTime_Blank);
stdCenterTime_SSP = std(CenterTime_SSP);
% Naive vs blank RH
CenterStats.NaiveBlank.tableSize = cat(1,CenterTime_Control,CenterTime_Blank);
CenterStats.NaiveBlank.Table = table('Size',[size(CenterStats.NaiveBlank.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Percent'});
CenterStats.NaiveBlank.Table.Group = cat(1,naiveGroup,blankGroup);
CenterStats.NaiveBlank.Table.Percent = cat(1,CenterTime_Control,CenterTime_Blank);
CenterStats.NaiveBlank.FitFormula = 'Percent ~ 1 + Group';
CenterStats.NaiveBlank.Stats = fitglme(CenterStats.NaiveBlank.Table,CenterStats.NaiveBlank.FitFormula);
% Blank vs SSP RH
CenterStats.BlankSSP.tableSize = cat(1,CenterTime_Blank,CenterTime_SSP);
CenterStats.BlankSSP.Table = table('Size',[size(CenterStats.BlankSSP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Percent'});
CenterStats.BlankSSP.Table.Group = cat(1,blankGroup,sspGroup);
CenterStats.BlankSSP.Table.Percent = cat(1,CenterTime_Blank,CenterTime_SSP);
CenterStats.BlankSSP.FitFormula = 'Percent ~ 1 + Group';
CenterStats.BlankSSP.Stats = fitglme(CenterStats.BlankSSP.Table,CenterStats.BlankSSP.FitFormula);

%% arousal state probability
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_ArousalStateProb_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'awakePercent','nremPercent','remPercent'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_ArousalStateProb_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        ephysStateData.(group).dummCheck = 1;
        ephysStateData.(group).group = {};
        ephysStateData.(group).animalID = {};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            if isfield(ephysStateData.(group),dataType) == false
                ephysStateData.(group).(dataType) = [];
            end
            ephysStateData.(group).(dataType) = cat(1,ephysStateData.(group).(dataType),Results_ArousalStateProb_Ephys.(group).(animalID).(dataType));
        end
        ephysStateData.(group).group = cat(1,ephysStateData.(group).group,group);
        ephysStateData.(group).animalID = cat(1,ephysStateData.(group).animalID,animalID);
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        ephysStateData.(group).(['mean_' dataType]) = mean(ephysStateData.(group).(dataType),1);
        ephysStateData.(group).(['std_' dataType]) = std(ephysStateData.(group).(dataType),0,1);
    end
    ephysStateData.(group).meanPercs = cat(1,ephysStateData.(group).mean_awakePercent,ephysStateData.(group).mean_nremPercent,ephysStateData.(group).mean_remPercent);
end
% statistics - unpaired ttest
% [awakeStats1.h,awakeStats1.p] = vartest(ephysStateData.Naive.awakePercent,ephysStateData.Blank_SAP.awakePercent);
% [awakeStats2.h,awakeStats2.p] = vartest(ephysStateData.Blank_SAP.awakePercent,ephysStateData.SSP_SAP.awakePercent);
% [nremStats1.h,nremStats1.p] = vartest(ephysStateData.Naive.nremPercent,ephysStateData.Blank_SAP.nremPercent);
% [nremStats2.h,nremStats2.p] = vartest(ephysStateData.Blank_SAP.nremPercent,ephysStateData.SSP_SAP.nremPercent);
% [remStats1.h,remStats1.p] = vartest(ephysStateData.Naive.remPercent,ephysStateData.Blank_SAP.remPercent);
% [remStats2.h,remStats2.p] = vartest(ephysStateData.Blank_SAP.remPercent,ephysStateData.SSP_SAP.remPercent);
% 
% statistics - unpaired ttest
[awakeStats1.h,awakeStats1.p] = ttest2(ephysStateData.Naive.awakePercent,ephysStateData.Blank_SAP.awakePercent);
[awakeStats2.h,awakeStats2.p] = ttest2(ephysStateData.Blank_SAP.awakePercent,ephysStateData.SSP_SAP.awakePercent);
[nremStats1.h,nremStats1.p] = ttest2(ephysStateData.Naive.nremPercent,ephysStateData.Blank_SAP.nremPercent);
[nremStats2.h,nremStats2.p] = ttest2(ephysStateData.Blank_SAP.nremPercent,ephysStateData.SSP_SAP.nremPercent);
[remStats1.h,remStats1.p] = ttest2(ephysStateData.Naive.remPercent,ephysStateData.Blank_SAP.remPercent);
[remStats2.h,remStats2.p] = ttest2(ephysStateData.Blank_SAP.remPercent,ephysStateData.SSP_SAP.remPercent);

% Blank.Awake = ephysStateData.Blank_SAP.awakePercent;
% Blank.NREM = ephysStateData.Blank_SAP.nremPercent;
% Blank.REM = ephysStateData.Blank_SAP.remPercent;
% 
% SP.Awake = ephysStateData.SSP_SAP.awakePercent;
% SP.NREM = ephysStateData.SSP_SAP.nremPercent;
% SP.REM = ephysStateData.SSP_SAP.remPercent;
% 
% Uninj.Awake = ephysStateData.Naive.awakePercent;
% Uninj.NREM = ephysStateData.Naive.nremPercent;
% Uninj.REM = ephysStateData.Naive.remPercent;
% 
% statePercentages = table(Blank.Awake,SP.Awake,Uninj.Awake);
% 
% [p,h,stats] = ttest2(ephysStateData.Blank_SAP.awakePercent,ephysStateData.SSP_SAP.awakePercent)
% [p,h,stats] = ranksum(ephysStateData.Blank_SAP.awakePercent,ephysStateData.SSP_SAP.awakePercent)

%% whisking behavior
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_WhiskBehav_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'whiskDurationSec','whiskDurationPerc'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_WhiskBehav_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        whiskingData.(group).dummCheck = 1;
        whiskingData.(group).group = {};
        whiskingData.(group).animalID = {};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            if isfield(whiskingData.(group),dataType) == false
                whiskingData.(group).(dataType) = [];
            end
        end
        % concatenate data across animals
        whiskingData.(group).whiskDurationSec = cat(1,whiskingData.(group).whiskDurationSec,Results_WhiskBehav_Ephys.(group).(animalID).whiskDurationSec/60);
        whiskingData.(group).whiskDurationPerc = cat(1,whiskingData.(group).whiskDurationPerc,Results_WhiskBehav_Ephys.(group).(animalID).whiskDurationPerc);
        whiskingData.(group).group = cat(1,whiskingData.(group).group,group);
        whiskingData.(group).animalID = cat(1,whiskingData.(group).animalID,animalID);
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        whiskingData.(group).(['mean_' dataType]) = mean(whiskingData.(group).(dataType),1);
        whiskingData.(group).(['std_' dataType]) = std(whiskingData.(group).(dataType),0,1);
    end
end
% statistics - unpaired ttest
[whiskPercStats1.h,whiskPercStats1.p] = ttest2(whiskingData.Naive.whiskDurationPerc,whiskingData.Blank_SAP.whiskDurationPerc);
[whiskPercStats2.h,whiskPercStats2.p] = ttest2(whiskingData.Blank_SAP.whiskDurationPerc,whiskingData.SSP_SAP.whiskDurationPerc);

%% pupil size
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PupilArea_Ephys';
load(resultsStruct);
cd(rootFolder)
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
behaviors = {'Rest','Whisk','Stim','NREM','REM'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PupilArea_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            pupilSizeData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(pupilSizeData.(group).(dataType),behavior) == false
                    pupilSizeData.(group).(dataType).(behavior).data = [];
                    pupilSizeData.(group).(dataType).(behavior).group = {};
                    pupilSizeData.(group).(dataType).(behavior).animalID = {};
                end
                if isempty(Results_PupilArea_Ephys.(group).(animalID).(dataType).(behavior).eventMeans) == false
                    pupilSizeData.(group).(dataType).(behavior).data = cat(1,pupilSizeData.(group).(dataType).(behavior).data,mean(Results_PupilArea_Ephys.(group).(animalID).(dataType).(behavior).eventMeans));
                    pupilSizeData.(group).(dataType).(behavior).group = cat(1,pupilSizeData.(group).(dataType).(behavior).group,group);
                    pupilSizeData.(group).(dataType).(behavior).animalID = cat(1,pupilSizeData.(group).(dataType).(behavior).animalID,animalID);
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
            pupilSizeData.(group).(dataType).(behavior).meanData = mean(pupilSizeData.(group).(dataType).(behavior).data,1);
            pupilSizeData.(group).(dataType).(behavior).stdData = std(pupilSizeData.(group).(dataType).(behavior).data,0,1);
        end
    end
end
% statistics - unpaired ttest
[pupilSizeStats.Rest1.h,pupilSizeStats.Rest1.p] = ttest2(pupilSizeData.Naive.zDiameter.Rest.data,pupilSizeData.Blank_SAP.zDiameter.Rest.data);
[pupilSizeStats.Rest2.h,pupilSizeStats.Rest2.p] = ttest2(pupilSizeData.Blank_SAP.zDiameter.Rest.data,pupilSizeData.SSP_SAP.zDiameter.Rest.data);
[pupilSizeStats.NREM1.h,pupilSizeStats.NREM1.p] = ttest2(pupilSizeData.Naive.zDiameter.NREM.data,pupilSizeData.Blank_SAP.zDiameter.NREM.data);
[pupilSizeStats.NREM2.h,pupilSizeStats.NREM2.p] = ttest2(pupilSizeData.Blank_SAP.zDiameter.NREM.data,pupilSizeData.SSP_SAP.zDiameter.NREM.data);
[pupilSizeStats.REM1.h,pupilSizeStats.REM1.p] = ttest2(pupilSizeData.Naive.zDiameter.REM.data,pupilSizeData.Blank_SAP.zDiameter.REM.data);
[pupilSizeStats.REM2.h,pupilSizeStats.REM2.p] = ttest2(pupilSizeData.Blank_SAP.zDiameter.REM.data,pupilSizeData.SSP_SAP.zDiameter.REM.data);

%% stimulus evoked pupil size
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PupilEvoked_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
solenoids = {'LPadSol','RPadSol','AudSol'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PupilEvoked_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(solenoids)
                solenoid = solenoids{1,dd};
                pupilEvokedData.(group).(dataType).(solenoid).dummCheck = 1;
                if isfield(pupilEvokedData.(group).(dataType).(solenoid),'data') == false
                    pupilEvokedData.(group).(dataType).(solenoid).data = [];
                    pupilEvokedData.(group).(dataType).(solenoid).peak = [];
                    pupilEvokedData.(group).(dataType).(solenoid).group = {};
                    pupilEvokedData.(group).(dataType).(solenoid).animalID = {};
                end
                pupilEvokedData.(group).(dataType).(solenoid).data = cat(1,pupilEvokedData.(group).(dataType).(solenoid).data,Results_PupilEvoked_Ephys.(group).(animalID).(dataType).Stim.(solenoid).mean);
                pupilEvokedData.(group).(dataType).(solenoid).peak = cat(1,pupilEvokedData.(group).(dataType).(solenoid).peak,max(Results_PupilEvoked_Ephys.(group).(animalID).(dataType).Stim.(solenoid).mean));
                pupilEvokedData.(group).(dataType).(solenoid).group = cat(1,pupilEvokedData.(group).(dataType).(solenoid).group,group);
                pupilEvokedData.(group).(dataType).(solenoid).animalID = cat(1,pupilEvokedData.(group).(dataType).(solenoid).animalID,animalID);
            end
        end
    end
end
% pair stimulation with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            [comparison] = FindSolenoidComparison('Both',solenoid);
            pupilEvokedData.(group).(dataType).(comparison).peak = pupilEvokedData.(group).(dataType).(solenoid).peak;
            pupilEvokedData.(group).(dataType).(comparison).meanPeak = mean(pupilEvokedData.(group).(dataType).(solenoid).peak);
            pupilEvokedData.(group).(dataType).(comparison).stdPeak = std(pupilEvokedData.(group).(dataType).(solenoid).peak,0,1);
            pupilEvokedData.(group).(dataType).(comparison).meanData = mean(pupilEvokedData.(group).(dataType).(solenoid).data,1);
            pupilEvokedData.(group).(dataType).(comparison).stdErrData = std(pupilEvokedData.(group).(dataType).(solenoid).data,0,1)./sqrt(size(pupilEvokedData.(group).(dataType).(solenoid).data,1));
        end
    end
end
% statistics - unpaired ttest
[pupilEvokedStats1.h,pupilEvokedStats1.p] = ttest2(pupilEvokedData.Naive.zDiameter.contra.peak,pupilEvokedData.Blank_SAP.zDiameter.contra.peak);
[pupilEvokedStats2.h,pupilEvokedStats2.p] = ttest2(pupilEvokedData.Blank_SAP.zDiameter.contra.peak,pupilEvokedData.SSP_SAP.zDiameter.contra.peak);

%% blinking
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PupilBlinkInterval_Ephys';
load(resultsStruct);
cd(rootFolder)
groups = {'Naive','SSP_SAP','Blank_SAP'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    blinkData.(group).data = [];
    animalIDs = fieldnames(Results_PupilBlinkInterval_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        % concatenate data across animals
        blinkData.(group).data = cat(1,blinkData.(group).data,mean(Results_PupilBlinkInterval_Ephys.(group).(animalID).interBlinkInterval));
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    blinkData.(group).meanData = mean(blinkData.(group).data,1);
    blinkData.(group).stdData = std(blinkData.(group).data,0,1);
end
% statistics - unpaired ttest
[BlinkingStats1.h,BlinkingStats1.p] = ttest2(blinkData.Naive.data,blinkData.Blank_SAP.data);
[BlinkingStats2.h,BlinkingStats2.p] = ttest2(blinkData.Blank_SAP.data,blinkData.SSP_SAP.data);

%% Figure S2
FigS2 = figure('Name','Figure S2');

subplot(3,3,1)
plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'*-r')
hold on;
A(1) =  mid(1,1);
A(2) =  mid(1,2);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
p1 = plot(A(1),A(2),'ko','MarkerSize',10);
for i=1:1:1200
    B(1) =  mid(i,1);
    B(2) =  mid(i,2);
    A(1) =  tailTrunk(i,1);
    A(2) =  tailTrunk(i,2);
    hold on;
    vectarrow(A,B)
end
hold on
A(1) =  mid(1201,1);
A(2) =  mid(1201,2);
p2 = plot(A(1),A(2),'g*','MarkerSize',10);
axis equal
set(gca,'xticklabel',{[]},'yticklabel',{[]})
legend([p1,p2],'start','end')
ylabel('30 cm wide')
xlabel('60 cm long')
axis([-50,1000,-50,600])
xticks([])
yticks([])
title([animalID 'Vector Plot'])

subplot(3,3,2)
cplot = bar(1,meanDistance_Control,0.1);
cplot.FaceColor = colors('sapphire');
hold on
plot(1,Distance_Control,'ob','MarkerSize',8,'MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','k');
errorbar(1,meanDistance_Control,stdDistance_Control,'-k','CapSize',18,'LineWidth',3);
kplot= bar(1.25,meanDistance_Blank,0.1);
kplot.FaceColor = colors('north texas green');
plot(1.25,Distance_Blank,'go','MarkerSize',8,'MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','k');
errorbar(1.25,meanDistance_Blank,stdDistance_Blank,'-k','CapSize',18,'LineWidth',3);
nplot = bar(1.5,meanDistance_SSP,0.1);
nplot.FaceColor = colors('electric purple');
plot(1.5,Distance_SSP,'or','MarkerSize',8,'MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','k');
errorbar(1.5,meanDistance_SSP,stdDistance_SSP,'-k','CapSize',18,'LineWidth',3);
ylabel('Distance Travelled (cm)')
xlim([0.75 1.75])
xticks([1 1.25 1.5]);
xticklabels({'Naive','Blank-SAP','SSP-SAP'})
set(gca,'box','off')
axis square

subplot(3,3,3)
cplot = bar(1,meanCenterTime_Control,0.1);
cplot.FaceColor = colors('sapphire');
hold on
plot(1,CenterTime_Control,'ob','MarkerSize',8,'MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','k');
errorbar(1,meanCenterTime_Control,stdCenterTime_Control,'-k','CapSize',18,'LineWidth',3);
kplot= bar(1.25,meanCenterTime_Blank,0.1);
kplot.FaceColor = colors('north texas green');
plot(1.25,CenterTime_Blank,'go','MarkerSize',8,'MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','k');
errorbar(1.25,meanCenterTime_Blank,stdCenterTime_Blank,'-k','CapSize',18,'LineWidth',3);
nplot = bar(1.5,meanCenterTime_SSP,0.1);
nplot.FaceColor = colors('electric purple');
plot(1.5,CenterTime_SSP,'or','MarkerSize',8,'MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','k');
errorbar(1.5,meanCenterTime_SSP,stdCenterTime_SSP,'-k','CapSize',18,'LineWidth',3);
ylabel('Percentage Center Time')
xlim([0.75 1.75])
xticks([1 1.25 1.5]);
xticklabels({'Naive','Blank-SAP','SSP-SAP'})
set(gca,'box','off')
axis square

subplot(3,3,4)
p1 = pie(ephysStateData.Naive.meanPercs);
pText = findobj(p1,'Type','text');
percentValues = get(pText,'String');
txt = {'Awake: ';'NREM: ';'REM: '};
combinedtxt = strcat(txt,percentValues);
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title('Naive')

subplot(3,3,5)
p1 = pie(ephysStateData.Blank_SAP.meanPercs);
pText = findobj(p1,'Type','text');
percentValues = get(pText,'String');
txt = {'Awake: ';'NREM: ';'REM: '};
combinedtxt = strcat(txt,percentValues);
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title('Blank_SAP')

subplot(3,3,6)
p1 = pie(ephysStateData.SSP_SAP.meanPercs);
pText = findobj(p1,'Type','text');
percentValues = get(pText,'String');
txt = {'Awake: ';'NREM: ';'REM: '};
combinedtxt = strcat(txt,percentValues);
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title('SSP_SAP')

subplot(3,4,9)
s1 = scatter(ones(1,length(whiskingData.Naive.whiskDurationPerc))*1,whiskingData.Naive.whiskDurationPerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0.25);
hold on;
e1 = errorbar(1,whiskingData.Naive.mean_whiskDurationPerc,whiskingData.Naive.std_whiskDurationPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(whiskingData.Blank_SAP.whiskDurationPerc))*2,whiskingData.Blank_SAP.whiskDurationPerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,whiskingData.Blank_SAP.mean_whiskDurationPerc,whiskingData.Blank_SAP.std_whiskDurationPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(whiskingData.SSP_SAP.whiskDurationPerc))*3,whiskingData.SSP_SAP.whiskDurationPerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(3,whiskingData.SSP_SAP.mean_whiskDurationPerc,whiskingData.SSP_SAP.std_whiskDurationPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('Time spent whisking (%)')
legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
axis square
xlim([0,4])

subplot(3,4,10)
timeVector = (0:12*30)/30 - 2;
plot(timeVector,pupilEvokedData.Naive.zDiameter.contra.meanData,'color',colors('sapphire'),'LineWidth',2);
hold on;
plot(timeVector,pupilEvokedData.Naive.zDiameter.contra.meanData + pupilEvokedData.Naive.zDiameter.contra.stdErrData,'color',colors('sapphire'),'LineWidth',0.25)
plot(timeVector,pupilEvokedData.Naive.zDiameter.contra.meanData - pupilEvokedData.Naive.zDiameter.contra.stdErrData,'color',colors('sapphire'),'LineWidth',0.25)
plot(timeVector,pupilEvokedData.Blank_SAP.zDiameter.contra.meanData,'color',colors('north texas green'),'LineWidth',2);
plot(timeVector,pupilEvokedData.Blank_SAP.zDiameter.contra.meanData + pupilEvokedData.Blank_SAP.zDiameter.contra.stdErrData,'color',colors('north texas green'),'LineWidth',0.25)
plot(timeVector,pupilEvokedData.Blank_SAP.zDiameter.contra.meanData - pupilEvokedData.Blank_SAP.zDiameter.contra.stdErrData,'color',colors('north texas green'),'LineWidth',0.25)
plot(timeVector,pupilEvokedData.SSP_SAP.zDiameter.contra.meanData,'color',colors('electric purple'),'LineWidth',2);
plot(timeVector,pupilEvokedData.SSP_SAP.zDiameter.contra.meanData + pupilEvokedData.SSP_SAP.zDiameter.contra.stdErrData,'color',colors('electric purple'),'LineWidth',0.25)
plot(timeVector,pupilEvokedData.SSP_SAP.zDiameter.contra.meanData - pupilEvokedData.SSP_SAP.zDiameter.contra.stdErrData,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\Delta z-units')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
axis square

subplot(3,4,11)
scatter(ones(1,length(pupilSizeData.Naive.zDiameter.Rest.data))*1,pupilSizeData.Naive.zDiameter.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0);
hold on;
e1 = errorbar(1,pupilSizeData.Naive.zDiameter.Rest.meanData,pupilSizeData.Naive.zDiameter.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(pupilSizeData.Blank_SAP.zDiameter.Rest.data))*2,pupilSizeData.Blank_SAP.zDiameter.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0);
e2 = errorbar(2,pupilSizeData.Blank_SAP.zDiameter.Rest.meanData,pupilSizeData.Blank_SAP.zDiameter.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(pupilSizeData.SSP_SAP.zDiameter.Rest.data))*3,pupilSizeData.SSP_SAP.zDiameter.Rest.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0);
e3 = errorbar(3,pupilSizeData.SSP_SAP.zDiameter.Rest.meanData,pupilSizeData.SSP_SAP.zDiameter.Rest.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

scatter(ones(1,length(pupilSizeData.Naive.zDiameter.NREM.data))*4,pupilSizeData.Naive.zDiameter.NREM.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0);
hold on;
e1 = errorbar(4,pupilSizeData.Naive.zDiameter.NREM.meanData,pupilSizeData.Naive.zDiameter.NREM.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(pupilSizeData.Blank_SAP.zDiameter.NREM.data))*5,pupilSizeData.Blank_SAP.zDiameter.NREM.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0);
e2 = errorbar(5,pupilSizeData.Blank_SAP.zDiameter.NREM.meanData,pupilSizeData.Blank_SAP.zDiameter.NREM.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(pupilSizeData.SSP_SAP.zDiameter.NREM.data))*6,pupilSizeData.SSP_SAP.zDiameter.NREM.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0);
e3 = errorbar(6,pupilSizeData.SSP_SAP.zDiameter.NREM.meanData,pupilSizeData.SSP_SAP.zDiameter.NREM.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

scatter(ones(1,length(pupilSizeData.Naive.zDiameter.REM.data))*7,pupilSizeData.Naive.zDiameter.REM.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0);
hold on;
e1 = errorbar(7,pupilSizeData.Naive.zDiameter.REM.meanData,pupilSizeData.Naive.zDiameter.REM.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(pupilSizeData.Blank_SAP.zDiameter.REM.data))*8,pupilSizeData.Blank_SAP.zDiameter.REM.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0);
e2 = errorbar(8,pupilSizeData.Blank_SAP.zDiameter.REM.meanData,pupilSizeData.Blank_SAP.zDiameter.REM.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(pupilSizeData.SSP_SAP.zDiameter.REM.data))*9,pupilSizeData.SSP_SAP.zDiameter.REM.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0);
e3 = errorbar(9,pupilSizeData.SSP_SAP.zDiameter.REM.meanData,pupilSizeData.SSP_SAP.zDiameter.REM.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('z-units')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,10])
set(gca,'box','off')

subplot(3,4,12)
xInds = ones(1,length(blinkData.Naive.data));
scatter(xInds*1,blinkData.Naive.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0.25);
hold on;
e1 = errorbar(1,blinkData.Naive.meanData,blinkData.Naive.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(blinkData.Blank_SAP.data));
scatter(xInds*2,blinkData.Blank_SAP.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,blinkData.Blank_SAP.meanData,blinkData.Blank_SAP.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(blinkData.SSP_SAP.data));
scatter(xInds*3,blinkData.SSP_SAP.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(3,blinkData.SSP_SAP.meanData,blinkData.SSP_SAP.stdData,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('IBI (s)')
xlim([0,4])
set(gca,'box','off')
set(gca,'xtick',[])
axis square

%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figs & Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS2,[dirpath 'FigS2']);
    % statistical diary
    diaryFile = [dirpath 'FigS2_Stats.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on

    % open field distance traveled
    disp('======================================================================================================================')
    disp('Total distance traveled, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Naive ' num2str(meanDistance_Control) ' +/- ' num2str(stdDistance_Control)]); disp(' ')
    disp(['Blank-SAP ' num2str(meanDistance_Blank) ' +/- ' num2str(stdDistance_Blank)]); disp(' ')
    disp(['SSP-SAP ' num2str(meanDistance_SSP) ' +/- ' num2str(stdDistance_SSP)]); disp(' ')
    disp('Naive vs. Blank GLME'); disp(' ')
    disp(DistanceStats.NaiveBlank.Stats)
    disp('Blank vs. SSP GLME'); disp(' ')
    disp(DistanceStats.BlankSSP.Stats)

    % open field % time spent in center
    disp('======================================================================================================================')
    disp('Center Time, n = 12-23 mice per group, mean +/- StD'); disp(' ')
    disp(['Naive: ' num2str(meanCenterTime_Control) ' +/- ' num2str(stdCenterTime_Control)]); disp(' ')
    disp(['Blank-SAP: ' num2str(meanCenterTime_Blank) ' +/- ' num2str(stdCenterTime_Blank)]); disp(' ')
    disp(['SSP-SAP: ' num2str(meanCenterTime_SSP) ' +/- ' num2str(stdCenterTime_SSP)]); disp(' ')
    disp('Naive vs. Blank GLME'); disp(' ')
    disp(CenterStats.NaiveBlank.Stats)
    disp('Blank vs. SSP GLME'); disp(' ')
    disp(CenterStats.BlankSSP.Stats)

    % arousal-state percentages
    disp('======================================================================================================================')
    disp('Arousal-state %s, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Naive Awake: ' num2str(ephysStateData.Naive.mean_awakePercent) ' +/- ' num2str(ephysStateData.Naive.std_awakePercent)]); disp(' ')
    disp(['Naive NREM: ' num2str(ephysStateData.Naive.mean_nremPercent) ' +/- ' num2str(ephysStateData.Naive.std_nremPercent)]); disp(' ')
    disp(['Naive REM: ' num2str(ephysStateData.Naive.mean_remPercent) ' +/- ' num2str(ephysStateData.Naive.std_remPercent)]); disp(' ')
    disp(['Blank-SAP Awake: ' num2str(ephysStateData.Blank_SAP.mean_awakePercent) ' +/- ' num2str(ephysStateData.Blank_SAP.std_awakePercent)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(ephysStateData.Blank_SAP.mean_nremPercent) ' +/- ' num2str(ephysStateData.Blank_SAP.std_nremPercent)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(ephysStateData.Blank_SAP.mean_remPercent) ' +/- ' num2str(ephysStateData.Blank_SAP.std_remPercent)]); disp(' ')
    disp(['SSP-SAP Awake: ' num2str(ephysStateData.SSP_SAP.mean_awakePercent) ' +/- ' num2str(ephysStateData.SSP_SAP.std_awakePercent)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(ephysStateData.SSP_SAP.mean_nremPercent) ' +/- ' num2str(ephysStateData.SSP_SAP.std_nremPercent)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(ephysStateData.SSP_SAP.mean_remPercent) ' +/- ' num2str(ephysStateData.SSP_SAP.std_remPercent)]); disp(' ')
    disp(['Naive vs. Blank Awake ttest p = ' num2str(awakeStats1.p)]); disp(' ')
    disp(['Blank vs. SAP Awake ttest p = ' num2str(awakeStats2.p)]); disp(' ')
    disp(['Naive vs. Blank NREM ttest p = ' num2str(nremStats1.p)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(nremStats2.p)]); disp(' ')
    disp(['Naive vs. Blank REM ttest p = ' num2str(remStats1.p)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(remStats2.p)]); disp(' ')

    % time spent whisking
    disp('======================================================================================================================')
    disp('Time spent whisking, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Naive: ' num2str(whiskingData.Naive.mean_whiskDurationPerc) ' +/- ' num2str(whiskingData.Naive.std_whiskDurationPerc)]); disp(' ')
    disp(['Blank-SAP: ' num2str(whiskingData.Blank_SAP.mean_whiskDurationPerc) ' +/- ' num2str(whiskingData.Blank_SAP.std_whiskDurationPerc)]); disp(' ')
    disp(['SSP-SAP: ' num2str(whiskingData.SSP_SAP.mean_whiskDurationPerc) ' +/- ' num2str(whiskingData.SSP_SAP.std_whiskDurationPerc)]); disp(' ')
    disp(['Naive vs. Blank ttest p = ' num2str(whiskPercStats1.p)]); disp(' ')
    disp(['Blank vs. SAP ttest p = ' num2str(whiskPercStats2.p)]); disp(' ')

    % pupil size
    disp('======================================================================================================================')
    disp(['Naive Awake: ' num2str(pupilSizeData.Naive.zDiameter.Rest.meanData) ' +/- ' num2str(pupilSizeData.Naive.zDiameter.Rest.stdData)]); disp(' ')
    disp(['Naive NREM: ' num2str(pupilSizeData.Naive.zDiameter.NREM.meanData) ' +/- ' num2str(pupilSizeData.Naive.zDiameter.NREM.stdData)]); disp(' ')
    disp(['Naive REM: ' num2str(pupilSizeData.Naive.zDiameter.REM.meanData) ' +/- ' num2str(pupilSizeData.Naive.zDiameter.REM.stdData)]); disp(' ')
    disp(['Blank-SAP Awake: ' num2str(pupilSizeData.Blank_SAP.zDiameter.Rest.meanData) ' +/- ' num2str(pupilSizeData.Blank_SAP.zDiameter.Rest.stdData)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(pupilSizeData.Blank_SAP.zDiameter.NREM.meanData) ' +/- ' num2str(pupilSizeData.Blank_SAP.zDiameter.NREM.stdData)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(pupilSizeData.Blank_SAP.zDiameter.REM.meanData) ' +/- ' num2str(pupilSizeData.Blank_SAP.zDiameter.REM.stdData)]); disp(' ')
    disp(['SSP-SAP Awake: ' num2str(pupilSizeData.SSP_SAP.zDiameter.Rest.meanData) ' +/- ' num2str(pupilSizeData.SSP_SAP.zDiameter.Rest.stdData)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(pupilSizeData.SSP_SAP.zDiameter.NREM.meanData) ' +/- ' num2str(pupilSizeData.SSP_SAP.zDiameter.NREM.stdData)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(pupilSizeData.SSP_SAP.zDiameter.REM.meanData) ' +/- ' num2str(pupilSizeData.SSP_SAP.zDiameter.REM.stdData)]); disp(' ')
    disp(['Naive vs. Blank Awake ttest p = ' num2str(pupilSizeStats.Rest1.p)]); disp(' ')
    disp(['Blank vs. SAP Awake ttest p = ' num2str(pupilSizeStats.Rest2.p)]); disp(' ')
    disp(['Naive vs. Blank NREM ttest p = ' num2str(pupilSizeStats.NREM1.p)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(pupilSizeStats.NREM2.p)]); disp(' ')
    disp(['Naive vs. Blank REM ttest p = ' num2str(pupilSizeStats.REM1.p)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(pupilSizeStats.REM2.p)]); disp(' ')

    % pupil evoked response
    disp('======================================================================================================================')
    disp('Pupil evoked response, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Naive: ' num2str(pupilEvokedData.Naive.zDiameter.contra.meanPeak) ' +/- ' num2str(pupilEvokedData.Naive.zDiameter.contra.stdPeak)]); disp(' ')
    disp(['Blank-SAP: ' num2str(pupilEvokedData.Blank_SAP.zDiameter.contra.meanPeak) ' +/- ' num2str(pupilEvokedData.Blank_SAP.zDiameter.contra.stdPeak)]); disp(' ')
    disp(['SSP-SAP: ' num2str(pupilEvokedData.SSP_SAP.zDiameter.contra.meanPeak) ' +/- ' num2str(pupilEvokedData.SSP_SAP.zDiameter.contra.stdPeak)]); disp(' ')
    disp(['Naive vs. Blank ttest p = ' num2str(pupilEvokedStats1.p)]); disp(' ')
    disp(['Blank vs. SAP ttest p = ' num2str(pupilEvokedStats2.p)]); disp(' ')

    % inter-blink-interval
    disp('======================================================================================================================')
    disp('Inter-blink-interval, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Naive: ' num2str(blinkData.Naive.meanData) ' +/- ' num2str(blinkData.Naive.stdData)]); disp(' ')
    disp(['Blank-SAP: ' num2str(blinkData.Blank_SAP.meanData) ' +/- ' num2str(blinkData.Blank_SAP.stdData)]); disp(' ')
    disp(['SSP-SAP: ' num2str(blinkData.SSP_SAP.meanData) ' +/- ' num2str(blinkData.SSP_SAP.stdData)]); disp(' ')
    disp(['Naive vs. Blank ttest p = ' num2str(BlinkingStats1.p)]); disp(' ')
    disp(['Blank vs. SAP ttest p = ' num2str(BlinkingStats2.p)]); disp(' ')

    diary off
end