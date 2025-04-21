function [] = Fig6_eLife2025(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% Ephys bilateral coherence
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_BilatCoher_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f','binC','binf','group','animalID'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_BilatCoher_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            ephysBilatCoherData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(ephysBilatCoherData.(group).(dataType),behavior) == false
                    for ee = 1:length(variables)
                        variable = variables{1,ee};
                        if any(strcmp(variable,{'group','animalID','binf'})) == true
                            ephysBilatCoherData.(group).(dataType).(behavior).(variable) = {};
                        else
                            ephysBilatCoherData.(group).(dataType).(behavior).(variable) = [];
                        end
                    end
                end
                if isempty(Results_BilatCoher_Ephys.(group).(animalID).(dataType).(behavior).C) == false
                    ephysBilatCoherData.(group).(dataType).(behavior).C = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).C,Results_BilatCoher_Ephys.(group).(animalID).(dataType).(behavior).C');
                    ephysBilatCoherData.(group).(dataType).(behavior).f = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).f,Results_BilatCoher_Ephys.(group).(animalID).(dataType).(behavior).f);
                    freqBand = round(Results_BilatCoher_Ephys.(group).(animalID).(dataType).(behavior).f,2);
                    frequencyList = unique(freqBand);
                    for qq = 1:length(frequencyList)/2
                        freqIdx = find(freqBand == frequencyList(1,qq));
                        ephysBilatCoherData.(group).(dataType).(behavior).binC = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).binC,mean(Results_BilatCoher_Ephys.(group).(animalID).(dataType).(behavior).C(freqIdx)));
                        ephysBilatCoherData.(group).(dataType).(behavior).binf = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).binf,num2str(mean(freqBand(freqIdx))));
                        ephysBilatCoherData.(group).(dataType).(behavior).group = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).group,group);
                        ephysBilatCoherData.(group).(dataType).(behavior).animalID = cat(1,ephysBilatCoherData.(group).(dataType).(behavior).animalID,animalID);
                    end
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                if any(strcmp(variable,{'C','f'})) == true
                    ephysBilatCoherData.(group).(dataType).(behavior).(['mean_' variable]) = mean(ephysBilatCoherData.(group).(dataType).(behavior).(variable),1);
                    ephysBilatCoherData.(group).(dataType).(behavior).(['stdErr_' variable]) = std(ephysBilatCoherData.(group).(dataType).(behavior).(variable),0,1)./sqrt(size(ephysBilatCoherData.(group).(dataType).(behavior).(variable),1));
                end
            end
        end
    end
end
% GLME comparing coherence at all frequencies
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        ephysBilatCoherStats.(dataType).(behavior).tableSize = cat(1,ephysBilatCoherData.Blank_SAP.(dataType).(behavior).binC,ephysBilatCoherData.SSP_SAP.(dataType).(behavior).binC);
        ephysBilatCoherStats.(dataType).(behavior).Table = table('Size',[size(ephysBilatCoherStats.(dataType).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'AnimalID','Treatment','Frequency','Coherence'});
        ephysBilatCoherStats.(dataType).(behavior).Table.AnimalID = cat(1,ephysBilatCoherData.Blank_SAP.(dataType).(behavior).animalID,ephysBilatCoherData.SSP_SAP.(dataType).(behavior).animalID);
        ephysBilatCoherStats.(dataType).(behavior).Table.Treatment = cat(1,ephysBilatCoherData.Blank_SAP.(dataType).(behavior).group,ephysBilatCoherData.SSP_SAP.(dataType).(behavior).group);
        ephysBilatCoherStats.(dataType).(behavior).Table.Frequency = cat(1,ephysBilatCoherData.Blank_SAP.(dataType).(behavior).binf,ephysBilatCoherData.SSP_SAP.(dataType).(behavior).binf);
        ephysBilatCoherStats.(dataType).(behavior).Table.Coherence = cat(1,ephysBilatCoherData.Blank_SAP.(dataType).(behavior).binC,ephysBilatCoherData.SSP_SAP.(dataType).(behavior).binC);
        ephysBilatCoherStats.(dataType).(behavior).FitFormula = 'Coherence ~ 1 + Treatment + (1|Frequency) + (1|AnimalID)';
        ephysBilatCoherStats.(dataType).(behavior).Stats = fitglme(ephysBilatCoherStats.(dataType).(behavior).Table,ephysBilatCoherStats.(dataType).(behavior).FitFormula);
    end
end

%% GCaMP bilateral coherence
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_BilatCoher_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'C','f','binC','binf','group','animalID'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_BilatCoher_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            gcampBilatCoherData.(group).(dataType).dummCheck = 1;
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                if isfield(gcampBilatCoherData.(group).(dataType),behavior) == false
                    for ee = 1:length(variables)
                        variable = variables{1,ee};
                        if any(strcmp(variable,{'group','animalID','binf'})) == true
                            gcampBilatCoherData.(group).(dataType).(behavior).(variable) = {};
                        else
                            gcampBilatCoherData.(group).(dataType).(behavior).(variable) = [];
                        end
                    end
                end
                if isempty(Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).C) == false
                    gcampBilatCoherData.(group).(dataType).(behavior).C = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).C,Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).C');
                    gcampBilatCoherData.(group).(dataType).(behavior).f = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).f,Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).f);
                    freqBand = round(Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).f,2);
                    frequencyList = unique(freqBand);
                    for qq = 1:length(frequencyList)/2
                        freqIdx = find(freqBand == frequencyList(1,qq));
                        gcampBilatCoherData.(group).(dataType).(behavior).binC = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).binC,mean(Results_BilatCoher_GCaMP.(group).(animalID).(dataType).(behavior).C(freqIdx)));
                        gcampBilatCoherData.(group).(dataType).(behavior).binf = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).binf,num2str(mean(freqBand(freqIdx))));
                        gcampBilatCoherData.(group).(dataType).(behavior).group = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).group,group);
                        gcampBilatCoherData.(group).(dataType).(behavior).animalID = cat(1,gcampBilatCoherData.(group).(dataType).(behavior).animalID,animalID);
                    end
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                if any(strcmp(variable,{'C','f'})) == true
                    gcampBilatCoherData.(group).(dataType).(behavior).(['mean_' variable]) = mean(gcampBilatCoherData.(group).(dataType).(behavior).(variable),1);
                    gcampBilatCoherData.(group).(dataType).(behavior).(['stdErr_' variable]) = std(gcampBilatCoherData.(group).(dataType).(behavior).(variable),0,1)./sqrt(size(gcampBilatCoherData.(group).(dataType).(behavior).(variable),1));
                end
            end
        end
    end
end
% GLME comparing coherence at all frequencies
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        gcampBilatCoherStats.(dataType).(behavior).tableSize = cat(1,gcampBilatCoherData.Blank_SAP.(dataType).(behavior).binC,gcampBilatCoherData.SSP_SAP.(dataType).(behavior).binC);
        gcampBilatCoherStats.(dataType).(behavior).Table = table('Size',[size(gcampBilatCoherStats.(dataType).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'AnimalID','Treatment','Frequency','Coherence'});
        gcampBilatCoherStats.(dataType).(behavior).Table.AnimalID = cat(1,gcampBilatCoherData.Blank_SAP.(dataType).(behavior).animalID,gcampBilatCoherData.SSP_SAP.(dataType).(behavior).animalID);
        gcampBilatCoherStats.(dataType).(behavior).Table.Treatment = cat(1,gcampBilatCoherData.Blank_SAP.(dataType).(behavior).group,gcampBilatCoherData.SSP_SAP.(dataType).(behavior).group);
        gcampBilatCoherStats.(dataType).(behavior).Table.Frequency = cat(1,gcampBilatCoherData.Blank_SAP.(dataType).(behavior).binf,gcampBilatCoherData.SSP_SAP.(dataType).(behavior).binf);
        gcampBilatCoherStats.(dataType).(behavior).Table.Coherence = cat(1,gcampBilatCoherData.Blank_SAP.(dataType).(behavior).binC,gcampBilatCoherData.SSP_SAP.(dataType).(behavior).binC);
        gcampBilatCoherStats.(dataType).(behavior).FitFormula = 'Coherence ~ 1 + Treatment + (1|Frequency) + (1|AnimalID)';
        gcampBilatCoherStats.(dataType).(behavior).Stats = fitglme(gcampBilatCoherStats.(dataType).(behavior).Table,gcampBilatCoherStats.(dataType).(behavior).FitFormula);
    end
end

%% Figure 6
Fig6 = figure('Name','Fig. 6');

% Ephys bilateral HbT - Rest
subplot(4,3,1);
p1 = semilogx(ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_C + ephysBilatCoherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Rest.mean_C - ephysBilatCoherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('black'),'LineWidth',0.25);
p2 = semilogx(ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_C + ephysBilatCoherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Rest.mean_C - ephysBilatCoherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbT Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
legend([p1,p2],'Blank-SAP','SSP-SAP')
xlim([0.1,0.5])
ylim([0.6,1])
set(gca,'box','off')
axis square

% Ephys bilateral HbT - Alert
subplot(4,3,2);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_C + ephysBilatCoherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Alert.mean_C - ephysBilatCoherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_C + ephysBilatCoherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Alert.mean_C - ephysBilatCoherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbT Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.01,0.5])
ylim([0.6,1])
set(gca,'box','off')
axis square

% Ephys bilateral HbT - Asleep
subplot(4,3,3);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_C + ephysBilatCoherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.HbT.Asleep.mean_C - ephysBilatCoherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_C + ephysBilatCoherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.HbT.Asleep.mean_C - ephysBilatCoherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat HbT Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.01,0.5])
ylim([0.6,1])
set(gca,'box','off')
axis square

% Ephys bilateral gammaBandPower - Rest
subplot(4,3,4);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_C + ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.mean_C - ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_C + ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.mean_C - ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat gammaBandPower Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.1,0.5])
ylim([0,1])
set(gca,'box','off')
axis square

% Ephys bilateral gammaBandPower - Alert
subplot(4,3,5);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_C + ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.mean_C - ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_C + ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.mean_C - ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat gammaBandPower Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.01,0.5])
ylim([0,1])
set(gca,'box','off')
axis square

% Ephys bilateral gammaBandPower - Asleep
subplot(4,3,6);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_C + ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.mean_C - ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_C,'color',colors('electric purple'),'LineWidth',2);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_C + ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
semilogx(ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_f,ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.mean_C - ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.stdErr_C,'color',colors('electric purple'),'LineWidth',0.25);
title('bilat gammaBandPower Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.01,0.5])
ylim([0,1])
set(gca,'box','off')
axis square

% GCaMP bilateral HbT - Rest
subplot(4,3,7);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_C + gcampBilatCoherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Rest.mean_C - gcampBilatCoherData.Blank_SAP.HbT.Rest.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_C + gcampBilatCoherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Rest.mean_C - gcampBilatCoherData.SSP_SAP.HbT.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
title('bilat HbT Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.1,0.5])
ylim([0.7,1])
set(gca,'box','off')
axis square

% GCaMP bilateral HbT - Alert
subplot(4,3,8);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_C + gcampBilatCoherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Alert.mean_C - gcampBilatCoherData.Blank_SAP.HbT.Alert.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_C + gcampBilatCoherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Alert.mean_C - gcampBilatCoherData.SSP_SAP.HbT.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
title('bilat HbT Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.01,0.5])
ylim([0.7,1])
set(gca,'box','off')
axis square

% GCaMP bilateral HbT - Asleep
subplot(4,3,9);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_C + gcampBilatCoherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.HbT.Asleep.mean_C - gcampBilatCoherData.Blank_SAP.HbT.Asleep.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_C + gcampBilatCoherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.HbT.Asleep.mean_C - gcampBilatCoherData.SSP_SAP.HbT.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
title('bilat HbT Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.01,0.5])
ylim([0.7,1])
set(gca,'box','off')
axis square

% GCaMP bilateral GCaMP - Rest
subplot(4,3,10);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_C + gcampBilatCoherData.Blank_SAP.GCaMP.Rest.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Rest.mean_C - gcampBilatCoherData.Blank_SAP.GCaMP.Rest.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_C,'color',colors('north texas green'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_C + gcampBilatCoherData.SSP_SAP.GCaMP.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Rest.mean_C - gcampBilatCoherData.SSP_SAP.GCaMP.Rest.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
title('bilat GCaMP Rest')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.1,0.5])
ylim([0.6,1])
set(gca,'box','off')
axis square

% GCaMP bilateral GCaMP - Alert
subplot(4,3,11);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_C + gcampBilatCoherData.Blank_SAP.GCaMP.Alert.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Alert.mean_C - gcampBilatCoherData.Blank_SAP.GCaMP.Alert.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_C,'color',colors('north texas green'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_C + gcampBilatCoherData.SSP_SAP.GCaMP.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Alert.mean_C - gcampBilatCoherData.SSP_SAP.GCaMP.Alert.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
title('bilat GCaMP Alert')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.01,0.5])
ylim([0.6,1])
set(gca,'box','off')
axis square

% GCaMP bilateral GCaMP - Asleep
subplot(4,3,12);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_C,'color',colors('black'),'LineWidth',2);
hold on
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_C + gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.mean_C - gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.stdErr_C,'color',colors('black'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_C,'color',colors('north texas green'),'LineWidth',2);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_C + gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
semilogx(gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_f,gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.mean_C - gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.stdErr_C,'color',colors('north texas green'),'LineWidth',0.25);
title('bilat GCaMP Asleep')
xlabel('Freq (Hz)')
ylabel('Coherence')
xlim([0.01,0.5])
ylim([0.6,1])
set(gca,'box','off')
axis square

%% Save figure and stats
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figs/Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig6,[dirpath 'Fig6']);
    diaryFile = [dirpath 'Fig6_Stats.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    % statistical diary
    diary(diaryFile)
    diary on

    comparisons = 1;
    alphaA = 0.05/comparisons;
    alphaB = 0.01/comparisons;
    alphaC = 0.001/comparisons;

    % HbT coherence (Rest)
    disp('======================================================================================================================')
    disp('HbT coherence (Rest), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(ephysBilatCoherData.Blank_SAP.HbT.Rest.binC)) ' +/- ' num2str(std(ephysBilatCoherData.Blank_SAP.HbT.Rest.binC,0,1)./sqrt(size(ephysBilatCoherData.Blank_SAP.HbT.Rest.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(ephysBilatCoherData.SSP_SAP.HbT.Rest.binC)) ' +/- ' num2str(std(ephysBilatCoherData.SSP_SAP.HbT.Rest.binC,0,1)./sqrt(size(ephysBilatCoherData.SSP_SAP.HbT.Rest.binC,1)))]); disp(' ')
    disp('GLME statistics for Rest HbT [Ephys]')
    disp(ephysBilatCoherStats.HbT.Rest.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % HbT coherence (Alert)
    disp('======================================================================================================================')
    disp('HbT coherence (Alert), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(ephysBilatCoherData.Blank_SAP.HbT.Alert.binC)) ' +/- ' num2str(std(ephysBilatCoherData.Blank_SAP.HbT.Alert.binC,0,1)./sqrt(size(ephysBilatCoherData.Blank_SAP.HbT.Alert.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(ephysBilatCoherData.SSP_SAP.HbT.Alert.binC)) ' +/- ' num2str(std(ephysBilatCoherData.SSP_SAP.HbT.Alert.binC,0,1)./sqrt(size(ephysBilatCoherData.SSP_SAP.HbT.Alert.binC,1)))]); disp(' ')
    disp('GLME statistics for Alert HbT [Ephys]')
    disp(ephysBilatCoherStats.HbT.Alert.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % HbT coherence (Asleep)
    disp('======================================================================================================================')
    disp('HbT coherence (Asleep), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(ephysBilatCoherData.Blank_SAP.HbT.Asleep.binC)) ' +/- ' num2str(std(ephysBilatCoherData.Blank_SAP.HbT.Asleep.binC,0,1)./sqrt(size(ephysBilatCoherData.Blank_SAP.HbT.Asleep.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(ephysBilatCoherData.SSP_SAP.HbT.Asleep.binC)) ' +/- ' num2str(std(ephysBilatCoherData.SSP_SAP.HbT.Asleep.binC,0,1)./sqrt(size(ephysBilatCoherData.SSP_SAP.HbT.Asleep.binC,1)))]); disp(' ')
    disp('GLME statistics for Asleep HbT [Ephys]')
    disp(ephysBilatCoherStats.HbT.Asleep.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % Gamma coherence (Rest)
    disp('======================================================================================================================')
    disp('Gamma coherence (Rest), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.binC)) ' +/- ' num2str(std(ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.binC,0,1)./sqrt(size(ephysBilatCoherData.Blank_SAP.gammaBandPower.Rest.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.binC)) ' +/- ' num2str(std(ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.binC,0,1)./sqrt(size(ephysBilatCoherData.SSP_SAP.gammaBandPower.Rest.binC,1)))]); disp(' ')
    disp('GLME statistics for Rest gamma [Ephys]')
    disp(ephysBilatCoherStats.gammaBandPower.Rest.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % Gamma coherence (Alert)
    disp('======================================================================================================================')
    disp('HbT coherence (Alert), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.binC)) ' +/- ' num2str(std(ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.binC,0,1)./sqrt(size(ephysBilatCoherData.Blank_SAP.gammaBandPower.Alert.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.binC)) ' +/- ' num2str(std(ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.binC,0,1)./sqrt(size(ephysBilatCoherData.SSP_SAP.gammaBandPower.Alert.binC,1)))]); disp(' ')
    disp('GLME statistics for Alert gamma [Ephys]')
    disp(ephysBilatCoherStats.gammaBandPower.Alert.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % Gamma coherence (Asleep)
    disp('======================================================================================================================')
    disp('HbT coherence (Asleep), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.binC)) ' +/- ' num2str(std(ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.binC,0,1)./sqrt(size(ephysBilatCoherData.Blank_SAP.gammaBandPower.Asleep.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.binC)) ' +/- ' num2str(std(ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.binC,0,1)./sqrt(size(ephysBilatCoherData.SSP_SAP.gammaBandPower.Asleep.binC,1)))]); disp(' ')
    disp('GLME statistics for Asleep gamma [Ephys]')
    disp(ephysBilatCoherStats.gammaBandPower.Asleep.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % HbT coherence (Rest)
    disp('======================================================================================================================')
    disp('HbT coherence (Rest), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampBilatCoherData.Blank_SAP.HbT.Rest.binC)) ' +/- ' num2str(std(gcampBilatCoherData.Blank_SAP.HbT.Rest.binC,0,1)./sqrt(size(gcampBilatCoherData.Blank_SAP.HbT.Rest.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampBilatCoherData.SSP_SAP.HbT.Rest.binC)) ' +/- ' num2str(std(gcampBilatCoherData.SSP_SAP.HbT.Rest.binC,0,1)./sqrt(size(gcampBilatCoherData.SSP_SAP.HbT.Rest.binC,1)))]); disp(' ')
    disp('GLME statistics for Rest HbT [GCaMP]')
    disp(gcampBilatCoherStats.HbT.Rest.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % HbT coherence (Alert)
    disp('======================================================================================================================')
    disp('HbT coherence (Alert), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampBilatCoherData.Blank_SAP.HbT.Alert.binC)) ' +/- ' num2str(std(gcampBilatCoherData.Blank_SAP.HbT.Alert.binC,0,1)./sqrt(size(gcampBilatCoherData.Blank_SAP.HbT.Alert.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampBilatCoherData.SSP_SAP.HbT.Alert.binC)) ' +/- ' num2str(std(gcampBilatCoherData.SSP_SAP.HbT.Alert.binC,0,1)./sqrt(size(gcampBilatCoherData.SSP_SAP.HbT.Alert.binC,1)))]); disp(' ')
    disp('GLME statistics for Alert HbT [GCaMP]')
    disp(gcampBilatCoherStats.HbT.Alert.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % HbT coherence (Asleep)
    disp('======================================================================================================================')
    disp('HbT coherence (Asleep), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampBilatCoherData.Blank_SAP.HbT.Asleep.binC)) ' +/- ' num2str(std(gcampBilatCoherData.Blank_SAP.HbT.Asleep.binC,0,1)./sqrt(size(gcampBilatCoherData.Blank_SAP.HbT.Asleep.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampBilatCoherData.SSP_SAP.HbT.Asleep.binC)) ' +/- ' num2str(std(gcampBilatCoherData.SSP_SAP.HbT.Asleep.binC,0,1)./sqrt(size(gcampBilatCoherData.SSP_SAP.HbT.Asleep.binC,1)))]); disp(' ')
    disp('GLME statistics for Asleep HbT [GCaMP]')
    disp(gcampBilatCoherStats.HbT.Asleep.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % GCaMP coherence (Rest)
    disp('======================================================================================================================')
    disp('HbT coherence (Rest), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampBilatCoherData.Blank_SAP.GCaMP.Rest.binC)) ' +/- ' num2str(std(gcampBilatCoherData.Blank_SAP.GCaMP.Rest.binC,0,1)./sqrt(size(gcampBilatCoherData.Blank_SAP.GCaMP.Rest.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampBilatCoherData.SSP_SAP.GCaMP.Rest.binC)) ' +/- ' num2str(std(gcampBilatCoherData.SSP_SAP.GCaMP.Rest.binC,0,1)./sqrt(size(gcampBilatCoherData.SSP_SAP.GCaMP.Rest.binC,1)))]); disp(' ')
    disp('GLME statistics for Rest F/F [GCaMP]')
    disp(gcampBilatCoherStats.GCaMP.Rest.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % GCaMP coherence (Alert)
    disp('======================================================================================================================')
    disp('HbT coherence (Alert), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampBilatCoherData.Blank_SAP.GCaMP.Alert.binC)) ' +/- ' num2str(std(gcampBilatCoherData.Blank_SAP.GCaMP.Alert.binC,0,1)./sqrt(size(gcampBilatCoherData.Blank_SAP.GCaMP.Alert.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampBilatCoherData.SSP_SAP.GCaMP.Alert.binC)) ' +/- ' num2str(std(gcampBilatCoherData.SSP_SAP.GCaMP.Alert.binC,0,1)./sqrt(size(gcampBilatCoherData.SSP_SAP.GCaMP.Alert.binC,1)))]); disp(' ')
    disp('GLME statistics for Alert F/F [GCaMP]')
    disp(gcampBilatCoherStats.GCaMP.Alert.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % GCaMP coherence (Asleep)
    disp('======================================================================================================================')
    disp('HbT coherence (Asleep), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.binC)) ' +/- ' num2str(std(gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.binC,0,1)./sqrt(size(gcampBilatCoherData.Blank_SAP.GCaMP.Asleep.binC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.binC)) ' +/- ' num2str(std(gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.binC,0,1)./sqrt(size(gcampBilatCoherData.SSP_SAP.GCaMP.Asleep.binC,1)))]); disp(' ')
    disp('GLME statistics for Asleep F/F [GCaMP]')
    disp(gcampBilatCoherStats.GCaMP.Asleep.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    diary off
end