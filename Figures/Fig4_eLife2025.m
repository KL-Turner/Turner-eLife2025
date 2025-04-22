function [] = Fig4_eLife2025(rootFolder,saveState,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% LFP power
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PowerSpec_LFP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Alert','Asleep','All'};
variables = {'S','f','deltaS'};
dimensions = [2,1,1];
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PowerSpec_LFP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                data.(group).(hemisphere).(behavior).dummCheck = 1;
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    dimension = dimensions(ee);
                    if isfield(data.(group).(hemisphere).(behavior),(variable)) == false
                        data.(group).(hemisphere).(behavior).(variable) = [];
                        data.(group).(hemisphere).(behavior).group = {};
                        data.(group).(hemisphere).(behavior).animalID = {};
                        data.(group).(hemisphere).(behavior).hemisphere = {};
                        data.(group).(hemisphere).(behavior).behavior = {};
                    end
                    % pull data if field isn't empty
                    if isempty(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S) == false
                        if strcmp(variable,'deltaS') == true
                            index = find(round(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).f,2) == 4);
                            deltaIndex = index(end);
                            data.(group).(hemisphere).(behavior).(variable) = cat(dimension,data.(group).(hemisphere).(behavior).(variable),mean(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S(1:deltaIndex)));
                            % for stats
                            data.(group).(hemisphere).(behavior).group = cat(1,data.(group).(hemisphere).(behavior).group,group);
                            data.(group).(hemisphere).(behavior).animalID = cat(1,data.(group).(hemisphere).(behavior).animalID,animalID);
                            data.(group).(hemisphere).(behavior).hemisphere = cat(1,data.(group).(hemisphere).(behavior).hemisphere,hemisphere);
                            data.(group).(hemisphere).(behavior).behavior = cat(1,data.(group).(hemisphere).(behavior).behavior,behavior);
                        else
                            data.(group).(hemisphere).(behavior).(variable) = cat(dimension,data.(group).(hemisphere).(behavior).(variable),Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).(variable));
                        end
                    end
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
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                dimension = dimensions(dd);
                data.(group).(hemisphere).(behavior).(['mean_' variable]) = mean(data.(group).(hemisphere).(behavior).(variable),dimension,'omitnan');
                data.(group).(hemisphere).(behavior).(['stdErr_' variable]) = std(data.(group).(hemisphere).(behavior).(variable),0,dimension,'omitnan')./sqrt(size(data.(group).(hemisphere).(behavior).(variable),dimension));
            end
        end
    end
end

%% statistics - generalized linear mixed effects model
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        lfpStats.(hemisphere).(behavior).tableSize = cat(1,data.Blank_SAP.(hemisphere).(behavior).deltaS,data.SSP_SAP.(hemisphere).(behavior).deltaS,data.Naive.(hemisphere).(behavior).deltaS);
        lfpStats.(hemisphere).(behavior).Table = table('Size',[size(lfpStats.(hemisphere).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'group','animalID','behavior','deltaS'});
        lfpStats.(hemisphere).(behavior).Table.group = cat(1,data.Blank_SAP.(hemisphere).(behavior).group,data.SSP_SAP.(hemisphere).(behavior).group,data.Naive.(hemisphere).(behavior).group);
        lfpStats.(hemisphere).(behavior).Table.animalID = cat(1,data.Blank_SAP.(hemisphere).(behavior).animalID,data.SSP_SAP.(hemisphere).(behavior).animalID,data.Naive.(hemisphere).(behavior).animalID);
        lfpStats.(hemisphere).(behavior).Table.behavior = cat(1,data.Blank_SAP.(hemisphere).(behavior).behavior,data.SSP_SAP.(hemisphere).(behavior).behavior,data.Naive.(hemisphere).(behavior).behavior);
        lfpStats.(hemisphere).(behavior).Table.deltaS = cat(1,data.Blank_SAP.(hemisphere).(behavior).deltaS,data.SSP_SAP.(hemisphere).(behavior).deltaS,data.Naive.(hemisphere).(behavior).deltaS);
        lfpStats.(hemisphere).(behavior).FitFormula = 'deltaS ~ 1 + group + behavior + (1|animalID)';
        lfpStats.(hemisphere).(behavior).Stats = fitglme(lfpStats.(hemisphere).(behavior).Table,lfpStats.(hemisphere).(behavior).FitFormula);
    end
end

%% IOS ephys brief whisker stim
path = [rootFolder delim 'Results_Turner'];
cd(path)
resultsStruct = 'Results_Evoked_Ephys';
load(resultsStruct);
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','cortMUA','cortGam','cortS','T','F','timeVector','group','animalID','snipHbT','snipGam'}; % temporary
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(solenoids)
                solenoid = solenoids{1,dd};
                iosEphysData.(group).(hemisphere).(solenoid).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(iosEphysData.(group).(hemisphere).(solenoid),dataType) == false
                        if any(strcmp(dataType,{'group','animalID'})) == true
                            iosEphysData.(group).(hemisphere).(solenoid).(dataType) = {};
                        else
                            iosEphysData.(group).(hemisphere).(solenoid).(dataType) = [];
                        end
                    end
                end
                iosEphysData.(group).(hemisphere).(solenoid).HbT = cat(1,iosEphysData.(group).(hemisphere).(solenoid).HbT,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).HbT);
                iosEphysData.(group).(hemisphere).(solenoid).cortMUA = cat(1,iosEphysData.(group).(hemisphere).(solenoid).cortMUA,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortMUA);
                iosEphysData.(group).(hemisphere).(solenoid).cortGam = cat(1,iosEphysData.(group).(hemisphere).(solenoid).cortGam,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortGam);
                iosEphysData.(group).(hemisphere).(solenoid).cortS = cat(3,iosEphysData.(group).(hemisphere).(solenoid).cortS,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortLFP);
                iosEphysData.(group).(hemisphere).(solenoid).T = cat(1,iosEphysData.(group).(hemisphere).(solenoid).T,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).T);
                iosEphysData.(group).(hemisphere).(solenoid).F = cat(1,iosEphysData.(group).(hemisphere).(solenoid).F,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).F);
                iosEphysData.(group).(hemisphere).(solenoid).timeVector = cat(1,iosEphysData.(group).(hemisphere).(solenoid).timeVector,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).timeVector);
                iosEphysData.(group).(hemisphere).(solenoid).group = cat(1,iosEphysData.(group).(hemisphere).(solenoid).group,group);
                iosEphysData.(group).(hemisphere).(solenoid).animalID = cat(1,iosEphysData.(group).(hemisphere).(solenoid).animalID,animalID);
                % temporary
                iosEphysData.(group).(hemisphere).(solenoid).snipHbT = cat(1,iosEphysData.(group).(hemisphere).(solenoid).snipHbT,max(Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).HbT(60:150)));
                iosEphysData.(group).(hemisphere).(solenoid).snipGam = cat(1,iosEphysData.(group).(hemisphere).(solenoid).snipGam,max(Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortGam(60:150)));
            end
        end
    end
end
% pair stimulation with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            [comparison] = FindSolenoidComparison(hemisphere,solenoid);
            iosEphysData.(group).(hemisphere).(comparison).group = {};
            iosEphysData.(group).(hemisphere).(comparison).animalID = {};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if any(strcmp(dataType,{'group','animalID'})) == true
                    iosEphysData.(group).(hemisphere).(comparison).(dataType) = iosEphysData.(group).(hemisphere).(solenoid).(dataType);
                elseif strcmp(dataType,'cortS')
                    iosEphysData.(group).(hemisphere).(comparison).(dataType) = iosEphysData.(group).(hemisphere).(solenoid).(dataType);
                    iosEphysData.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(iosEphysData.(group).(hemisphere).(solenoid).(dataType),3).*100;
                else
                    iosEphysData.(group).(hemisphere).(comparison).(dataType) = iosEphysData.(group).(hemisphere).(solenoid).(dataType);
                    iosEphysData.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(iosEphysData.(group).(hemisphere).(solenoid).(dataType),1);
                    iosEphysData.(group).(hemisphere).(comparison).(['stdErr_' dataType]) = std(iosEphysData.(group).(hemisphere).(solenoid).(dataType),1)./sqrt(size(iosEphysData.(group).(hemisphere).(solenoid).(dataType),1));
                    if strcmp(dataType,'HbT') == true
                        for ee = 1:size(iosEphysData.(group).(hemisphere).(comparison).HbT,1)
                            startIdx = find(iosEphysData.(group).(hemisphere).(solenoid).timeVector(ee,:) == 2);
                            endIdx =  find(iosEphysData.(group).(hemisphere).(solenoid).timeVector(ee,:) == 4);
                            hbtSnip = iosEphysData.(group).(hemisphere).(comparison).HbT(ee,:);
                            iosEphysData.(group).(hemisphere).(comparison).AUC(ee,1) = mean(hbtSnip(startIdx:endIdx));
                        end
                    end
                end
            end
        end
    end
end
% statistics - generalized linear mixed effects model
iosEphysStats.tableSize = cat(1,iosEphysData.Blank_SAP.RH.contra.AUC,iosEphysData.SSP_SAP.RH.contra.AUC);
iosEphysStats.Table = table('Size',[size(iosEphysStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
iosEphysStats.Table.Mouse = cat(1,iosEphysData.Blank_SAP.RH.contra.animalID,iosEphysData.SSP_SAP.RH.contra.animalID);
iosEphysStats.Table.Group = cat(1,iosEphysData.Blank_SAP.RH.contra.group,iosEphysData.SSP_SAP.RH.contra.group);
iosEphysStats.Table.AUC = cat(1,iosEphysData.Blank_SAP.RH.contra.AUC,iosEphysData.SSP_SAP.RH.contra.AUC);
iosEphysStats.FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
iosEphysStats.Stats = fitglme(iosEphysStats.Table,iosEphysStats.FitFormula);

%% IOS GCaMP
path = [rootFolder delim 'Results_Turner'];
cd(path)
resultsStruct = 'Results_Evoked_GCaMP';
load(resultsStruct);
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','HbO','HbR','GCaMP','timeVector','group','animalID','snipHbT','snipGCaMP'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(solenoids)
                solenoid = solenoids{1,dd};
                gcampdata.(group).(hemisphere).(solenoid).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(gcampdata.(group).(hemisphere).(solenoid),dataType) == false
                        if any(strcmp(dataType,{'group','animalID'})) == true
                            gcampdata.(group).(hemisphere).(solenoid).(dataType) = {};
                        else
                            gcampdata.(group).(hemisphere).(solenoid).(dataType) = [];
                        end
                    end
                end
                gcampdata.(group).(hemisphere).(solenoid).HbT = cat(1,gcampdata.(group).(hemisphere).(solenoid).HbT,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbT.mean);
                gcampdata.(group).(hemisphere).(solenoid).HbO = cat(1,gcampdata.(group).(hemisphere).(solenoid).HbO,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbO.mean);
                gcampdata.(group).(hemisphere).(solenoid).HbR = cat(1,gcampdata.(group).(hemisphere).(solenoid).HbR,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbR.mean);
                gcampdata.(group).(hemisphere).(solenoid).GCaMP = cat(1,gcampdata.(group).(hemisphere).(solenoid).GCaMP,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).GCaMP.mean*100);
                gcampdata.(group).(hemisphere).(solenoid).timeVector = cat(1,gcampdata.(group).(hemisphere).(solenoid).timeVector,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbT.timeVector);
                gcampdata.(group).(hemisphere).(solenoid).group = cat(1,gcampdata.(group).(hemisphere).(solenoid).group,group);
                gcampdata.(group).(hemisphere).(solenoid).animalID = cat(1,gcampdata.(group).(hemisphere).(solenoid).animalID,animalID);
                % temporary
                gcampdata.(group).(hemisphere).(solenoid).snipHbT = cat(1,gcampdata.(group).(hemisphere).(solenoid).snipHbT,mean(Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbT.mean(40:70)));
                gcampdata.(group).(hemisphere).(solenoid).snipGCaMP = cat(1,gcampdata.(group).(hemisphere).(solenoid).snipGCaMP,mean(Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).GCaMP.mean(35:85)));
            end
        end
    end
end
% pair stimulation with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            [comparison] = FindSolenoidComparison(hemisphere,solenoid);
            gcampdata.(group).(hemisphere).(comparison).group = {};
            gcampdata.(group).(hemisphere).(comparison).animalID = {};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if any(strcmp(dataType,{'group','animalID'})) == true
                    gcampdata.(group).(hemisphere).(comparison).(dataType) = gcampdata.(group).(hemisphere).(solenoid).(dataType);
                else
                    gcampdata.(group).(hemisphere).(comparison).(dataType) = gcampdata.(group).(hemisphere).(solenoid).(dataType);
                    gcampdata.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(gcampdata.(group).(hemisphere).(solenoid).(dataType),1);
                    gcampdata.(group).(hemisphere).(comparison).(['stdErr_' dataType]) = std(gcampdata.(group).(hemisphere).(solenoid).(dataType),1)./sqrt(size(gcampdata.(group).(hemisphere).(solenoid).(dataType),1));
                    if strcmp(dataType,{'timeVector','snipHbT','snipGCaMP'}') == false
                        for ee = 1:size(gcampdata.(group).(hemisphere).(comparison).(dataType),1)
                            if strcmp(dataType,'GCaMP') == true
                                startIdx = find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 2);
                                endIdx =  find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 5);
                            else
                                startIdx = find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 1.5);
                                endIdx =  find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 6.5);
                            end
                            dataSnip = gcampdata.(group).(hemisphere).(comparison).(dataType)(ee,:);
                            gcampdata.(group).(hemisphere).(comparison).(['AUC_' dataType])(ee,1) = mean(dataSnip(startIdx:endIdx));
                        end
                    end
                end
            end
        end
    end
end
% statistics - generalized linear mixed effects model
specDT = {'HbT','HbO','HbR','GCaMP'};
for aa = 1:length(specDT)
    gcampStats.(specDT{1,aa}).tableSize = cat(1,gcampdata.Blank_SAP.RH.contra.(['AUC_' specDT{1,aa}]),gcampdata.SSP_SAP.RH.contra.(['AUC_' specDT{1,aa}]));
    gcampStats.(specDT{1,aa}).Table = table('Size',[size(gcampStats.(specDT{1,aa}).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
    gcampStats.(specDT{1,aa}).Table.Mouse = cat(1,gcampdata.Blank_SAP.RH.contra.animalID,gcampdata.SSP_SAP.RH.contra.animalID);
    gcampStats.(specDT{1,aa}).Table.Group = cat(1,gcampdata.Blank_SAP.RH.contra.group,gcampdata.SSP_SAP.RH.contra.group);
    gcampStats.(specDT{1,aa}).Table.AUC = cat(1,gcampdata.Blank_SAP.RH.contra.(['AUC_' specDT{1,aa}]),gcampdata.SSP_SAP.RH.contra.(['AUC_' specDT{1,aa}]));
    gcampStats.(specDT{1,aa}).FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
    gcampStats.(specDT{1,aa}).Stats = fitglme(gcampStats.(specDT{1,aa}).Table,gcampStats.(specDT{1,aa}).FitFormula);
end

%% Gamma-HbT stats
x = [iosEphysData.Blank_SAP.RH.contra.snipGam; iosEphysData.SSP_SAP.RH.contra.snipGam];
y = [iosEphysData.Blank_SAP.RH.contra.snipHbT; iosEphysData.SSP_SAP.RH.contra.snipHbT];
group = [ones(size(iosEphysData.Blank_SAP.RH.contra.snipGam)); 2 * ones(size(iosEphysData.SSP_SAP.RH.contra.snipGam))]; % 1 for first dataset, 2 for second
tbl = table(x, y, group);
tbl.group = categorical(tbl.group); % Convert group to categorical variable
gammaMdl = fitlm(tbl, 'y ~ x*group - 1');

%% GCaMP-HbT stats
x = [gcampdata.Blank_SAP.RH.contra.snipGCaMP; gcampdata.SSP_SAP.RH.contra.snipGCaMP];
y = [gcampdata.Blank_SAP.RH.contra.snipHbT; gcampdata.SSP_SAP.RH.contra.snipHbT];
group = [ones(size(gcampdata.Blank_SAP.RH.contra.snipGCaMP)); 2 * ones(size(gcampdata.SSP_SAP.RH.contra.snipGCaMP))]; % 1 for first dataset, 2 for second
tbl = table(x, y, group);
tbl.group = categorical(tbl.group); % Convert group to categorical variable
gcampMdl = fitlm(tbl, 'y ~ x*group - 1');

%% Figure 4
Fig4 = figure('Name','Fig. 4');

% Ephys alert LFP
ax1 = subplot(2,3,1);
loglog(data.Blank_SAP.RH.Alert.mean_f,data.Blank_SAP.RH.Alert.mean_S,'color',colors('black'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.RH.Alert.mean_f,data.Blank_SAP.RH.Alert.mean_S + data.Blank_SAP.RH.Alert.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(data.Blank_SAP.RH.Alert.mean_f,data.Blank_SAP.RH.Alert.mean_S - data.Blank_SAP.RH.Alert.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.Alert.mean_f,data.SSP_SAP.RH.Alert.mean_S,'color',colors('sapphire'),'LineWidth',2);
loglog(data.SSP_SAP.RH.Alert.mean_f,data.SSP_SAP.RH.Alert.mean_S + data.SSP_SAP.RH.Alert.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.Alert.mean_f,data.SSP_SAP.RH.Alert.mean_S - data.SSP_SAP.RH.Alert.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);

% Ephys asleep LFP
ax2 = subplot(2,3,2);
loglog(data.Blank_SAP.RH.Asleep.mean_f,data.Blank_SAP.RH.Asleep.mean_S,'color',colors('black'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.RH.Asleep.mean_f,data.Blank_SAP.RH.Asleep.mean_S + data.Blank_SAP.RH.Asleep.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(data.Blank_SAP.RH.Asleep.mean_f,data.Blank_SAP.RH.Asleep.mean_S - data.Blank_SAP.RH.Asleep.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.Asleep.mean_f,data.SSP_SAP.RH.Asleep.mean_S,'color',colors('sapphire'),'LineWidth',2);
loglog(data.SSP_SAP.RH.Asleep.mean_f,data.SSP_SAP.RH.Asleep.mean_S + data.SSP_SAP.RH.Asleep.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.Asleep.mean_f,data.SSP_SAP.RH.Asleep.mean_S - data.SSP_SAP.RH.Asleep.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);

% Ephys all LFP
ax3 = subplot(2,3,3);
loglog(data.Blank_SAP.RH.All.mean_f,data.Blank_SAP.RH.All.mean_S,'color',colors('black'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.RH.All.mean_f,data.Blank_SAP.RH.All.mean_S + data.Blank_SAP.RH.All.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(data.Blank_SAP.RH.All.mean_f,data.Blank_SAP.RH.All.mean_S - data.Blank_SAP.RH.All.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.All.mean_f,data.SSP_SAP.RH.All.mean_S,'color',colors('sapphire'),'LineWidth',2);
loglog(data.SSP_SAP.RH.All.mean_f,data.SSP_SAP.RH.All.mean_S + data.SSP_SAP.RH.All.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.All.mean_f,data.SSP_SAP.RH.All.mean_S - data.SSP_SAP.RH.All.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);
linkaxes([ax1,ax2,ax3],'xy')

% Gamma-HbT
subplot(2,3,4)
s1 = scatter(iosEphysData.Blank_SAP.RH.contra.snipGam,iosEphysData.Blank_SAP.RH.contra.snipHbT,'filled','MarkerFaceColor',colors('black'));
hold on;
p = iosEphysData.Blank_SAP.RH.contra.snipGam(:)\iosEphysData.Blank_SAP.RH.contra.snipHbT(:);
pp = p*iosEphysData.Blank_SAP.RH.contra.snipGam;
plot(iosEphysData.Blank_SAP.RH.contra.snipGam,pp,'color',colors('black'));
s2 = scatter(iosEphysData.SSP_SAP.RH.contra.snipGam,iosEphysData.SSP_SAP.RH.contra.snipHbT,'filled','MarkerFaceColor',colors('sapphire'));
p = iosEphysData.SSP_SAP.RH.contra.snipGam(:)\iosEphysData.SSP_SAP.RH.contra.snipHbT(:);
pp = p*iosEphysData.SSP_SAP.RH.contra.snipGam;
plot(iosEphysData.SSP_SAP.RH.contra.snipGam,pp,'color',colors('sapphire'));
xlabel('Gamma (%)')
ylabel('\DeltaHbT (\muM)')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
axis square

% GCaMP-HbT
subplot(2,3,5)
s1 = scatter(gcampdata.Blank_SAP.RH.contra.snipGCaMP*100,gcampdata.Blank_SAP.RH.contra.snipHbT,'filled','MarkerFaceColor',colors('black'));
hold on;
p = gcampdata.Blank_SAP.RH.contra.snipGCaMP(:)*100\gcampdata.Blank_SAP.RH.contra.snipHbT(:);
pp = p*gcampdata.Blank_SAP.RH.contra.snipGCaMP*100;
plot(gcampdata.Blank_SAP.RH.contra.snipGCaMP*100,pp,'color',colors('black'));
s2 = scatter(gcampdata.SSP_SAP.RH.contra.snipGCaMP*100,gcampdata.SSP_SAP.RH.contra.snipHbT,'filled','MarkerFaceColor',colors('north texas green'));
p = gcampdata.SSP_SAP.RH.contra.snipGCaMP(:)*100\gcampdata.SSP_SAP.RH.contra.snipHbT(:);
pp = p*gcampdata.SSP_SAP.RH.contra.snipGCaMP*100;
plot(gcampdata.SSP_SAP.RH.contra.snipGCaMP*100,pp,'color',colors('north texas green'));
xlabel('GCaMP (%)')
ylabel('\DeltaHbT (\muM)')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([0,16])
ylim([0,40])
axis square

%% Save figure and stats
if saveState == true
    dirpath = [rootFolder delim 'MATLAB Figs & Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig4,[dirpath 'Fig4']);
    diaryFile = [dirpath 'Fig4_Stats.txt'];
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
    
    % LFP delta power (alert)
    disp('======================================================================================================================')
    disp('LFP Alert delta power: Blank (N = 9, 4M/5F); SSP (N = 9, 5M/4F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(data.Blank_SAP.RH.Alert.deltaS,'omitnan')) ' +/- ' num2str(std(data.Blank_SAP.RH.Alert.deltaS,0,1,'omitnan')./sqrt(size(data.Blank_SAP.RH.Alert.deltaS,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(data.SSP_SAP.RH.Alert.deltaS,'omitnan')) ' +/- ' num2str(std(data.SSP_SAP.RH.Alert.deltaS,0,1,'omitnan')./sqrt(size(data.SSP_SAP.RH.Alert.deltaS,1)))]); disp(' ')
    disp('GLME statistics for LFP delta (1:4 Hz)')
    disp(lfpStats.RH.Alert.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % LFP delta power (asleep)
    disp('======================================================================================================================')
    disp('LFP Asleep delta power: Blank (N = 7, 3M/4F); SSP (N = 7, 4M/3F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(data.Blank_SAP.RH.Asleep.deltaS,'omitnan')) ' +/- ' num2str(std(data.Blank_SAP.RH.Asleep.deltaS,0,1,'omitnan')./sqrt(size(data.Blank_SAP.RH.Asleep.deltaS,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(data.SSP_SAP.RH.Asleep.deltaS,'omitnan')) ' +/- ' num2str(std(data.SSP_SAP.RH.Asleep.deltaS,0,1,'omitnan')./sqrt(size(data.SSP_SAP.RH.Asleep.deltaS,1)))]); disp(' ')
    disp('GLME statistics for LFP delta (1:4 Hz)')
    disp(lfpStats.RH.Asleep.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % LFP delta power (all)
    disp('======================================================================================================================')
    disp('LFP All delta power: Blank (N = 9, 4M/5F); SSP (N = 9, 5M/4F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(data.Blank_SAP.RH.All.deltaS,'omitnan')) ' +/- ' num2str(std(data.Blank_SAP.RH.All.deltaS,0,1,'omitnan')./sqrt(size(data.Blank_SAP.RH.All.deltaS,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(data.SSP_SAP.RH.All.deltaS,'omitnan')) ' +/- ' num2str(std(data.SSP_SAP.RH.All.deltaS,0,1,'omitnan')./sqrt(size(data.SSP_SAP.RH.All.deltaS,1)))]); disp(' ')
    disp('GLME statistics for LFP delta (1:4 Hz)')
    disp(lfpStats.RH.All.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % Gamma-HbT
    disp('======================================================================================================================')
    disp(gammaMdl.Coefficients)
    disp(['Blank slope: ' num2str(mean(iosEphysData.Blank_SAP.RH.contra.snipHbT./iosEphysData.Blank_SAP.RH.contra.snipGam,1)) ' +/- ' num2str(std((iosEphysData.Blank_SAP.RH.contra.snipHbT./iosEphysData.Blank_SAP.RH.contra.snipGam),0,1))])
    disp(['SP slope: ' num2str(mean(iosEphysData.SSP_SAP.RH.contra.snipHbT./iosEphysData.SSP_SAP.RH.contra.snipGam,1)) ' +/ ' num2str(std((iosEphysData.SSP_SAP.RH.contra.snipHbT./iosEphysData.SSP_SAP.RH.contra.snipGam),0,1))])

    % GCaMP-HbT
    disp('======================================================================================================================')
    disp(gcampMdl.Coefficients)
    disp(['Blank slope: ' num2str(mean(gcampdata.Blank_SAP.RH.contra.snipHbT./gcampdata.Blank_SAP.RH.contra.snipGCaMP,1)) ' +/- ' num2str(std((gcampdata.Blank_SAP.RH.contra.snipHbT./gcampdata.Blank_SAP.RH.contra.snipGCaMP),0,1))])
    disp(['SP slope: ' num2str(mean(gcampdata.SSP_SAP.RH.contra.snipHbT./gcampdata.SSP_SAP.RH.contra.snipGCaMP,1)) ' +/ ' num2str(std((gcampdata.SSP_SAP.RH.contra.snipHbT./gcampdata.SSP_SAP.RH.contra.snipGCaMP),0,1))])

    diary off
end
