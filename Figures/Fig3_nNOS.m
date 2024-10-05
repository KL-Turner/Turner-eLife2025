function [] = Fig3_nNOS(rootFolder,saveState,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

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

%% Running spectroscopy
path = [rootFolder delim 'Results_Zhang'];
cd(path)
experiment.Blank = {'T192';'T205';'T206';'T208';'T209';'T211';'T225'};
runningBlankGroup(1:length(experiment.Blank)) = {'Blank_SAP'};
experiment.SSP = {'T200';'T212';'T213';'T216';'T217';'T218';'T219'};
runningSSPGroup(1:length(experiment.SSP)) = {'SSP_SAP'};
% set groups based on toxin injection
runningGroup.Blank.HR = [];
runningGroup.Blank.HbT = [];
runningGroup.Blank.HbD = [];
runningGroup.Blank.RestVar_HR = [];
runningGroup.Blank.RestVar_HbT = [];
runningGroup.Blank.RestVar_HbD = [];
runningGroup.SSP = runningGroup.Blank;
% fieldnames
fields = fieldnames(experiment);
for a0 = 1:numel(fields)
    expCondition = fields{a0};
    for a1 = 1:numel(experiment.(expCondition))
        animalID = experiment.(expCondition){a1};
        searchfolder = fullfile(path,'Results',animalID);
        ind_targetFile = getfilenames(searchfolder, [expCondition,'-SAP', '*.mat']);
        if ~isempty(ind_targetFile)
            out.(animalID) = genFigure_individual_SAP(ind_targetFile);
        end
        runningGroup.(expCondition).HR = [runningGroup.(expCondition).HR;nanmean(out.(animalID).HR.LTA,1)]; %#ok<*NANMEAN>
        runningGroup.(expCondition).HbT = [runningGroup.(expCondition).HbT;nanmean(out.(animalID).HbT.PC.LTA,1)];
        runningGroup.(expCondition).HbD = [runningGroup.(expCondition).HbD;nanmean(out.(animalID).HbD.PC.LTA,1)];
        runningGroup.(expCondition).RestVar_HR = [runningGroup.(expCondition).RestVar_HR; nanmean(out.(animalID).RestVar.PC.HR)];
        runningGroup.(expCondition).RestVar_HbT = [runningGroup.(expCondition).RestVar_HbT; nanmean(out.(animalID).RestVar.PC.HbT)];
        runningGroup.(expCondition).RestVar_HbD = [runningGroup.(expCondition).RestVar_HbD; nanmean(out.(animalID).RestVar.PC.HbD)];
    end
end
time = -89/30:1/30:5;
blankRunningMean = mean(runningGroup.Blank.HbT,1);
sapRunningMean = mean(runningGroup.SSP.HbT,1);
blankRunningStdErr = std(runningGroup.Blank.HbT,[],1)/sqrt(size(runningGroup.Blank.HbT,1));
sapRunningStdErr = std(runningGroup.SSP.HbT,[],1)/sqrt(size(runningGroup.SSP.HbT,1));
qzGroups = {'Blank','SSP'};
for aa = 1:length(qzGroups)
    qzGroup = qzGroups{1,aa};
    for ee = 1:size(runningGroup.(qzGroup).HbT,1)
        startIdx = find(time == 1.5);
        endIdx =  find(time == 5);
        hbtSnip = runningGroup.(qzGroup).HbT(ee,:);
        runningStats.(qzGroup).AUC(ee,1) = mean(hbtSnip(startIdx:endIdx));
    end
end
% statistics - generalized linear mixed effects model
runningStats.tableSize = cat(1,runningStats.Blank.AUC,runningStats.SSP.AUC);
runningStats.Table = table('Size',[size(runningStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
runningStats.Table.Mouse = cat(1,experiment.Blank,experiment.SSP);
runningStats.Table.Group = cat(1,runningBlankGroup',runningSSPGroup');
runningStats.Table.AUC = cat(1,runningStats.Blank.AUC,runningStats.SSP.AUC);
runningStats.FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
runningStats.Stats = fitglme(runningStats.Table,runningStats.FitFormula);

%% Two photon long whisker stim
path = [rootFolder delim 'Results_Turner'];
cd(path)
resultsStruct = 'Results_Evoked_2P';
load(resultsStruct);
% loop variables
groups = {'SSP_SAP','Blank_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
dataTypes = {'diameter','baseline','timeVector','count','group','animalID'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_2P.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Evoked_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            for dd = 1:length(solenoids)
                solenoid = solenoids{1,dd};
                twoPdata.(group).(solenoid).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(twoPdata.(group).(solenoid),dataType) == false
                        if any(strcmp(dataType,{'group','animalID'})) == true
                            twoPdata.(group).(solenoid).(dataType) = {};
                        else
                            twoPdata.(group).(solenoid).(dataType) = [];
                        end
                    end
                end
                twoPdata.(group).(solenoid).diameter = cat(1,twoPdata.(group).(solenoid).diameter,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).diameter*100);
                twoPdata.(group).(solenoid).baseline = cat(1,twoPdata.(group).(solenoid).baseline,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).baseline);
                twoPdata.(group).(solenoid).count = cat(1,twoPdata.(group).(solenoid).count,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).count);
                twoPdata.(group).(solenoid).timeVector = cat(1,twoPdata.(group).(solenoid).timeVector,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).timeVector);
                twoPdata.(group).(solenoid).group = cat(1,twoPdata.(group).(solenoid).group,group);
                twoPdata.(group).(solenoid).animalID = cat(1,twoPdata.(group).(solenoid).animalID,animalID);
            end
        end
    end
end
% pair stimulation types with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(solenoids)
        solenoid = solenoids{1,bb};
        [comparison] = FindSolenoidComparison('RH',solenoid);
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            twoPdata.(group).(comparison).(dataType) = twoPdata.(group).(solenoid).(dataType);
            if any(strcmp(dataType,{'group','animalID'})) == false
                twoPdata.(group).(comparison).(['mean_'  dataType]) = mean(twoPdata.(group).(solenoid).(dataType),1);
                twoPdata.(group).(comparison).(['stdErr_' dataType]) = std(twoPdata.(group).(solenoid).(dataType),0,1)./sqrt(size(twoPdata.(group).(solenoid).(dataType),1));
            end
        end
    end
end
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:size(twoPdata.(group).contra.diameter,1)
        startIdx = find(twoPdata.(group).contra.timeVector(ee,:) == 3);
        endIdx =  find(twoPdata.(group).contra.timeVector(ee,:) == 7);
        diameterSnip = twoPdata.(group).contra.diameter(ee,:);
        twoPdata.(group).contra.AUC(ee,1) = mean(diameterSnip(startIdx:endIdx));
    end
end
% statistics - generalized linear mixed effects model
diameterStats.tableSize = cat(1,twoPdata.Blank_SAP.contra.AUC,twoPdata.SSP_SAP.contra.AUC);
diameterStats.Table = table('Size',[size(diameterStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
diameterStats.Table.Mouse = cat(1,twoPdata.Blank_SAP.contra.animalID,twoPdata.SSP_SAP.contra.animalID);
diameterStats.Table.Group = cat(1,twoPdata.Blank_SAP.contra.group,twoPdata.SSP_SAP.contra.group);
diameterStats.Table.AUC = cat(1,twoPdata.Blank_SAP.contra.AUC,twoPdata.SSP_SAP.contra.AUC);
diameterStats.FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
diameterStats.Stats = fitglme(diameterStats.Table,diameterStats.FitFormula);

%% IOS long whisker stim
path = [rootFolder delim 'Results_Turner'];
cd(path)
resultsStruct = 'Results_Evoked_Pulse';
load(resultsStruct);
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
dataTypes = {'HbT','timeVector','animalID','group'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_Pulse.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            pulseData.(group).(solenoid).dummCheck = 1;
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if isfield(pulseData.(group).(solenoid),dataType) == false
                    if any(strcmp(dataType,{'group','animalID'})) == true
                        pulseData.(group).(solenoid).(dataType) = {};
                    else
                        pulseData.(group).(solenoid).(dataType) = [];
                    end
                end
            end
            pulseData.(group).(solenoid).HbT = cat(1,pulseData.(group).(solenoid).HbT,Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).HbT);
            pulseData.(group).(solenoid).timeVector = cat(1,pulseData.(group).(solenoid).timeVector,Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).timeVector);
            pulseData.(group).(solenoid).group = cat(1,pulseData.(group).(solenoid).group,group);
            pulseData.(group).(solenoid).animalID = cat(1,pulseData.(group).(solenoid).animalID,animalID);
        end
    end
end
% pair stimulation types with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(solenoids)
        solenoid = solenoids{1,bb};
        [comparison] = FindSolenoidComparison('RH',solenoid);
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            pulseData.(group).(comparison).(dataType) = pulseData.(group).(solenoid).(dataType);
            if any(strcmp(dataType,{'group','animalID'})) == false
                pulseData.(group).(comparison).(['mean_'  dataType]) = mean(pulseData.(group).(solenoid).(dataType),1);
                pulseData.(group).(comparison).(['stdErr_' dataType]) = std(pulseData.(group).(solenoid).(dataType),1)./sqrt(size(pulseData.(group).(solenoid).(dataType),1));
            end
        end
    end
end
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:size(pulseData.(group).contra.HbT,1)
        startIdx = find(pulseData.(group).contra.timeVector(ee,:) == 1.5);
        endIdx =  find(pulseData.(group).contra.timeVector(ee,:) == 6.5);
        pulseSnip = pulseData.(group).contra.HbT(ee,:);
        pulseData.(group).contra.AUC(ee,1) = mean(pulseSnip(startIdx:endIdx));
    end
end
% statistics - generalized linear mixed effects model
pulseStats.tableSize = cat(1,pulseData.Blank_SAP.contra.AUC,pulseData.SSP_SAP.contra.AUC);
pulseStats.Table = table('Size',[size(pulseStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
pulseStats.Table.Mouse = cat(1,pulseData.Blank_SAP.contra.animalID,pulseData.SSP_SAP.contra.animalID);
pulseStats.Table.Group = cat(1,pulseData.Blank_SAP.contra.group,pulseData.SSP_SAP.contra.group);
pulseStats.Table.AUC = cat(1,pulseData.Blank_SAP.contra.AUC,pulseData.SSP_SAP.contra.AUC);
pulseStats.FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
pulseStats.Stats = fitglme(pulseStats.Table,pulseStats.FitFormula);

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

%% Figure 3
Fig3 = figure('Name','Fig. 3');

% Ephys stimulation
subplot(2,3,1)
p1 = plot(iosEphysData.Blank_SAP.RH.contra.mean_timeVector,iosEphysData.Blank_SAP.RH.contra.mean_HbT,'color',colors('black'),'LineWidth',2);
hold on
plot(iosEphysData.Blank_SAP.RH.contra.mean_timeVector,iosEphysData.Blank_SAP.RH.contra.mean_HbT + iosEphysData.Blank_SAP.RH.contra.stdErr_HbT,'color',colors('black'),'LineWidth',0.25)
plot(iosEphysData.Blank_SAP.RH.contra.mean_timeVector,iosEphysData.Blank_SAP.RH.contra.mean_HbT - iosEphysData.Blank_SAP.RH.contra.stdErr_HbT,'color',colors('black'),'LineWidth',0.25)
p2 = plot(iosEphysData.SSP_SAP.RH.contra.mean_timeVector,iosEphysData.SSP_SAP.RH.contra.mean_HbT,'color',colors('electric purple'),'LineWidth',2);
plot(iosEphysData.SSP_SAP.RH.contra.mean_timeVector,iosEphysData.SSP_SAP.RH.contra.mean_HbT + iosEphysData.SSP_SAP.RH.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
plot(iosEphysData.SSP_SAP.RH.contra.mean_timeVector,iosEphysData.SSP_SAP.RH.contra.mean_HbT - iosEphysData.SSP_SAP.RH.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
axis square
axis tight
xlim([-2,5]);

% GCaMP HbT
subplot(2,3,2);
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbT,'color',colors('black'),'LineWidth',2);
hold on;
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbT + gcampdata.Blank_SAP.RH.contra.stdErr_HbT,'color',colors('black'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbT - gcampdata.Blank_SAP.RH.contra.stdErr_HbT,'color',colors('black'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbT,'color',colors('electric purple'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbT + gcampdata.SSP_SAP.RH.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbT - gcampdata.SSP_SAP.RH.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,10]);

% GCaMP F/F
subplot(2,3,3)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_GCaMP,'color',colors('black'),'LineWidth',2);
hold on;
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_GCaMP + gcampdata.Blank_SAP.RH.contra.stdErr_GCaMP,'color',colors('black'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_GCaMP - gcampdata.Blank_SAP.RH.contra.stdErr_GCaMP,'color',colors('black'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_GCaMP,'color',colors('north texas green'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_GCaMP + gcampdata.SSP_SAP.RH.contra.stdErr_GCaMP,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_GCaMP - gcampdata.SSP_SAP.RH.contra.stdErr_GCaMP,'color',colors('north texas green'),'LineWidth',0.25)
ylabel('\DeltaF/F (%)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,10]);

% GCaMP HbO / HbR
subplot(2,3,4);
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbO,'color',colors('black'),'LineWidth',2);
hold on;
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbO + gcampdata.Blank_SAP.RH.contra.stdErr_HbO,'color',colors('black'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbO - gcampdata.Blank_SAP.RH.contra.stdErr_HbO,'color',colors('black'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbO,'color',colors('deep carrot orange'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbO + gcampdata.SSP_SAP.RH.contra.stdErr_HbO,'color',colors('deep carrot orange'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbO - gcampdata.SSP_SAP.RH.contra.stdErr_HbO,'color',colors('deep carrot orange'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbR,'color',colors('black'),'LineWidth',2);
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbR + gcampdata.Blank_SAP.RH.contra.stdErr_HbR,'color',colors('black'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbR - gcampdata.Blank_SAP.RH.contra.stdErr_HbR,'color',colors('black'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbR,'color',colors('vegas gold'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbR + gcampdata.SSP_SAP.RH.contra.stdErr_HbR,'color',colors('vegas gold'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbR - gcampdata.SSP_SAP.RH.contra.stdErr_HbR,'color',colors('vegas gold'),'LineWidth',0.25)
ylabel('\DeltaHbO-R (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,10]);

% 2p stimulation
subplot(2,3,5)
plot(twoPdata.Blank_SAP.contra.mean_timeVector,twoPdata.Blank_SAP.contra.mean_diameter,'color',colors('black'),'LineWidth',2);
hold on;
plot(twoPdata.Blank_SAP.contra.mean_timeVector,twoPdata.Blank_SAP.contra.mean_diameter + twoPdata.Blank_SAP.contra.stdErr_diameter,'color',colors('black'),'LineWidth',0.25)
plot(twoPdata.Blank_SAP.contra.mean_timeVector,twoPdata.Blank_SAP.contra.mean_diameter - twoPdata.Blank_SAP.contra.stdErr_diameter,'color',colors('black'),'LineWidth',0.25)
plot(twoPdata.SSP_SAP.contra.mean_timeVector,twoPdata.SSP_SAP.contra.mean_diameter,'color',colors('dark candy apple red'),'LineWidth',2);
plot(twoPdata.SSP_SAP.contra.mean_timeVector,twoPdata.SSP_SAP.contra.mean_diameter + twoPdata.SSP_SAP.contra.stdErr_diameter,'color',colors('dark candy apple red'),'LineWidth',0.25)
plot(twoPdata.SSP_SAP.contra.mean_timeVector,twoPdata.SSP_SAP.contra.mean_diameter - twoPdata.SSP_SAP.contra.stdErr_diameter,'color',colors('dark candy apple red'),'LineWidth',0.25)
ylabel('\DeltaD/D (%)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,10]);

%% Save figure and stats
if saveState == true
    dirpath = [rootFolder delim 'MATLAB Figs/Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig3,[dirpath 'Fig3']);
    diaryFile = [dirpath 'Fig3_Stats.txt'];
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
   
    % IOS brief stim AUC
    disp('======================================================================================================================')
    disp('AUC 0.1 s stim: Blank (N = 9, 4M/5F); SSP (N = 9, 5M/4F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(iosEphysData.Blank_SAP.RH.contra.AUC)) ' +/- ' num2str(std(iosEphysData.Blank_SAP.RH.contra.AUC,0,1)./sqrt(size(iosEphysData.Blank_SAP.RH.contra.AUC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(iosEphysData.SSP_SAP.RH.contra.AUC)) ' +/- ' num2str(std(iosEphysData.SSP_SAP.RH.contra.AUC,0,1)./sqrt(size(iosEphysData.SSP_SAP.RH.contra.AUC,1)))]); disp(' ')
    disp('GLME statistics for IOS AUC (t = 2:4 sec)')
    disp(iosEphysStats.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % IOS long stim HbT
    disp('======================================================================================================================')
    disp('IOS 5s stim HbT: Blank (N = 7, 3M/4F); SSP (N = 8, 4M/4F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampdata.Blank_SAP.RH.contra.AUC_HbT)) ' +/- ' num2str(std(gcampdata.Blank_SAP.RH.contra.AUC_HbT,0,1)./sqrt(size(gcampdata.Blank_SAP.RH.contra.AUC_HbT,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampdata.SSP_SAP.RH.contra.AUC_HbT)) ' +/- ' num2str(std(gcampdata.SSP_SAP.RH.contra.AUC_HbT,0,1)./sqrt(size(gcampdata.SSP_SAP.RH.contra.AUC_HbT,1)))]); disp(' ')
    disp('GLME statistics for GCaMP (HbT) (t = 1.5:6.5 sec)')
    disp(gcampStats.HbT.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % IOS long stim GCaMP
    disp('======================================================================================================================')
    disp('IOS 5s stim GCaMP: Blank (N = 7, 3M/4F); SSP (N = 8, 4M/4F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampdata.Blank_SAP.RH.contra.AUC_GCaMP)) ' +/- ' num2str(std(gcampdata.Blank_SAP.RH.contra.AUC_GCaMP,0,1)./sqrt(size(gcampdata.Blank_SAP.RH.contra.AUC_GCaMP,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampdata.SSP_SAP.RH.contra.AUC_GCaMP)) ' +/- ' num2str(std(gcampdata.SSP_SAP.RH.contra.AUC_GCaMP,0,1)./sqrt(size(gcampdata.SSP_SAP.RH.contra.AUC_GCaMP,1)))]); disp(' ')
    disp('GLME statistics for GCaMP (GCaMP) (t = 2:5 sec)')
    disp(gcampStats.GCaMP.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % IOS long stim HbO
    disp('======================================================================================================================')
    disp('IOS 5s stim HbO: Blank (N = 7, 3M/4F); SSP (N = 8, 4M/4F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampdata.Blank_SAP.RH.contra.AUC_HbO)) ' +/- ' num2str(std(gcampdata.Blank_SAP.RH.contra.AUC_HbO,0,1)./sqrt(size(gcampdata.Blank_SAP.RH.contra.AUC_HbO,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampdata.SSP_SAP.RH.contra.AUC_HbO)) ' +/- ' num2str(std(gcampdata.SSP_SAP.RH.contra.AUC_HbO,0,1)./sqrt(size(gcampdata.SSP_SAP.RH.contra.AUC_HbO,1)))]); disp(' ')
    disp('GLME statistics for GCaMP (HbO) (t = 1.5:6.5 sec)')
    disp(gcampStats.HbO.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % IOS long stim HbR
    disp('======================================================================================================================')
    disp('IOS 5s stim HbR: Blank (N = 7, 3M/4F); SSP (N = 8, 4M/4F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampdata.Blank_SAP.RH.contra.AUC_HbR)) ' +/- ' num2str(std(gcampdata.Blank_SAP.RH.contra.AUC_HbR,0,1)./sqrt(size(gcampdata.Blank_SAP.RH.contra.AUC_HbR,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampdata.SSP_SAP.RH.contra.AUC_HbR)) ' +/- ' num2str(std(gcampdata.SSP_SAP.RH.contra.AUC_HbR,0,1)./sqrt(size(gcampdata.SSP_SAP.RH.contra.AUC_HbR,1)))]); disp(' ')
    disp('GLME statistics for GCaMP (HbR) (t = 1.5:6.5 sec)')
    disp(gcampStats.HbR.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % IOS long stim (pulse non-GCaMP)
    disp('======================================================================================================================')
    disp('IOS 5s stim AUC: Blank (N = 8, 5M/3F); SSP (N = 7, 2M/5F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(pulseData.Blank_SAP.contra.AUC)) ' +/- ' num2str(std(pulseData.Blank_SAP.contra.AUC,0,1)./sqrt(size(pulseData.Blank_SAP.contra.AUC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(pulseData.SSP_SAP.contra.AUC)) ' +/- ' num2str(std(pulseData.SSP_SAP.contra.AUC,0,1)./sqrt(size(pulseData.SSP_SAP.contra.AUC,1)))]); disp(' ')
    disp('GLME statistics for IOS pulse AUC (t = 1.5:6.5 sec)')
    disp(pulseStats.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % IOS running
    disp('======================================================================================================================')
    disp('Running AUC: Blank (N = 7, 5M/2F); SSP (N = 7, 2M/5F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(runningStats.Blank.AUC)) ' +/- ' num2str(std(runningStats.Blank.AUC,0,1)./sqrt(size(runningStats.Blank.AUC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(runningStats.SSP.AUC)) ' +/- ' num2str(std(runningStats.SSP.AUC,0,1)./sqrt(size(runningStats.SSP.AUC,1)))]); disp(' ')
    disp('GLME statistics for running AUC (t = 1.5:2.5 sec)')
    disp(runningStats.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % 2P long stim
    disp('======================================================================================================================')
    disp('2P 5s stim AUC: Blank (N = 9,  5M/4F, n = 81); SSP (N = 7, 2M/5F, n = 70); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(twoPdata.Blank_SAP.contra.AUC)) ' +/- ' num2str(std(twoPdata.Blank_SAP.contra.AUC,0,1)./sqrt(size(twoPdata.Blank_SAP.contra.AUC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(twoPdata.SSP_SAP.contra.AUC)) ' +/- ' num2str(std(twoPdata.SSP_SAP.contra.AUC,0,1)./sqrt(size(twoPdata.SSP_SAP.contra.AUC,1)))]); disp(' ')
    disp('GLME statistics for two photon AUC (t = 3:7 sec)')
    disp(diameterStats.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    diary off
end
