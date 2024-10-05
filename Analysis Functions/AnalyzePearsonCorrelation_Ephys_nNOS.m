function [Results_PearsonCorr_Ephys] = AnalyzePearsonCorrelation_Ephys(animalID,group,set,rootFolder,delim,Results_PearsonCorr_Ephys)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% character list of ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% find and load RestData struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% find and load ManualDecisions struct
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% find and load EventData struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
% find and load RestingBaselines strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData strut
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% find and load Forest_ScoringResults struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% parameters
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.Whisk = 7;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% criteria for whisking
WhiskCriteria.Fieldname = {'duration','puffDistance'};
WhiskCriteria.Comparison = {'gt','gt'};
WhiskCriteria.Value = {5,5};
WhiskStimCriteria.Fieldname = {'puffDistance'};
WhiskStimCriteria.Comparison = {'gt'};
WhiskStimCriteria.Value = {5};
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestStimCriteria.Fieldname = {'stimDistances'};
RestStimCriteria.Comparison = {'gt'};
RestStimCriteria.Value = {5};
% loop variables
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
for a = 1:length(dataTypes)
    dataType = dataTypes{1,a};
    %% Rest
    clear LH_finalRestData RH_finalRestData LH_ProcRestData RH_ProcRestData rest_R
    if strcmp(dataType,'HbT') == true
        samplingRate = RestData.(dataType).LH.samplingRate;
        [restLogical] = FilterEvents_IOS(RestData.(dataType).LH,RestCriteria);
        [stimLogical] = FilterEvents_IOS(RestData.(dataType).LH,RestStimCriteria);
        combRestLogical = logical(restLogical.*stimLogical);
        restFileIDs = RestData.(dataType).LH.fileIDs(combRestLogical,:);
        restEventTimes = RestData.(dataType).LH.eventTimes(combRestLogical,:);
        restDurations = RestData.(dataType).LH.durations(combRestLogical,:);
        LH_RestingData = RestData.(dataType).LH.data(combRestLogical,:);
        RH_RestingData = RestData.(dataType).RH.data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    else
        samplingRate = RestData.cortical_LH.(dataType).samplingRate;
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [stimLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestStimCriteria);
        combRestLogical = logical(restLogical.*stimLogical);
        restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
        LH_RestingData = RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
        RH_RestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    end
    % lowpass filter
    [z,p,k] = butter(4,0.5/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % filter, detrend, and truncate data to minimum length to match events
    for gg = 1:length(LH_finalRestData)
        LH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,LH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
        RH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,RH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
    end
    % analyze correlation coefficient of resting epochs
    for n = 1:length(LH_ProcRestData)
        rest_CC = corrcoef(LH_ProcRestData{n,1},RH_ProcRestData{n,1});
        rest_R(n,1) = rest_CC(2,1);
    end
    % save results
    Results_PearsonCorr_Ephys.(group).(animalID).(dataType).Rest.R = rest_R;
    %% Whisk
    clear LH_finalWhiskData RH_finalWhiskData LH_ProcWhiskData RH_ProcWhiskData whisk_R
    if strcmp(dataType,'HbT') == true
        [whiskLogical] = FilterEvents_IOS(EventData.(dataType).LH.whisk,WhiskCriteria);
        [stimLogical] = FilterEvents_IOS(EventData.(dataType).LH.whisk,WhiskStimCriteria);
        combWhiskLogical = logical(whiskLogical.*stimLogical);
        whiskFileIDs = EventData.(dataType).LH.whisk.fileIDs(combWhiskLogical,:);
        whiskEventTimes = EventData.(dataType).LH.whisk.eventTime(combWhiskLogical,:);
        whiskDurations = EventData.(dataType).LH.whisk.duration(combWhiskLogical,:);
        LH_whiskData = EventData.(dataType).LH.whisk.data(combWhiskLogical,:);
        RH_whiskData = EventData.(dataType).RH.whisk.data(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    else
        [whiskLogical] = FilterEvents_IOS(EventData.cortical_LH.(dataType).whisk,WhiskCriteria);
        [stimLogical] = FilterEvents_IOS(EventData.cortical_LH.(dataType).whisk,WhiskStimCriteria);
        combWhiskLogical = logical(whiskLogical.*stimLogical);
        whiskFileIDs = EventData.cortical_LH.(dataType).whisk.fileIDs(combWhiskLogical,:);
        whiskEventTimes = EventData.cortical_LH.(dataType).whisk.eventTime(combWhiskLogical,:);
        whiskDurations = EventData.cortical_LH.(dataType).whisk.duration(combWhiskLogical,:);
        LH_whiskData = EventData.cortical_LH.(dataType).whisk.NormData(combWhiskLogical,:);
        RH_whiskData = EventData.cortical_RH.(dataType).whisk.NormData(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [LH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    end
    % filter, detrend, and take data from whisk onset through 5 seconds
    for gg = 1:size(LH_finalWhiskData,1)
        LH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,LH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        RH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,RH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
    end
    % analyze correlation coefficient between epochs
    for n = 1:size(LH_ProcWhiskData,1)
        whisk_CC = corrcoef(LH_ProcWhiskData(n,:),RH_ProcWhiskData(n,:));
        whisk_R(n,1) = whisk_CC(2,1);
    end
    % save results
    Results_PearsonCorr_Ephys.(group).(animalID).(dataType).Whisk.R = whisk_R;
    %% Alert
    clear LH_AlertData RH_AlertData LH_ProcAlertData RH_ProcAlertData scoringLabels alert_R
    LH_AlertData = [];
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        scoringLabelsA = scoringLabels(1:60);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsA,'Not Sleep')) >= 48
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                if strcmp(dataType,'HbT') == true
                    LH_AlertData{zz,1} = ProcData.data.(dataType).LH(1:300*samplingRate);
                    RH_AlertData{zz,1} = ProcData.data.(dataType).RH(1:300*samplingRate);
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AlertData{zz,1} = (ProcData.data.cortical_LH.(dataType)(1:300*samplingRate) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                        RH_AlertData{zz,1} = (ProcData.data.cortical_RH.(dataType)(1:300*samplingRate) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
        end
        scoringLabelsB = scoringLabels(61:120);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsB,'Not Sleep')) >= 48
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                if strcmp(dataType,'HbT') == true
                    LH_AlertData{zz,1} = ProcData.data.(dataType).LH(300*samplingRate + 1:600*samplingRate);
                    RH_AlertData{zz,1} = ProcData.data.(dataType).RH(300*samplingRate + 1:600*samplingRate);
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AlertData{zz,1} = (ProcData.data.cortical_LH.(dataType)(300*samplingRate + 1:600*samplingRate) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                        RH_AlertData{zz,1} = (ProcData.data.cortical_RH.(dataType)(300*samplingRate + 1:600*samplingRate) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
        end
        scoringLabelsC = scoringLabels(121:180);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsC,'Not Sleep')) >= 48
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                if strcmp(dataType,'HbT') == true
                    LH_AlertData{zz,1} = ProcData.data.(dataType).LH(600*samplingRate + 1:end);
                    RH_AlertData{zz,1} = ProcData.data.(dataType).RH(600*samplingRate + 1:end);
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AlertData{zz,1} = (ProcData.data.cortical_LH.(dataType)(600*samplingRate + 1:end) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                        RH_AlertData{zz,1} = (ProcData.data.cortical_RH.(dataType)(600*samplingRate + 1:end) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
        end
    end
    if isempty(LH_AlertData) == false
        % filter and detrend data
        for gg = 1:length(LH_AlertData)
            LH_ProcAlertData{gg,1} = detrend(filtfilt(sos,g,LH_AlertData{gg,1}),'constant');
            RH_ProcAlertData{gg,1} = detrend(filtfilt(sos,g,RH_AlertData{gg,1}),'constant');
        end
        % analyze correlation coefficient between epochs
        for n = 1:length(LH_ProcAlertData)
            alert_CC = corrcoef(LH_ProcAlertData{n,1},RH_ProcAlertData{n,1});
            alert_R(n,1) = alert_CC(2,1);
        end
        % save results
        Results_PearsonCorr_Ephys.(group).(animalID).(dataType).Alert.R = alert_R;
    else
        % save results
        Results_PearsonCorr_Ephys.(group).(animalID).(dataType).Alert.R = [];
    end
    %% Asleep
    clear LH_AsleepData RH_AsleepData LH_ProcAsleepData RH_ProcAsleepData scoringLabels asleep_R
    LH_AsleepData = [];
    zz= 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        scoringLabelsA = scoringLabels(1:60);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsA,'Not Sleep')) <= 12
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                if strcmp(dataType,'HbT') == true
                    LH_AsleepData{zz,1} = ProcData.data.(dataType).LH(1:300*samplingRate);
                    RH_AsleepData{zz,1} = ProcData.data.(dataType).RH(1:300*samplingRate);
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AsleepData{zz,1} = (ProcData.data.cortical_LH.(dataType)(1:300*samplingRate) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                        RH_AsleepData{zz,1} = (ProcData.data.cortical_RH.(dataType)(1:300*samplingRate) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
        end
        scoringLabelsB = scoringLabels(61:120);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsB,'Not Sleep')) <= 12
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                if strcmp(dataType,'HbT') == true
                    LH_AsleepData{zz,1} = ProcData.data.(dataType).LH(300*samplingRate + 1:600*samplingRate);
                    RH_AsleepData{zz,1} = ProcData.data.(dataType).RH(300*samplingRate + 1:600*samplingRate);
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AsleepData{zz,1} = (ProcData.data.cortical_LH.(dataType)(300*samplingRate + 1:600*samplingRate) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                        RH_AsleepData{zz,1} = (ProcData.data.cortical_RH.(dataType)(300*samplingRate + 1:600*samplingRate) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
        end
        scoringLabelsC = scoringLabels(121:180);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsC,'Not Sleep')) <= 12
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                if strcmp(dataType,'HbT') == true
                    LH_AsleepData{zz,1} = ProcData.data.(dataType).LH(600*samplingRate + 1:end);
                    RH_AsleepData{zz,1} = ProcData.data.(dataType).RH(600*samplingRate + 1:end);
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        LH_AsleepData{zz,1} = (ProcData.data.cortical_LH.(dataType)(600*samplingRate + 1:end) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                        RH_AsleepData{zz,1} = (ProcData.data.cortical_RH.(dataType)(600*samplingRate + 1:end) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
        end
    end
    if isempty(LH_AsleepData) == false
        % filter and detrend data
        for gg = 1:length(LH_AsleepData)
            LH_ProcAsleepData{gg,1} = detrend(filtfilt(sos,g,LH_AsleepData{gg,1}),'constant');
            RH_ProcAsleepData{gg,1} = detrend(filtfilt(sos,g,RH_AsleepData{gg,1}),'constant');
        end
        % analyze correlation coefficient between epochs
        for n = 1:length(LH_ProcAsleepData)
            asleep_CC = corrcoef(LH_ProcAsleepData{n,1},RH_ProcAsleepData{n,1});
            asleep_R(n,1) = asleep_CC(2,1);
        end
        % save results
        Results_PearsonCorr_Ephys.(group).(animalID).(dataType).Asleep.R = asleep_R;
    else
        % save results
        Results_PearsonCorr_Ephys.(group).(animalID).(dataType).Asleep.R = [];
    end
    %% All
    clear LH_AllData RH_AllData LH_ProcAllData RH_ProcAllData all_R
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,~] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        load(procDataFileID,'-mat')
        stims = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(stims) == true
            if strcmp(dataType,'HbT') == true
                LH_AllData{zz,1} = ProcData.data.(dataType).LH(1:300*samplingRate);
                RH_AllData{zz,1} = ProcData.data.(dataType).RH(1:300*samplingRate);
                zz = zz + 1;
                LH_AllData{zz,1} = ProcData.data.(dataType).LH(300*samplingRate + 1:600*samplingRate);
                RH_AllData{zz,1} = ProcData.data.(dataType).RH(300*samplingRate + 1:600*samplingRate);
                zz = zz + 1;
                LH_AllData{zz,1} = ProcData.data.(dataType).LH(600*samplingRate + 1:end);
                RH_AllData{zz,1} = ProcData.data.(dataType).RH(600*samplingRate + 1:end);
                zz = zz + 1;
            else
                motionArtifact = ProcData.notes.motionArtifact;
                if motionArtifact == false
                    LH_AllData{zz,1} = (ProcData.data.cortical_LH.(dataType)(1:300*samplingRate) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                    RH_AllData{zz,1} = (ProcData.data.cortical_RH.(dataType)(1:300*samplingRate) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                    zz = zz + 1;
                    LH_AllData{zz,1} = (ProcData.data.cortical_LH.(dataType)(300*samplingRate + 1:600*samplingRate) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                    RH_AllData{zz,1} = (ProcData.data.cortical_RH.(dataType)(300*samplingRate + 1:600*samplingRate) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                    zz = zz + 1;
                    LH_AllData{zz,1} = (ProcData.data.cortical_LH.(dataType)(600*samplingRate + 1:end) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                    RH_AllData{zz,1} = (ProcData.data.cortical_RH.(dataType)(600*samplingRate + 1:end) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                    zz = zz + 1;
                end
            end
        end
    end
    % filter and detrend data
    for gg = 1:length(LH_AllData)
        LH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,LH_AllData{gg,1}),'constant');
        RH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,RH_AllData{gg,1}),'constant');
    end
    % analyze correlation coefficient between epochs
    for n = 1:length(LH_ProcAllData)
        all_CC = corrcoef(LH_ProcAllData{n,1},RH_ProcAllData{n,1});
        all_R(n,1) = all_CC(2,1);
    end
    % save results
    Results_PearsonCorr_Ephys.(group).(animalID).(dataType).All.R = all_R;
    %% NREM
    clear LH_nremData RH_nremData nrem_R
    if strcmp(dataType,'HbT') == true
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    else
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    end
    % filter, detrend, and truncate data to data to minimum length to match events
    for j = 1:length(LH_nremData)
        LH_nremData{j,1} = detrend(filtfilt(sos,g,LH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{j,1} = detrend(filtfilt(sos,g,RH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
    end
    % analyze correlation coefficient between epochs
    for n = 1:length(LH_nremData)
        nrem_CC = corrcoef(LH_nremData{n,1},RH_nremData{n,1});
        nrem_R(n,1) = nrem_CC(2,1);
    end
    % save results
    Results_PearsonCorr_Ephys.(group).(animalID).(dataType).NREM.R = nrem_R;
    %% REM
    clear LH_remData RH_remData rem_R
    if strcmp(dataType,'HbT') == true
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    else
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    end
    % filter, detrend, and truncate data to data to minimum length to match events
    for m = 1:length(LH_remData)
        LH_remData{m,1} = detrend(filtfilt(sos,g,LH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        RH_remData{m,1} = detrend(filtfilt(sos,g,RH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
    end
    % analyze correlation coefficient between epochs
    for n = 1:length(LH_remData)
        rem_CC = corrcoef(LH_remData{n,1},RH_remData{n,1});
        rem_R(n,1) = rem_CC(2,1);
    end
    % save results
    Results_PearsonCorr_Ephys.(group).(animalID).(dataType).REM.R = rem_R;
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PearsonCorr_Ephys.mat','Results_PearsonCorr_Ephys')
cd([rootFolder delim 'Data'])