function [Results_PearsonCorr_GCaMP] = AnalyzePearsonCorrelation_GCaMP_eLife2025(animalID,group,set,rootFolder,delim,Results_PearsonCorr_GCaMP)
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
% find and load RestData.mat struct
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
dataTypes = {'HbT','HbO','HbR','GCaMP'};
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% Rest
    clear LH_finalRestData RH_finalRestData fLH_finalRestData fRH_finalRestData LH_ProcRestData RH_ProcRestData fLH_ProcRestData fRH_ProcRestData rest_R rest_fR
    samplingRate = RestData.(dataType).LH.samplingRate;
    [restLogical] = FilterEvents_IOS_eLife2025(RestData.(dataType).LH,RestCriteria);
    [stimLogical] = FilterEvents_IOS_eLife2025(RestData.(dataType).LH,RestStimCriteria);
    combRestLogical = logical(restLogical.*stimLogical);
    restFileIDs = RestData.(dataType).LH.fileIDs(combRestLogical,:);
    restEventTimes = RestData.(dataType).LH.eventTimes(combRestLogical,:);
    restDurations = RestData.(dataType).LH.durations(combRestLogical,:);
    LH_RestingData = RestData.(dataType).LH.data(combRestLogical,:);
    RH_RestingData = RestData.(dataType).RH.data(combRestLogical,:);
    % keep only the data that occurs within the manually-approved alert regions
    [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2025(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2025(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    fLH_RestingData = RestData.(dataType).fLH.data(combRestLogical,:);
    fRH_RestingData = RestData.(dataType).fRH.data(combRestLogical,:);
    % keep only the data that occurs within the manually-approved alert regions
    [fLH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2025(fLH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [fRH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2025(fRH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    % lowpass filter
    [z,p,k] = butter(4,0.5/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % filter, detrend, and truncate data to minimum length to match events
    for bb = 1:length(LH_finalRestData)
        LH_ProcRestData{bb,1} = detrend(filtfilt(sos,g,LH_finalRestData{bb,1}(1:params.minTime.Rest*samplingRate)),'constant');
        RH_ProcRestData{bb,1} = detrend(filtfilt(sos,g,RH_finalRestData{bb,1}(1:params.minTime.Rest*samplingRate)),'constant');
        fLH_ProcRestData{bb,1} = detrend(filtfilt(sos,g,fLH_finalRestData{bb,1}(1:params.minTime.Rest*samplingRate)),'constant');
        fRH_ProcRestData{bb,1} = detrend(filtfilt(sos,g,fRH_finalRestData{bb,1}(1:params.minTime.Rest*samplingRate)),'constant');
    end
    % analyze correlation coefficient of resting epochs
    for bb = 1:length(LH_ProcRestData)
        rest_CC = corrcoef(LH_ProcRestData{bb,1},RH_ProcRestData{bb,1});
        rest_R(bb,1) = rest_CC(2,1);
        rest_fCC = corrcoef(fLH_ProcRestData{bb,1},fRH_ProcRestData{bb,1});
        rest_fR(bb,1) = rest_fCC(2,1);
    end
    % save results
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Rest.R = rest_R;
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Rest.fR = rest_fR;
    %% Whisk
    clear LH_finalWhiskData RH_finalWhiskData fLH_finalWhiskData fRH_finalWhiskData LH_ProcWhiskData RH_ProcWhiskData fLH_ProcWhiskData fRH_ProcWhiskData whisk_R whisk_fR
    [whiskLogical] = FilterEvents_IOS_eLife2025(EventData.(dataType).LH.whisk,WhiskCriteria);
    [stimLogical] = FilterEvents_IOS_eLife2025(EventData.(dataType).LH.whisk,WhiskStimCriteria);
    combWhiskLogical = logical(whiskLogical.*stimLogical);
    whiskFileIDs = EventData.(dataType).LH.whisk.fileIDs(combWhiskLogical,:);
    whiskEventTimes = EventData.(dataType).LH.whisk.eventTime(combWhiskLogical,:);
    whiskDurations = EventData.(dataType).LH.whisk.duration(combWhiskLogical,:);
    LH_whiskData = EventData.(dataType).LH.whisk.data(combWhiskLogical,:);
    RH_whiskData = EventData.(dataType).RH.whisk.data(combWhiskLogical,:);
    % keep only the data that occurs within the manually-approved alert regions
    [LH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS_eLife2025(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS_eLife2025(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    fLH_whiskData = EventData.(dataType).fLH.whisk.data(combWhiskLogical,:);
    fRH_whiskData = EventData.(dataType).fRH.whisk.data(combWhiskLogical,:);
    % keep only the data that occurs within the manually-approved alert regions
    [fLH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS_eLife2025(fLH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    [fRH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS_eLife2025(fRH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % filter, detrend, and take data from whisk onset through 5 seconds
    for bb = 1:size(LH_finalWhiskData,1)
        LH_ProcWhiskData(bb,:) = detrend(filtfilt(sos,g,LH_finalWhiskData(bb,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        RH_ProcWhiskData(bb,:) = detrend(filtfilt(sos,g,RH_finalWhiskData(bb,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        fLH_ProcWhiskData(bb,:) = detrend(filtfilt(sos,g,fLH_finalWhiskData(bb,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        fRH_ProcWhiskData(bb,:) = detrend(filtfilt(sos,g,fRH_finalWhiskData(bb,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
    end
    % analyze correlation coefficient between epochs
    for bb = 1:size(LH_ProcWhiskData,1)
        whisk_CC = corrcoef(LH_ProcWhiskData(bb,:),RH_ProcWhiskData(bb,:));
        whisk_R(bb,1) = whisk_CC(2,1);
        whisk_fCC = corrcoef(fLH_ProcWhiskData(bb,:),fRH_ProcWhiskData(bb,:));
        whisk_fR(bb,1) = whisk_fCC(2,1);
    end
    % save results
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Whisk.R = whisk_R;
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Whisk.fR = whisk_fR;
    %% Alert
    clear LH_AlertData RH_AlertData LH_ProcAlertData RH_ProcAlertData fLH_AlertData fRH_AlertData fLH_ProcAlertData fRH_ProcAlertData alert_R alert_fR scoringLabels
    LH_AlertData = [];
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,~,fileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(fileID,ScoringResults.fileIDs{cc,1}) == true
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
                LH_AlertData{zz,1} = ProcData.data.(dataType).LH(1:300*samplingRate);
                RH_AlertData{zz,1} = ProcData.data.(dataType).RH(1:300*samplingRate);
                fLH_AlertData{zz,1} = ProcData.data.(dataType).fLH(1:300*samplingRate);
                fRH_AlertData{zz,1} = ProcData.data.(dataType).fRH(1:300*samplingRate);
                zz = zz + 1;
            end
        end
        scoringLabelsB = scoringLabels(61:120);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsB,'Not Sleep')) >= 48
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                LH_AlertData{zz,1} = ProcData.data.(dataType).LH(300*samplingRate + 1:600*samplingRate);
                RH_AlertData{zz,1} = ProcData.data.(dataType).RH(300*samplingRate + 1:600*samplingRate);
                fLH_AlertData{zz,1} = ProcData.data.(dataType).fLH(300*samplingRate + 1:600*samplingRate);
                fRH_AlertData{zz,1} = ProcData.data.(dataType).fRH(300*samplingRate + 1:600*samplingRate);
                zz = zz + 1;
            end
        end
        scoringLabelsC = scoringLabels(121:180);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsC,'Not Sleep')) >= 48
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                LH_AlertData{zz,1} = ProcData.data.(dataType).LH(600*samplingRate + 1:end);
                RH_AlertData{zz,1} = ProcData.data.(dataType).RH(600*samplingRate + 1:end);
                fLH_AlertData{zz,1} = ProcData.data.(dataType).fLH(600*samplingRate + 1:end);
                fRH_AlertData{zz,1} = ProcData.data.(dataType).fRH(600*samplingRate + 1:end);
                zz = zz + 1;
            end
        end
    end
    if isempty(LH_AlertData) == false
        % filter and detrend data
        for bb = 1:length(LH_AlertData)
            LH_ProcAlertData{bb,1} = detrend(filtfilt(sos,g,LH_AlertData{bb,1}),'constant');
            RH_ProcAlertData{bb,1} = detrend(filtfilt(sos,g,RH_AlertData{bb,1}),'constant');
            fLH_ProcAlertData{bb,1} = detrend(filtfilt(sos,g,fLH_AlertData{bb,1}),'constant');
            fRH_ProcAlertData{bb,1} = detrend(filtfilt(sos,g,fRH_AlertData{bb,1}),'constant');
        end
        % analyze correlation coefficient between epochs
        for bb = 1:length(LH_ProcAlertData)
            alert_CC = corrcoef(LH_ProcAlertData{bb,1},RH_ProcAlertData{bb,1});
            alert_R(bb,1) = alert_CC(2,1);
            alert_fCC = corrcoef(fLH_ProcAlertData{bb,1},fRH_ProcAlertData{bb,1});
            alert_fR(bb,1) = alert_fCC(2,1);
        end
        % save results
        Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Alert.R = alert_R;
        Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Alert.fR = alert_fR;
    else
        % save results
        Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Alert.R = [];
        Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Alert.fR = [];
    end
    %% Asleep
    clear LH_AsleepData RH_AsleepData LH_ProcAsleepData RH_ProcAsleepData fLH_AsleepData fRH_AsleepData fLH_ProcAsleepData fRH_ProcAsleepData asleep_R asleep_fR scoringLabels
    LH_AsleepData = [];
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,~,fileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(fileID,ScoringResults.fileIDs{cc,1}) == true
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
                LH_AsleepData{zz,1} = ProcData.data.(dataType).LH(1:300*samplingRate);
                RH_AsleepData{zz,1} = ProcData.data.(dataType).RH(1:300*samplingRate);
                fLH_AsleepData{zz,1} = ProcData.data.(dataType).fLH(1:300*samplingRate);
                fRH_AsleepData{zz,1} = ProcData.data.(dataType).fRH(1:300*samplingRate);
                zz = zz + 1;
            end
        end
        scoringLabelsB = scoringLabels(61:120);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsB,'Not Sleep')) <= 12
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                LH_AsleepData{zz,1} = ProcData.data.(dataType).LH(300*samplingRate + 1:600*samplingRate);
                RH_AsleepData{zz,1} = ProcData.data.(dataType).RH(300*samplingRate + 1:600*samplingRate);
                fLH_AsleepData{zz,1} = ProcData.data.(dataType).fLH(300*samplingRate + 1:600*samplingRate);
                fRH_AsleepData{zz,1} = ProcData.data.(dataType).fRH(300*samplingRate + 1:600*samplingRate);
                zz = zz + 1;
            end
        end
        scoringLabelsC = scoringLabels(121:180);
        % check labels to match arousal state
        if sum(strcmp(scoringLabelsC,'Not Sleep')) <= 12
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                LH_AsleepData{zz,1} = ProcData.data.(dataType).LH(600*samplingRate + 1:end);
                RH_AsleepData{zz,1} = ProcData.data.(dataType).RH(600*samplingRate + 1:end);
                fLH_AsleepData{zz,1} = ProcData.data.(dataType).fLH(600*samplingRate + 1:end);
                fRH_AsleepData{zz,1} = ProcData.data.(dataType).fRH(600*samplingRate + 1:end);
                zz = zz + 1;
            end
        end
    end
    if isempty(LH_AsleepData) == false
        % filter and detrend data
        for bb = 1:length(LH_AsleepData)
            LH_ProcAsleepData{bb,1} = detrend(filtfilt(sos,g,LH_AsleepData{bb,1}),'constant');
            RH_ProcAsleepData{bb,1} = detrend(filtfilt(sos,g,RH_AsleepData{bb,1}),'constant');
            fLH_ProcAsleepData{bb,1} = detrend(filtfilt(sos,g,fLH_AsleepData{bb,1}),'constant');
            fRH_ProcAsleepData{bb,1} = detrend(filtfilt(sos,g,fRH_AsleepData{bb,1}),'constant');
        end
        % analyze correlation coefficient between epochs
        for bb = 1:length(LH_ProcAsleepData)
            asleep_CC = corrcoef(LH_ProcAsleepData{bb,1},RH_ProcAsleepData{bb,1});
            asleep_R(bb,1) = asleep_CC(2,1);
            asleep_fCC = corrcoef(fLH_ProcAsleepData{bb,1},fRH_ProcAsleepData{bb,1});
            asleep_fR(bb,1) = asleep_fCC(2,1);
        end
        % save results
        Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Asleep.R = asleep_R;
        Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Asleep.fR = asleep_fR;
    else
        % save results
        Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Asleep.R = [];
        Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).Asleep.fR = [];
    end
    %% All
    clear LH_AllData RH_AllData LH_ProcAllData RH_ProcAllData fLH_AllData fRH_AllData fLH_ProcAllData fRH_ProcAllData all_R all_fR
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        load(procDataFileID,'-mat')
        stims = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(stims) == true
            LH_AllData{zz,1} = ProcData.data.(dataType).LH(1:300*samplingRate);
            RH_AllData{zz,1} = ProcData.data.(dataType).RH(1:300*samplingRate);
            fLH_AllData{zz,1} = ProcData.data.(dataType).fLH(1:300*samplingRate);
            fRH_AllData{zz,1} = ProcData.data.(dataType).fRH(1:300*samplingRate);
            zz = zz + 1;
            LH_AllData{zz,1} = ProcData.data.(dataType).LH(300*samplingRate + 1:600*samplingRate);
            RH_AllData{zz,1} = ProcData.data.(dataType).RH(300*samplingRate + 1:600*samplingRate);
            fLH_AllData{zz,1} = ProcData.data.(dataType).fLH(300*samplingRate + 1:600*samplingRate);
            fRH_AllData{zz,1} = ProcData.data.(dataType).fRH(300*samplingRate + 1:600*samplingRate);
            zz = zz + 1;
            LH_AllData{zz,1} = ProcData.data.(dataType).LH(600*samplingRate + 1:end);
            RH_AllData{zz,1} = ProcData.data.(dataType).RH(600*samplingRate + 1:end);
            fLH_AllData{zz,1} = ProcData.data.(dataType).fLH(600*samplingRate + 1:end);
            fRH_AllData{zz,1} = ProcData.data.(dataType).fRH(600*samplingRate + 1:end);
            zz = zz + 1;
        end
    end
    % filter and detrend data
    for bb = 1:length(LH_AllData)
        LH_ProcAllData{bb,1} = detrend(filtfilt(sos,g,LH_AllData{bb,1}),'constant');
        RH_ProcAllData{bb,1} = detrend(filtfilt(sos,g,RH_AllData{bb,1}),'constant');
        fLH_ProcAllData{bb,1} = detrend(filtfilt(sos,g,fLH_AllData{bb,1}),'constant');
        fRH_ProcAllData{bb,1} = detrend(filtfilt(sos,g,fRH_AllData{bb,1}),'constant');
    end
    % analyze correlation coefficient between epochs
    for bb = 1:length(LH_ProcAllData)
        all_CC = corrcoef(LH_ProcAllData{bb,1},RH_ProcAllData{bb,1});
        all_R(bb,1) = all_CC(2,1);
        all_fCC = corrcoef(fLH_ProcAllData{bb,1},fRH_ProcAllData{bb,1});
        all_fR(bb,1) = all_fCC(2,1);
    end
    % save results
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).All.R = all_R;
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).All.fR = all_fR;
    %% NREM
    clear LH_nremData RH_nremData fLH_nremData fRH_nremData nrem_R nrem_fR
    [LH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [RH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [fLH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).fLH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [fRH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).fRH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % filter, detrend, and truncate data to data to minimum length to match events
    for j = 1:length(LH_nremData)
        LH_nremData{j,1} = detrend(filtfilt(sos,g,LH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{j,1} = detrend(filtfilt(sos,g,RH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        fLH_nremData{j,1} = detrend(filtfilt(sos,g,fLH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        fRH_nremData{j,1} = detrend(filtfilt(sos,g,fRH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
    end
    % analyze correlation coefficient between epochs
    for bb = 1:length(LH_nremData)
        nrem_CC = corrcoef(LH_nremData{bb,1},RH_nremData{bb,1});
        nrem_R(bb,1) = nrem_CC(2,1);
        nrem_fCC = corrcoef(fLH_nremData{bb,1},fRH_nremData{bb,1});
        nrem_fR(bb,1) = nrem_fCC(2,1);
    end
    % save results
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).NREM.R = nrem_R;
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).NREM.fR = nrem_fR;
    %% REM
    clear LH_remData RH_remData fLH_remData fRH_remData rem_R rem_fR
    [LH_remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [RH_remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [fLH_remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).fLH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [fRH_remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).fRH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % filter, detrend, and truncate data to data to minimum length to match events
    for m = 1:length(LH_remData)
        LH_remData{m,1} = detrend(filtfilt(sos,g,LH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        RH_remData{m,1} = detrend(filtfilt(sos,g,RH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        fLH_remData{m,1} = detrend(filtfilt(sos,g,fLH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        fRH_remData{m,1} = detrend(filtfilt(sos,g,fRH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
    end
    % analyze correlation coefficient between epochs
    for bb = 1:length(LH_remData)
        rem_CC = corrcoef(LH_remData{bb,1},RH_remData{bb,1});
        rem_R(bb,1) = rem_CC(2,1);
        rem_fCC = corrcoef(fLH_remData{bb,1},fRH_remData{bb,1});
        rem_fR(bb,1) = rem_fCC(2,1);
    end
    % save results
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).REM.R = rem_R;
    Results_PearsonCorr_GCaMP.(group).(animalID).(dataType).REM.fR = rem_fR;
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PearsonCorr_GCaMP.mat','Results_PearsonCorr_GCaMP')
cd([rootFolder delim 'Data'])