function [Results_IntSig_GCaMP] = AnalyzeIntrinsicSignals_GCaMP(animalID,group,set,rootFolder,delim,Results_IntSig_GCaMP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
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
% parameters
modelType = 'Forest';
params.minTime.Rest = 10;
params.Offset = 2;
params.minTime.Whisk = params.Offset + 5;
params.minTime.Stim = params.Offset + 2;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% criteria for the whisking
WhiskCriteria.Fieldname = {'duration','puffDistance'};
WhiskCriteria.Comparison = {'gt','gt'};
WhiskCriteria.Value = {5,5};
WhiskStimCriteria.Fieldname = {'puffDistance'};
WhiskStimCriteria.Comparison = {'gt'};
WhiskStimCriteria.Value = {5};
% criteria for the resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestStimCriteria.Fieldname = {'stimDistances'};
RestStimCriteria.Comparison = {'gt'};
RestStimCriteria.Value = {5};
% criteria for the stimulation
StimCriteriaA.Value = {'RPadSol'};
StimCriteriaA.Fieldname = {'solenoidName'};
StimCriteriaA.Comparison = {'equal'};
StimCriteriaB.Value = {'LPadSol'};
StimCriteriaB.Fieldname = {'solenoidName'};
StimCriteriaB.Comparison = {'equal'};
% loop variables
hemispheres = {'LH','RH','fLH','fRH'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        % lowpass filter
        samplingRate = RestData.(dataType).(hemisphere).samplingRate;
        [z,p,k] = butter(4,1/(samplingRate/2),'low');
        [sos,g] = zp2sos(z,p,k);
        %% Rest
        [restLogical] = FilterEvents_IOS(RestData.(dataType).(hemisphere),RestCriteria);
        [stimLogical] = FilterEvents_IOS(RestData.(dataType).(hemisphere),RestStimCriteria);
        combRestLogical = logical(restLogical.*stimLogical);
        restFileIDs = RestData.(dataType).(hemisphere).fileIDs(combRestLogical,:);
        restEventTimes = RestData.(dataType).(hemisphere).eventTimes(combRestLogical,:);
        restDurations = RestData.(dataType).(hemisphere).durations(combRestLogical,:);
        restingData = RestData.(dataType).(hemisphere).data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [finalRestData,~,~,~] = RemoveInvalidData_IOS(restingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        % filter
        for gg = 1:length(finalRestData)
            procRestData{gg,1} = filtfilt(sos,g,finalRestData{gg,1});
            restMean(gg,1) = mean(procRestData{gg,1}(1:end));
        end
        % save results
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).Rest.indData = procRestData;
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).Rest.mean = restMean;
        %% Whisk
        [whiskLogical] = FilterEvents_IOS(EventData.(dataType).(hemisphere).whisk,WhiskCriteria);
        [stimLogical] = FilterEvents_IOS(EventData.(dataType).(hemisphere).whisk,WhiskStimCriteria);
        combWhiskLogical = logical(whiskLogical.*stimLogical);
        whiskFileIDs = EventData.(dataType).(hemisphere).whisk.fileIDs(combWhiskLogical,:);
        whiskEventTimes = EventData.(dataType).(hemisphere).whisk.eventTime(combWhiskLogical,:);
        whiskDurations = EventData.(dataType).(hemisphere).whisk.duration(combWhiskLogical,:);
        whiskData = EventData.(dataType).(hemisphere).whisk.data(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [finalWhiskData,~,~,~] = RemoveInvalidData_IOS(whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        % filter and mean-subtract 2 seconds prior to whisk
        for gg = 1:size(finalWhiskData,1)
            procWhiskData_temp = filtfilt(sos,g,finalWhiskData(gg,:));
            procWhiskData(gg,:) = procWhiskData_temp - mean(procWhiskData_temp(1:params.Offset*samplingRate));
            indWhisk{gg,1} = procWhiskData(gg,:);
            whiskCBVMean{gg,1} = mean(procWhiskData(gg,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        end
        % save results
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).Whisk.indData = indWhisk;
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).Whisk.mean = cell2mat(whiskCBVMean);
        %% Stim
        if any(strcmp(hemisphere,{'LH','fLH'})) == true
            StimCriteria = StimCriteriaA;
        elseif any(strcmp(hemisphere,{'RH','fRH'})) == true
            StimCriteria = StimCriteriaB;
        end
        stimFilter = FilterEvents_IOS(EventData.(dataType).(hemisphere).stim,StimCriteria);
        [stimFileIDs] = EventData.(dataType).LH.stim.fileIDs(stimFilter,:);
        [stimEventTimes] = EventData.(dataType).LH.stim.eventTime(stimFilter,:);
        stimDurations = zeros(length(stimEventTimes),1);
        [stimData] = EventData.(dataType).LH.stim.data(stimFilter,:);
        % keep only the data that occurs within the manually-approved awake regions
        [finalStimData,~,~,~] = RemoveInvalidData_IOS(stimData,stimFileIDs,stimDurations,stimEventTimes,ManualDecisions);
        % filter and mean-subtract 2 seconds prior to stimulus
        for gg = 1:size(finalStimData,1)
            procStimData_temp = filtfilt(sos,g,finalStimData(gg,:));
            procStimData(gg,:) = procStimData_temp - mean(procStimData_temp(1:params.Offset*samplingRate));
            indStim{gg,1} = procStimData(gg,:);
            stimMean(gg,1) = mean(procStimData(gg,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
        end
        % save results
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).Stim.indData = indStim;
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).Stim.mean = stimMean;
        %% NREM
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).(hemisphere),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        % filter and take mean during NREM epochs
        for nn = 1:length(nremData)
            indNREM{nn,1} = filtfilt(sos,g,nremData{nn,1});
            nremMean(nn,1) = mean(filtfilt(sos,g,nremData{nn,1}(1:end)));
        end
        % save results
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).NREM.indData = indNREM;
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).NREM.mean = nremMean;
        %% REM
        [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).(hemisphere),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        % filter and take mean during REM epochs
        for nn = 1:length(remData)
            indREM{nn,1} = filtfilt(sos,g,remData{nn,1});
            remMean(nn,1) = mean(filtfilt(sos,g,remData{nn,1}(1:end)));
        end
        % save results
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).REM.indData = indREM;
        Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).REM.mean = remMean;
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_IntSig_GCaMP.mat','Results_IntSig_GCaMP')
cd([rootFolder delim 'Data'])