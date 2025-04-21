function [Results_Evoked_Pulse] = AnalyzeEvokedResponses_Pulse_eLife2025(animalID,group,setName,rootFolder,delim,Results_Evoked_Pulse)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim setName delim animalID delim 'Imaging'];
cd(dataLocation)
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
% find and load ManualDecisions struct
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% find and load RestingBaselines struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load AllSpecStruct struct
allSpecStructFileStruct = dir('*_AllSpecStructB.mat');
allSpecStructFile = {allSpecStructFileStruct.name}';
allSpecStructFileID = char(allSpecStructFile);
load(allSpecStructFileID,'-mat')
% criteria for whisking
WhiskCriteriaA.Fieldname = {'duration','duration','puffDistance'};
WhiskCriteriaA.Comparison = {'gt','lt','gt'};
WhiskCriteriaA.Value = {0.5,2,5};
WhiskCriteriaB.Fieldname = {'duration','duration','puffDistance'};
WhiskCriteriaB.Comparison = {'gt','lt','gt'};
WhiskCriteriaB.Value = {2,5,5};
WhiskCriteriaC.Fieldname = {'duration','puffDistance'};
WhiskCriteriaC.Comparison = {'gt','gt'};
WhiskCriteriaC.Value = {5,5};
WhiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% criteria for stimulation
StimCriteriaA.Value = {'LPadSol'};
StimCriteriaA.Fieldname = {'solenoidName'};
StimCriteriaA.Comparison = {'equal'};
StimCriteriaB.Value = {'RPadSol'};
StimCriteriaB.Fieldname = {'solenoidName'};
StimCriteriaB.Comparison = {'equal'};
StimCriteriaC.Value = {'AudSol'};
StimCriteriaC.Fieldname = {'solenoidName'};
StimCriteriaC.Comparison = {'equal'};
stimCriteriaNames = {'stimCriteriaA','stimCriteriaB','stimCriteriaC'};
%% whisking
for bb = 1:length(WhiskCriteriaNames)
    samplingRate = EventData.CBV_HbT.adjBarrels.whisk.samplingRate;
    trialDuration_sec = EventData.CBV_HbT.adjBarrels.whisk.trialDuration_sec;
    timeVector = (0:(EventData.CBV_HbT.adjBarrels.whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV_HbT.adjBarrels.whisk.epoch.offset;
    offset = EventData.CBV_HbT.adjBarrels.whisk.epoch.offset;
    whiskCriteriaName = WhiskCriteriaNames{1,bb};
    if strcmp(whiskCriteriaName,'ShortWhisks') == true
        WhiskCriteria = WhiskCriteriaA;
    elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
        WhiskCriteria = WhiskCriteriaB;
    elseif strcmp(whiskCriteriaName,'LongWhisks') == true
        WhiskCriteria = WhiskCriteriaC;
    end
    % pull data from EventData.mat structure
    [whiskLogical] = FilterEvents_IOS_eLife2025(EventData.CBV_HbT.adjBarrels.whisk,WhiskCriteria);
    combWhiskLogical = logical(whiskLogical);
    [allWhiskHbTData] = EventData.CBV_HbT.adjBarrels.whisk.data(combWhiskLogical,:);
    [allWhiskFileIDs] = EventData.CBV_HbT.adjBarrels.whisk.fileIDs(combWhiskLogical,:);
    [allWhiskEventTimes] = EventData.CBV_HbT.adjBarrels.whisk.eventTime(combWhiskLogical,:);
    allWhiskDurations = EventData.CBV_HbT.adjBarrels.whisk.duration(combWhiskLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalWhiskHbTData,~,~,finalWhiskFileEventTimes] = RemoveInvalidData_IOS_eLife2025(allWhiskHbTData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
    % lowpass filter each whisking event and mean-subtract by the first 2 seconds
    clear procWhiskHbTData
    dd = 1;
    for cc = 1:size(finalWhiskHbTData,1)
        whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 2;
        whiskEndTime = whiskStartTime + 12;
        if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
            whiskHbTarray = finalWhiskHbTData(cc,:);
            filtWhiskHbTarray = sgolayfilt(whiskHbTarray,3,17);
            procWhiskHbTData(dd,:) = filtWhiskHbTarray - mean(filtWhiskHbTarray(1:(offset*samplingRate)));
            dd = dd + 1;
        end
    end
    meanWhiskHbTData = mean(procWhiskHbTData,1);
    % save results
    Results_Evoked_Pulse.(group).(animalID).Whisk.(whiskCriteriaName).indHbT = procWhiskHbTData;
    Results_Evoked_Pulse.(group).(animalID).Whisk.(whiskCriteriaName).HbT = meanWhiskHbTData;
    Results_Evoked_Pulse.(group).(animalID).Whisk.(whiskCriteriaName).timeVector = timeVector;
end
%% stimulation
for gg = 1:length(stimCriteriaNames)
    stimCriteriaName = stimCriteriaNames{1,gg};
    if strcmp(stimCriteriaName,'stimCriteriaA') == true
        StimCriteria = StimCriteriaA;
        solenoid = 'LPadSol';
    elseif strcmp(stimCriteriaName,'stimCriteriaB') == true
        StimCriteria = StimCriteriaB;
        solenoid = 'RPadSol';
    elseif strcmp(stimCriteriaName,'stimCriteriaC') == true
        StimCriteria = StimCriteriaC;
        solenoid = 'AudSol';
    end
    % pull data from EventData.mat structure
    allStimFilter = FilterEvents_IOS_eLife2025(EventData.CBV_HbT.adjBarrels.stim,StimCriteria);
    [allStimHbTData] = EventData.CBV_HbT.adjBarrels.stim.data(allStimFilter,:);
    [allStimFileIDs] = EventData.CBV_HbT.adjBarrels.stim.fileIDs(allStimFilter,:);
    [allStimEventTimes] = EventData.CBV_HbT.adjBarrels.stim.eventTime(allStimFilter,:);
    allStimDurations = zeros(length(allStimEventTimes),1);
    % keep only the data that occurs within the manually-approved awake regions
    [finalStimHbTData,~,~,finalStimFileEventTimes] = RemoveInvalidData_IOS_eLife2025(allStimHbTData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
    % lowpass filter each stim event and mean-subtract by the first 2 seconds
    clear procStimHbTData
    ii = 1;
    for hh = 1:size(finalStimHbTData,1)
        stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 2;
        stimEndTime = stimStartTime + 12;
        if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
            stimHbTarray = finalStimHbTData(hh,:);
            filtStimHbTarray = sgolayfilt(stimHbTarray,3,17);
            procStimHbTData(ii,:) = filtStimHbTarray - mean(filtStimHbTarray(1:(offset*samplingRate)));
            ii = ii + 1;
        end
    end
    meanStimHbTData = mean(procStimHbTData,1);
    % save results
    Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).indHbT = procStimHbTData;
    Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).HbT = meanStimHbTData;
    Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).timeVector = timeVector;
    Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).count = size(procStimHbTData,1);
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_Evoked_Pulse.mat','Results_Evoked_Pulse')
cd([rootFolder delim 'Data'])