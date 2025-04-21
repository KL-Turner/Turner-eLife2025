function [Results_Evoked_GCaMP] = AnalyzeEvokedResponses_GCaMP_eLife2025(animalID,group,set,rootFolder,delim,Results_Evoked_GCaMP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% find and load EventData struct
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
% loop variables
hemispheres = {'LH','RH','fLH','fRH'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        %% whisking
        samplingRate = EventData.(dataType).(hemisphere).whisk.samplingRate;
        trialDuration_sec = EventData.(dataType).(hemisphere).whisk.trialDuration_sec;
        timeVector = (0:(EventData.(dataType).(hemisphere).whisk.epoch.duration*samplingRate))/samplingRate - EventData.(dataType).(hemisphere).whisk.epoch.offset;
        offset = EventData.(dataType).(hemisphere).whisk.epoch.offset;
        for cc = 1:length(WhiskCriteriaNames)
            whiskCriteriaName = WhiskCriteriaNames{1,cc};
            if strcmp(whiskCriteriaName,'ShortWhisks') == true
                WhiskCriteria = WhiskCriteriaA;
            elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
                WhiskCriteria = WhiskCriteriaB;
            elseif strcmp(whiskCriteriaName,'LongWhisks') == true
                WhiskCriteria = WhiskCriteriaC;
            end
            % pull data from EventData.mat structure
            [whiskLogical] = FilterEvents_IOS_eLife2025(EventData.(dataType).(hemisphere).whisk,WhiskCriteria);
            combWhiskLogical = logical(whiskLogical);
            [whiskData] = EventData.(dataType).(hemisphere).whisk.data(combWhiskLogical,:);
            [whiskFileIDs] = EventData.(dataType).(hemisphere).whisk.fileIDs(combWhiskLogical,:);
            [whiskEventTimes] = EventData.(dataType).(hemisphere).whisk.eventTime(combWhiskLogical,:);
            whiskDurations = EventData.(dataType).(hemisphere).whisk.duration(combWhiskLogical,:);
            % keep only the data that occurs within the manually-approved awake regions
            [finalWhiskData,~,~,finalWhiskFileEventTimes] = RemoveInvalidData_IOS_eLife2025(whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
            % lowpass filter each whisking event and mean-subtract by the first 2 seconds
            clear procWhiskData
            dd = 1;
            for ee = 1:size(finalWhiskData,1)
                whiskStartTime = round(finalWhiskFileEventTimes(ee,1),1) - 2;
                whiskEndTime = whiskStartTime + 12;
                if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                    whiskArray = finalWhiskData(ee,:);
                    filtWhiskArray = sgolayfilt(whiskArray,3,17);
                    procWhiskData(dd,:) = filtWhiskArray - mean(filtWhiskArray(1:(offset*samplingRate)));
                    dd = dd + 1;
                end
            end
            meanWhiskData = mean(procWhiskData,1);
            % save results
            Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).(dataType).indData = procWhiskData;
            Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).(dataType).mean = meanWhiskData;
            Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).(dataType).timeVector = timeVector;
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
            stimFilter = FilterEvents_IOS_eLife2025(EventData.(dataType).(hemisphere).stim,StimCriteria);
            [stimData] = EventData.(dataType).(hemisphere).stim.data(stimFilter,:);
            [stimFileIDs] = EventData.(dataType).(hemisphere).stim.fileIDs(stimFilter,:);
            [stimEventTimes] = EventData.(dataType).(hemisphere).stim.eventTime(stimFilter,:);
            stimDurations = zeros(length(stimEventTimes),1);
            % keep only the data that occurs within the manually-approved awake regions
            [finalStimHbTData,~,~,finalStimFileEventTimes] = RemoveInvalidData_IOS_eLife2025(stimData,stimFileIDs,stimDurations,stimEventTimes,ManualDecisions);
            % lowpass filter each stim event and mean-subtract by the first 2 seconds
            clear procStimData
            ii = 1;
            for hh = 1:size(finalStimHbTData,1)
                stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 2;
                stimEndTime = stimStartTime + 12;
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                    stimArray = finalStimHbTData(hh,:);
                    filtStimArray = sgolayfilt(stimArray,3,17);
                    procStimData(ii,:) = filtStimArray - mean(filtStimArray(1:(offset*samplingRate)));
                    ii = ii + 1;
                end
            end
            meanStimData = mean(procStimData,1);
            % save results
            Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).(dataType).indData = procStimData;
            Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).(dataType).mean = meanStimData;
            Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).(dataType).timeVector = timeVector;
            Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).(dataType).count = size(procStimData,1);
        end
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_Evoked_GCaMP.mat','Results_Evoked_GCaMP')
cd([rootFolder delim 'Data'])