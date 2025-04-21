function [Results_Evoked_Ephys] = AnalyzeEvokedResponses_Ephys_eLife2025(animalID,group,setName,rootFolder,delim,Results_Evoked_Ephys)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim setName delim animalID delim 'Imaging'];
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
hemispheres = {'LH','RH'};
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    neuralDataType = ['cortical_' hemisphere];
    %% whisking
    samplingRate = EventData.HbT.(hemisphere).whisk.samplingRate;
    specSamplingRate = 10;
    trialDuration_sec = EventData.HbT.(hemisphere).whisk.trialDuration_sec;
    timeVector = (0:(EventData.HbT.(hemisphere).whisk.epoch.duration*samplingRate))/samplingRate - EventData.HbT.(hemisphere).whisk.epoch.offset;
    offset = EventData.HbT.(hemisphere).whisk.epoch.offset;
    for bb = 1:length(WhiskCriteriaNames)
        whiskCriteriaName = WhiskCriteriaNames{1,bb};
        if strcmp(whiskCriteriaName,'ShortWhisks') == true
            WhiskCriteria = WhiskCriteriaA;
        elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
            WhiskCriteria = WhiskCriteriaB;
        elseif strcmp(whiskCriteriaName,'LongWhisks') == true
            WhiskCriteria = WhiskCriteriaC;
        end
        % pull data from EventData.mat structure
        [whiskLogical] = FilterEvents_IOS_eLife2025(EventData.HbT.(hemisphere).whisk,WhiskCriteria);
        combWhiskLogical = logical(whiskLogical);
        [whiskHbTData] = EventData.HbT.(hemisphere).whisk.data(combWhiskLogical,:);
        [whiskCortMUAData] = EventData.(neuralDataType).muaPower.whisk.NormData(combWhiskLogical,:);
        [whiskHipMUAData] = EventData.hippocampus.muaPower.whisk.NormData(combWhiskLogical,:);
        [whiskCortGamData] = EventData.(neuralDataType).gammaBandPower.whisk.NormData(combWhiskLogical,:);
        [whiskHipGamData] = EventData.hippocampus.gammaBandPower.whisk.NormData(combWhiskLogical,:);
        [whiskFileIDs] = EventData.HbT.(hemisphere).whisk.fileIDs(combWhiskLogical,:);
        [whiskEventTimes] = EventData.HbT.(hemisphere).whisk.eventTime(combWhiskLogical,:);
        whiskDurations = EventData.HbT.(hemisphere).whisk.duration(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [finalWhiskHbTData,finalWhiskFileIDs,~,finalWhiskFileEventTimes] = RemoveInvalidData_IOS_eLife2025(whiskHbTData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [finalWhiskCortMUAData,~,~,~] = RemoveInvalidData_IOS_eLife2025(whiskCortMUAData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [finalWhiskHipMUAData,~,~,~] = RemoveInvalidData_IOS_eLife2025(whiskHipMUAData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [finalWhiskCortGamData,~,~,~] = RemoveInvalidData_IOS_eLife2025(whiskCortGamData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [finalWhiskHipGamData,~,~,~] = RemoveInvalidData_IOS_eLife2025(whiskHipGamData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procWhiskHbTData procWhiskCortMUAData procWhiskHipMUAData procWhiskCortGamData procWhiskHipGamData finalWhiskStartTimes finalWhiskEndTimes finalWhiskFiles
        dd = 1;
        for cc = 1:size(finalWhiskHbTData,1)
            whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 2;
            whiskEndTime = whiskStartTime + 12;
            finalWhiskFileID = finalWhiskFileIDs{cc,1};
            if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                whiskHbTarray = finalWhiskHbTData(cc,:);
                whiskCortMUAarray = finalWhiskCortMUAData(cc,:);
                whiskHipMUAarray = finalWhiskHipMUAData(cc,:);
                whiskCortGamArray = finalWhiskCortGamData(cc,:);
                whiskHipGamArray = finalWhiskHipGamData(cc,:);
                filtWhiskHbTarray = sgolayfilt(whiskHbTarray,3,17);
                filtWhiskCortMUAarray = sgolayfilt(whiskCortMUAarray,3,17);
                filtWhiskHipMUAarray = sgolayfilt(whiskHipMUAarray,3,17);
                filtWhiskCortGamArray = sgolayfilt(whiskCortGamArray,3,17);
                filtWhiskHipGamArray = sgolayfilt(whiskHipGamArray,3,17);
                procWhiskHbTData(dd,:) = filtWhiskHbTarray - mean(filtWhiskHbTarray(1:(offset*samplingRate)));
                procWhiskCortMUAData(dd,:) = filtWhiskCortMUAarray - mean(filtWhiskCortMUAarray(1:(offset*samplingRate)));
                procWhiskHipMUAData(dd,:) = filtWhiskHipMUAarray - mean(filtWhiskHipMUAarray(1:(offset*samplingRate)));
                procWhiskCortGamData(dd,:) = filtWhiskCortGamArray - mean(filtWhiskCortGamArray(1:(offset*samplingRate)));
                procWhiskHipGamData(dd,:) = filtWhiskHipGamArray - mean(filtWhiskHipGamArray(1:(offset*samplingRate)));
                finalWhiskStartTimes(dd,1) = whiskStartTime;
                finalWhiskEndTimes(dd,1) = whiskEndTime;
                finalWhiskFiles{dd,1} = finalWhiskFileID;
                dd = dd + 1;
            end
        end
        meanWhiskHbTData = mean(procWhiskHbTData,1);
        meanWhiskCortMUAData = mean(procWhiskCortMUAData,1)*100;
        meanWhiskHipMUAData = mean(procWhiskHipMUAData,1)*100;
        meanWhiskCortGamData = mean(procWhiskCortGamData,1)*100;
        meanWhiskHipGamData = mean(procWhiskHipGamData,1)*100;
        % extract LFP from spectrograms associated with the whisking indecies
        whiskCortZhold = [];
        whiskHipZhold = [];
        for ee = 1:length(finalWhiskFiles)
            % load normalized one-second bin data from each file
            whiskFileID = finalWhiskFiles{ee,1};
            whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
            whiskSpecField = neuralDataType;
            for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                    whiskCortS_Data = AllSpecData.(whiskSpecField).normS{ff,1};
                    whiskHipS_Data = AllSpecData.hippocampus.normS{ff,1};
                    F = AllSpecData.(whiskSpecField).F{ff,1};
                    T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                end
            end
            whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
            whiskStartTimeIndex = whiskStartTimeIndex(1);
            whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
            whiskDurationIndex = whiskDurationIndex(end);
            whiskCortS_Vals = whiskCortS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
            whiskHipS_Vals = whiskHipS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
            % mean subtract each row with detrend - transpose since detrend goes down columns
            transpWhiskCortS_Vals = whiskCortS_Vals';
            transpWhiskHipS_Vals = whiskHipS_Vals';
            dTWhiskCortS_Vals = transpWhiskCortS_Vals;
            dTWhiskCortS_Vals = dTWhiskCortS_Vals(1:12*specSamplingRate + 1,:);
            dTWhiskHipS_Vals = transpWhiskHipS_Vals;
            dTWhiskHipS_Vals = dTWhiskHipS_Vals(1:12*specSamplingRate + 1,:);
            % transpose back to original orientation
            whiskCortZhold = cat(3,whiskCortZhold,dTWhiskCortS_Vals');
            whiskHipZhold = cat(3,whiskHipZhold,dTWhiskHipS_Vals');
        end
        % figure time/frequency axis and average each S data matrix through time
        meanWhiskCortS = mean(whiskCortZhold,3);
        meanWhiskHipS = mean(whiskHipZhold,3);
        T2 = -2:(1/specSamplingRate):10;
        % save results
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).indHbT = procWhiskHbTData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).HbT = meanWhiskHbTData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).cortMUA = meanWhiskCortMUAData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).hipMUA = meanWhiskHipMUAData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).cortGam = meanWhiskCortGamData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).hipGam = meanWhiskHipGamData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).cortLFP = meanWhiskCortS;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).hipLFP = meanWhiskHipS;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).T = T2;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).F = F;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Whisk.(whiskCriteriaName).timeVector = timeVector;
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
        stimFilter = FilterEvents_IOS_eLife2025(EventData.HbT.(hemisphere).stim,StimCriteria);
        [stimHbTData] = EventData.HbT.(hemisphere).stim.data(stimFilter,:);
        [stimCortMUAData] = EventData.(neuralDataType).muaPower.stim.NormData(stimFilter,:);
        [stimHipMUAData] = EventData.hippocampus.muaPower.stim.NormData(stimFilter,:);
        [stimCortGamData] = EventData.(neuralDataType).gammaBandPower.stim.NormData(stimFilter,:);
        [stimHipGamData] = EventData.hippocampus.gammaBandPower.stim.NormData(stimFilter,:);
        [stimFileIDs] = EventData.HbT.(hemisphere).stim.fileIDs(stimFilter,:);
        [stimEventTimes] = EventData.HbT.(hemisphere).stim.eventTime(stimFilter,:);
        stimDurations = zeros(length(stimEventTimes),1);
        % keep only the data that occurs within the manually-approved awake regions
        [finalStimHbTData,finalStimFileIDs,~,finalStimFileEventTimes] = RemoveInvalidData_IOS_eLife2025(stimHbTData,stimFileIDs,stimDurations,stimEventTimes,ManualDecisions);
        [finalStimCortMUAData,~,~,~] = RemoveInvalidData_IOS_eLife2025(stimCortMUAData,stimFileIDs,stimDurations,stimEventTimes,ManualDecisions);
        [finalStimHipMUAData,~,~,~] = RemoveInvalidData_IOS_eLife2025(stimHipMUAData,stimFileIDs,stimDurations,stimEventTimes,ManualDecisions);
        [finalStimCortGamData,~,~,~] = RemoveInvalidData_IOS_eLife2025(stimCortGamData,stimFileIDs,stimDurations,stimEventTimes,ManualDecisions);
        [finalStimHipGamData,~,~,~] = RemoveInvalidData_IOS_eLife2025(stimHipGamData,stimFileIDs,stimDurations,stimEventTimes,ManualDecisions);
        % lowpass filter each stim event and mean-subtract by the first 2 seconds
        clear procStimHbTData procStimCortMUAData procStimHipMUAData procStimCortGamData procStimHipGamData finalStimStartTimes finalStimEndTimes finalStimFiles
        ii = 1;
        for hh = 1:size(finalStimHbTData,1)
            stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 2;
            stimEndTime = stimStartTime + 12;
            finalStimFileID = finalStimFileIDs{hh,1};
            if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                stimHbTarray = finalStimHbTData(hh,:);
                stimCortMUAarray = finalStimCortMUAData(hh,:);
                stimHipMUAarray = finalStimHipMUAData(hh,:);
                stimCortGamArray = finalStimCortGamData(hh,:);
                stimHipGamArray = finalStimHipGamData(hh,:);
                filtStimHbTarray = sgolayfilt(stimHbTarray,3,17);
                % filtStimCortMUAarray = sgolayfilt(stimCortMUAarray,3,17);
                filtStimCortMUAarray = stimCortMUAarray;
                filtStimHipMUAarray = sgolayfilt(stimHipMUAarray,3,17);
                % filtStimCortGamArray = sgolayfilt(stimCortGamArray,3,17);
                filtStimCortGamArray = stimCortGamArray;
                filtStimHipGamArray = sgolayfilt(stimHipGamArray,3,17);
                procStimHbTData(hh,:) = filtStimHbTarray - mean(filtStimHbTarray(1:(offset*samplingRate)));
                procStimCortMUAData(hh,:) = filtStimCortMUAarray - mean(filtStimCortMUAarray(1:(offset*samplingRate)));
                procStimHipMUAData(hh,:) = filtStimHipMUAarray - mean(filtStimHipMUAarray(1:(offset*samplingRate)));
                procStimCortGamData(hh,:) = filtStimCortGamArray - mean(filtStimCortGamArray(1:(offset*samplingRate)));
                procStimHipGamData(hh,:) = filtStimHipGamArray - mean(filtStimHipGamArray(1:(offset*samplingRate)));
                finalStimStartTimes(ii,1) = stimStartTime;
                finalStimEndTimes(ii,1) = stimEndTime;
                finalStimFiles{ii,1} = finalStimFileID;
                ii = ii + 1;
            end
        end
        meanStimHbTData = mean(procStimHbTData,1);
        meanStimCortMUAData = mean(procStimCortMUAData,1)*100;
        meanStimHipMUAData = mean(procStimHipMUAData,1)*100;
        meanStimCortGamData = mean(procStimCortGamData,1)*100;
        meanStimHipGamData = mean(procStimHipGamData,1)*100;
        % extract LFP from spectrograms associated with the stimuli indecies
        stimCortZhold = [];
        stimHipZhold = [];
        for jj = 1:length(finalStimFiles)
            % load normalized one-second bin data from each file
            stimFileID = finalStimFiles{jj,1};
            stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
            stimSpecField = neuralDataType;
            for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                    stimCortS_Data = AllSpecData.(stimSpecField).normS{kk,1};
                    stimHipS_Data = AllSpecData.hippocampus.normS{kk,1};
                end
            end
            stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
            stimStartTimeIndex = stimStartTimeIndex(1);
            stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
            stimDurationIndex = stimDurationIndex(end);
            stimCortS_Vals = stimCortS_Data(:,stimStartTimeIndex:stimDurationIndex);
            stimHipS_Vals = stimHipS_Data(:,stimStartTimeIndex:stimDurationIndex);
            % mean subtract each row with detrend
            transpStimCortS_Vals = stimCortS_Vals';
            transpStimHipS_Vals = stimHipS_Vals';
            dTStimCortS_Vals = transpStimCortS_Vals;
            dTStimCortS_Vals = dTStimCortS_Vals(1:12*specSamplingRate + 1,:);
            dTStimHipS_Vals = transpStimHipS_Vals;
            dTStimHipS_Vals = dTStimHipS_Vals(1:12*specSamplingRate + 1,:);
            % transpose back to original orientation
            stimCortZhold = cat(3,stimCortZhold,dTStimCortS_Vals');
            stimHipZhold = cat(3,stimHipZhold,dTStimHipS_Vals');
        end
        % figure time/frequency axis and average each S data matrix through time
        meanStimCortS = mean(stimCortZhold,3);
        meanStimHipS = mean(stimHipZhold,3);
        % save results
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).indHbT = procStimHbTData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).HbT = meanStimHbTData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortMUA = meanStimCortMUAData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).hipMUA = meanStimHipMUAData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortGam = meanStimCortGamData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).hipGam = meanStimHipGamData;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortLFP = meanStimCortS;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).hipLFP = meanStimHipS;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).T = T2;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).F = F;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).timeVector = timeVector;
        Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).count = size(procStimHbTData,1);
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_Evoked_Ephys.mat','Results_Evoked_Ephys')
cd([rootFolder delim 'Data'])