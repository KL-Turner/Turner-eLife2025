function [AnalysisResults] = AnalyzeEvokedResponses_HRF2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
dataTypes = {'adjLH','adjRH'};
%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % find and load EventData.mat struct
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
    % find and load manual baseline event information
    manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
    manualBaselineFile = {manualBaselineFileStruct.name}';
    manualBaselineFileID = char(manualBaselineFile);
    load(manualBaselineFileID)
    % find and load RestingBaselines.mat struct
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % find and load AllSpecStruct.mat struct
    allSpecStructFileStruct = dir('*_AllSpecStructB.mat');
    allSpecStructFile = {allSpecStructFileStruct.name}';
    allSpecStructFileID = char(allSpecStructFile);
    load(allSpecStructFileID)
    % forest ID sctruct
    forestScoringResultsID = 'Forest_ScoringResults.mat';
    load(forestScoringResultsID,'-mat')
    % criteria for whisking
    WhiskCriteriaA.Fieldname = {'duration','duration','puffDistance'};
    WhiskCriteriaA.Comparison = {'gt','lt','gt'};
    WhiskCriteriaA.Value = {0.5,2,5};
    %     WhiskCriteriaB.Fieldname = {'duration','duration','puffDistance'};
    %     WhiskCriteriaB.Comparison = {'gt','lt','gt'};
    %     WhiskCriteriaB.Value = {2,5,5};
    WhiskCriteriaB.Fieldname = {'duration','duration','restTime','puffDistance'};
    WhiskCriteriaB.Comparison = {'lt','gt','gt','gt'};
    WhiskCriteriaB.Value = {3,0.5,4,5};
    WhiskCriteriaC.Fieldname = {'duration','puffDistance'};
    WhiskCriteriaC.Comparison = {'gt','gt'};
    WhiskCriteriaC.Value = {5,5};
    WhiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
    % criteria for stimulation
    %     StimCriteriaA.Value = {'LPadSol'};
    %     StimCriteriaA.Fieldname = {'solenoidName'};
    %     StimCriteriaA.Comparison = {'equal'};
    %     StimCriteriaB.Value = {'RPadSol'};
    %     StimCriteriaB.Fieldname = {'solenoidName'};
    %     StimCriteriaB.Comparison = {'equal'};
    %     StimCriteriaC.Value = {'AudSol'};
    %     StimCriteriaC.Fieldname = {'solenoidName'};
    %     StimCriteriaC.Comparison = {'equal'};
    %     stimCriteriaNames = {'stimCriteriaA','stimCriteriaB','stimCriteriaC'};
    StimCriteriaA.Value = {'LPadSol',0.01,0.01,5};
    StimCriteriaA.Fieldname = {'solenoidName','whiskScore_Pre','movementScore_Pre','PuffDistance'};
    StimCriteriaA.Comparison  = {'equal','lt','lt','gt'};
    StimCriteriaB.Value = {'RPadSol',0.01,0.01,5};
    StimCriteriaB.Fieldname = {'solenoidName','whiskScore_Pre','movementScore_Pre','PuffDistance'};
    StimCriteriaB.Comparison  = {'equal','lt','lt','gt'};
    StimCriteriaC.Value = {'AudSol',0.01,0.01,5};
    StimCriteriaC.Fieldname = {'solenoidName','whiskScore_Pre','movementScore_Pre','PuffDistance'};
    StimCriteriaC.Comparison  = {'equal','lt','lt','gt'};
    stimCriteriaNames = {'stimCriteriaA','stimCriteriaB','stimCriteriaC'};
    %% analyze whisking-evoked responses
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        neuralDataType = ['cortical_' dataType(4:end)];
        % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
        samplingRate = EventData.CBV_HbT.(dataType).whisk.samplingRate;
        specSamplingRate = 10;
        trialDuration_sec = EventData.CBV_HbT.(dataType).whisk.trialDuration_sec;
        timeVector = (0:(EventData.CBV_HbT.(dataType).whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV_HbT.(dataType).whisk.epoch.offset;
        offset = EventData.CBV_HbT.(dataType).whisk.epoch.offset;
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
            [whiskLogical] = FilterEvents_HRF2020(EventData.CBV_HbT.(dataType).whisk,WhiskCriteria);
            combWhiskLogical = logical(whiskLogical);
            [allWhiskHbTData] = EventData.CBV_HbT.(dataType).whisk.data(combWhiskLogical,:);
            [allWhiskCBVData] = EventData.CBV.(dataType).whisk.NormData(combWhiskLogical,:);
            [allWhiskCorticalMUAData] = EventData.(neuralDataType).muaPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskHippocampalMUAData] = EventData.hippocampus.muaPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskCorticalGamData] = EventData.(neuralDataType).gammaBandPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskHippocampalGamData] = EventData.hippocampus.gammaBandPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskFileIDs] = EventData.CBV_HbT.(dataType).whisk.fileIDs(combWhiskLogical,:);
            [allWhiskEventTimes] = EventData.CBV_HbT.(dataType).whisk.eventTime(combWhiskLogical,:);
            allWhiskDurations = EventData.CBV_HbT.(dataType).whisk.duration(combWhiskLogical,:);
            % keep only the data that occurs within the manually-approved awake regions
            [finalWhiskHbTData,finalWhiskFileIDs,~,finalWhiskFileEventTimes] = RemoveInvalidData_HRF2020(allWhiskHbTData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCBVData,~,~,~] = RemoveInvalidData_HRF2020(allWhiskCBVData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCorticalMUAData,~,~,~] = RemoveInvalidData_HRF2020(allWhiskCorticalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskHippocampalMUAData,~,~,~] = RemoveInvalidData_HRF2020(allWhiskHippocampalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCorticalGamData,~,~,~] = RemoveInvalidData_HRF2020(allWhiskCorticalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskHippocampalGamData,~,~,~] = RemoveInvalidData_HRF2020(allWhiskHippocampalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            % lowpass filter each whisking event and mean-subtract by the first 2 seconds
            clear procWhiskHbTData procWhiskCBVData procWhiskCorticalMUAData procWhiskHippocampalMUAData procWhiskCorticalGamData procWhiskHippocampalGamData finalWhiskStartTimes finalWhiskEndTimes finalWhiskFiles
            dd = 1;
            for cc = 1:size(finalWhiskHbTData,1)
                whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 2;
                whiskEndTime = whiskStartTime + 12;
                finalWhiskFileID = finalWhiskFileIDs{cc,1};
                if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                    whiskHbTarray = finalWhiskHbTData(cc,:);
                    whiskCBVarray = finalWhiskCBVData(cc,:);
                    whiskCorticalMUAarray = finalWhiskCorticalMUAData(cc,:);
                    whiskHippocampalMUAarray = finalWhiskHippocampalMUAData(cc,:);
                    whiskCorticalGamArray = finalWhiskCorticalGamData(cc,:);
                    whiskHippocampalGamArray = finalWhiskHippocampalGamData(cc,:);
                    filtWhiskHbTarray = sgolayfilt(whiskHbTarray,3,17);
                    filtWhiskCBVarray = sgolayfilt(whiskCBVarray,3,17);
%                     filtWhiskCorticalMUAarray = sgolayfilt(whiskCorticalMUAarray,3,17);
%                     filtWhiskHippocampalMUAarray = sgolayfilt(whiskHippocampalMUAarray,3,17);
%                     filtWhiskCorticalGamArray = sgolayfilt(whiskCorticalGamArray,3,17);
%                     filtWhiskHippocampalGamArray = sgolayfilt(whiskHippocampalGamArray,3,17);
                    filtWhiskCorticalMUAarray = (whiskCorticalMUAarray);
                    filtWhiskHippocampalMUAarray = (whiskHippocampalMUAarray);
                    filtWhiskCorticalGamArray = (whiskCorticalGamArray);
                    filtWhiskHippocampalGamArray = (whiskHippocampalGamArray);
                    procWhiskHbTData(dd,:) = filtWhiskHbTarray - mean(filtWhiskHbTarray(1:(offset*samplingRate))); %#ok<*AGROW>
                    procWhiskCBVData(dd,:) = filtWhiskCBVarray - mean(filtWhiskCBVarray(1:(offset*samplingRate)));
                    procWhiskCorticalMUAData(dd,:) = filtWhiskCorticalMUAarray - mean(filtWhiskCorticalMUAarray(1:(offset*samplingRate)));
                    procWhiskHippocampalMUAData(dd,:) = filtWhiskHippocampalMUAarray - mean(filtWhiskHippocampalMUAarray(1:(offset*samplingRate)));
                    procWhiskCorticalGamData(dd,:) = filtWhiskCorticalGamArray - mean(filtWhiskCorticalGamArray(1:(offset*samplingRate)));
                    procWhiskHippocampalGamData(dd,:) = filtWhiskHippocampalGamArray - mean(filtWhiskHippocampalGamArray(1:(offset*samplingRate)));
                    finalWhiskStartTimes(dd,1) = whiskStartTime;
                    finalWhiskEndTimes(dd,1) = whiskEndTime;
                    finalWhiskFiles{dd,1} = finalWhiskFileID;
                    dd = dd + 1;
                end
            end
            meanWhiskHbTData = mean(procWhiskHbTData,1);
            stdWhiskHbTData = std(procWhiskHbTData,0,1);
            meanWhiskCBVData = mean(procWhiskCBVData,1)*100;
            stdWhiskCBVData = std(procWhiskCBVData,0,1)*100;
            meanWhiskCorticalMUAData = mean(procWhiskCorticalMUAData,1)*100;
            stdWhiskCorticalMUAData = std(procWhiskCorticalMUAData,0,1)*100;
            meanWhiskHippocampalMUAData = mean(procWhiskHippocampalMUAData,1)*100;
            stdWhiskHippocampalMUAData = std(procWhiskHippocampalMUAData,0,1)*100;
            meanWhiskCorticalGamData = mean(procWhiskCorticalGamData,1)*100;
            stdWhiskCorticalGamData = std(procWhiskCorticalGamData,0,1)*100;
            meanWhiskHippocampalGamData = mean(procWhiskHippocampalGamData,1)*100;
            stdWhiskHippocampalGamData = std(procWhiskHippocampalGamData,0,1)*100;
            % extract LFP from spectrograms associated with the whisking indecies
            whiskCorticalZhold = [];
            whiskHippocampalZhold = [];
            for ee = 1:length(finalWhiskFiles)
                % load normalized one-second bin data from each file
                whiskFileID = finalWhiskFiles{ee,1};
                whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
                whiskSpecField = neuralDataType;
                for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                    if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                        whiskCorticalS_Data = AllSpecData.(whiskSpecField).normS{ff,1};
                        whiskHippocampalS_Data = AllSpecData.hippocampus.normS{ff,1};
                        F = AllSpecData.(whiskSpecField).F{ff,1};
                        T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                    end
                end
                whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
                whiskStartTimeIndex = whiskStartTimeIndex(1);
                whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
                whiskDurationIndex = whiskDurationIndex(end);
                whiskCorticalS_Vals = whiskCorticalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
                whiskHippocampalS_Vals = whiskHippocampalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
                % mean subtract each row with detrend - transpose since detrend goes down columns
                transpWhiskCorticalS_Vals = whiskCorticalS_Vals';
                transpWhiskHippocampalS_Vals = whiskHippocampalS_Vals';
                dTWhiskCorticalS_Vals = transpWhiskCorticalS_Vals;
                dTWhiskCorticalS_Vals = dTWhiskCorticalS_Vals(1:12*specSamplingRate + 1,:);
                dTWhiskHippocampalS_Vals = transpWhiskHippocampalS_Vals;
                dTWhiskHippocampalS_Vals = dTWhiskHippocampalS_Vals(1:12*specSamplingRate + 1,:);
                % transpose back to original orientation
                whiskCorticalZhold = cat(3,whiskCorticalZhold,dTWhiskCorticalS_Vals');
                whiskHippocampalZhold = cat(3,whiskHippocampalZhold,dTWhiskHippocampalS_Vals');
            end
            % figure time/frequency axis and average each S data matrix through time
            meanWhiskCorticalS = mean(whiskCorticalZhold,3);
            meanWhiskHippocampalS = mean(whiskHippocampalZhold,3);
            T2 = -2:(1/specSamplingRate):10;
            % save results
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).CBV_HbT.HbT = meanWhiskHbTData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).CBV_HbT.HbTStD = stdWhiskHbTData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).CBV.CBV = meanWhiskCBVData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).CBV.CBVStD = stdWhiskCBVData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalData = meanWhiskHippocampalMUAData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalStD = stdWhiskHippocampalMUAData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).Gam.corticalData = meanWhiskCorticalGamData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).Gam.corticalStD = stdWhiskCorticalGamData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).Gam.hippocampalData = meanWhiskHippocampalGamData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).Gam.hippocampalStD = stdWhiskHippocampalGamData;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).timeVector = timeVector;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).LFP.corticalS = meanWhiskCorticalS;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).LFP.hippocampalS = meanWhiskHippocampalS;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).LFP.T = T2;
            AnalysisResults.Evoked.(animalID).Whisk.(dataType).(whiskCriteriaName).LFP.F = F;
        end
        %% analyze stimulus-evoked responses
        for gg = 1:length(stimCriteriaNames)
            stimCriteriaName = stimCriteriaNames{1,gg};
            if strcmp(stimCriteriaName,'stimCriteriaA') == true
                stimCriteria = StimCriteriaA;
                solenoid = 'LPadSol';
            elseif strcmp(stimCriteriaName,'stimCriteriaB') == true
                stimCriteria = StimCriteriaB;
                solenoid = 'RPadSol';
            elseif strcmp(stimCriteriaName,'stimCriteriaC') == true
                stimCriteria = StimCriteriaC;
                solenoid = 'AudSol';
            end
            % pull data from EventData.mat structure
            allStimFilter = FilterEvents_HRF2020(EventData.CBV_HbT.(dataType).stim,stimCriteria);
            [allStimHbTData] = EventData.CBV_HbT.(dataType).stim.data(allStimFilter,:);
            [allStimCBVData] = EventData.CBV.(dataType).stim.NormData(allStimFilter,:);
            [allStimCortMUAData] = EventData.(neuralDataType).muaPower.stim.NormData(allStimFilter,:);
            [allStimHipMUAData] = EventData.hippocampus.muaPower.stim.NormData(allStimFilter,:);
            [allStimCortGamData] = EventData.(neuralDataType).gammaBandPower.stim.NormData(allStimFilter,:);
            [allStimHipGamData] = EventData.hippocampus.gammaBandPower.stim.NormData(allStimFilter,:);
            [allStimFileIDs] = EventData.CBV_HbT.(dataType).stim.fileIDs(allStimFilter,:);
            [allStimEventTimes] = EventData.CBV_HbT.(dataType).stim.eventTime(allStimFilter,:);
            allStimDurations = zeros(length(allStimEventTimes),1);
            % keep only the data that occurs within the manually-approved awake regions
            [finalStimHbTData,finalStimFileIDs,~,finalStimFileEventTimes] = RemoveInvalidData_HRF2020(allStimHbTData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCBVData,~,~,~] = RemoveInvalidData_HRF2020(allStimCBVData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortMUAData,~,~,~] = RemoveInvalidData_HRF2020(allStimCortMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimHipMUAData,~,~,~] = RemoveInvalidData_HRF2020(allStimHipMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortGamData,~,~,~] = RemoveInvalidData_HRF2020(allStimCortGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimHipGamData,~,~,~] = RemoveInvalidData_HRF2020(allStimHipGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            % lowpass filter each stim event and mean-subtract by the first 2 seconds
            clear procStimHbTData procStimCBVData procStimCortMUAData procStimHipMUAData procStimCortGamData procStimHipGamData finalStimStartTimes finalStimEndTimes finalStimFiles
            ii = 1;
            for hh = 1:size(finalStimHbTData,1)
                stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 2;
                stimEndTime = stimStartTime + 12;
                finalStimFileID = finalStimFileIDs{hh,1};
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                    stimHbTarray = finalStimHbTData(hh,:);
                    stimCBVarray = finalStimCBVData(hh,:);
                    stimCortMUAarray = finalStimCortMUAData(hh,:);
                    stimHipMUAarray = finalStimHipMUAData(hh,:);
                    stimCortGamArray = finalStimCortGamData(hh,:);
                    stimHipGamArray = finalStimHipGamData(hh,:);
                    filtStimHbTarray = sgolayfilt(stimHbTarray,3,17);
                    filtStimCBVarray = sgolayfilt(stimCBVarray,3,17);
%                     filtStimCortMUAarray = sgolayfilt(stimCortMUAarray,3,17);
%                     filtStimHipMUAarray = sgolayfilt(stimHipMUAarray,3,17);
%                     filtStimCortGamArray = sgolayfilt(stimCortGamArray,3,17);
%                     filtStimHipGamArray = sgolayfilt(stimHipGamArray,3,17);
                    filtStimCortMUAarray = (stimCortMUAarray);
                    filtStimHipMUAarray = (stimHipMUAarray);
                    filtStimCortGamArray = (stimCortGamArray);
                    filtStimHipGamArray = (stimHipGamArray);
                    procStimHbTData(hh,:) = filtStimHbTarray - mean(filtStimHbTarray(1:(offset*samplingRate)));
                    procStimCBVData(hh,:) = filtStimCBVarray - mean(filtStimCBVarray(1:(offset*samplingRate)));
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
            stdStimHbTData = std(procStimHbTData,0,1);
            meanStimCBVData = mean(procStimCBVData,1)*100;
            stdStimCBVData = std(procStimCBVData,0,1)*100;
            meanStimCortMUAData = mean(procStimCortMUAData,1)*100;
            stdStimCortMUAData = std(procStimCortMUAData,0,1)*100;
            meanStimHipMUAData = mean(procStimHipMUAData,1)*100;
            stdStimHipMUAData = std(procStimHipMUAData,0,1)*100;
            meanStimCortGamData = mean(procStimCortGamData,1)*100;
            stdStimCortGamData = std(procStimCortGamData,0,1)*100;
            meanStimHipGamData = mean(procStimHipGamData,1)*100;
            stdStimHipGamData = std(procStimHipGamData,0,1)*100;
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
                        stimCorticalS_Data = AllSpecData.(stimSpecField).normS{kk,1};
                        stimHippocampalS_Data = AllSpecData.hippocampus.normS{kk,1};
                    end
                end
                stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
                stimStartTimeIndex = stimStartTimeIndex(1);
                stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
                stimDurationIndex = stimDurationIndex(end);
                stimCorticalS_Vals = stimCorticalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                stimHippocampalS_Vals = stimHippocampalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                % mean subtract each row with detrend
                transpStimCorticalS_Vals = stimCorticalS_Vals';
                transpStimHippocampalS_Vals = stimHippocampalS_Vals';
                dTStimCortS_Vals = transpStimCorticalS_Vals;
                dTStimCortS_Vals = dTStimCortS_Vals(1:12*specSamplingRate + 1,:);
                dTStimHipS_Vals = transpStimHippocampalS_Vals;
                dTStimHipS_Vals = dTStimHipS_Vals(1:12*specSamplingRate + 1,:);
                stimCortZhold = cat(3,stimCortZhold,dTStimCortS_Vals');
                stimHipZhold = cat(3,stimHipZhold,dTStimHipS_Vals');
            end
            % figure time/frequency axis and average each S data matrix through time
            meanStimCortS = mean(stimCortZhold,3);
            meanStimHipS = mean(stimHipZhold,3);
            % save results
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).count = size(procStimHipMUAData,1);
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).CBV_HbT.HbT = meanStimHbTData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).CBV_HbT.HbTStD = stdStimHbTData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).CBV.CBV = meanStimCBVData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).CBV.CBVStD = stdStimCBVData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).MUA.corticalData = meanStimCortMUAData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).MUA.corticalStD = stdStimCortMUAData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).MUA.hippocampalData = meanStimHipMUAData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).MUA.hippocampalStD = stdStimHipMUAData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).Gam.corticalData = meanStimCortGamData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).Gam.corticalStD = stdStimCortGamData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).Gam.hippocampalData = meanStimHipGamData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).Gam.hippocampalStD = stdStimHipGamData;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).timeVector = timeVector;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).LFP.corticalS = meanStimCortS;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).LFP.hippocampalS = meanStimHipS;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).LFP.T = T2;
            AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoid).LFP.F = F;
        end
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end

