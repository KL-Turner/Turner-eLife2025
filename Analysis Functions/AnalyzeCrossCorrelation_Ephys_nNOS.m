function [Results_CrossCorr_Ephys] = AnalyzeCrossCorrelation_Ephys(animalID,group,set,rootFolder,delim,Results_CrossCorr_Ephys)
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
% find and load RestingBaselines struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load AllSpecStruct struct
allSpecStructFileStruct = dir('*_AllSpecStructC.mat');
allSpecStructFile = {allSpecStructFileStruct.name}';
allSpecStructFileID = char(allSpecStructFile);
load(allSpecStructFileID,'-mat')
% find and load SleepData struct
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
params.minTime.NREM = 30;
params.minTime.REM = 60;
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'stimDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
hemispheres = {'LH','RH'};
dataTypes = {'deltaBandPower','gammaBandPower','muaPower'};
% go through each valid data type for arousal-based cross-correlation analysis
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for qq = 1:length(dataTypes)
        dataType = dataTypes{1,qq};
        neuralDataType = ['cortical_' hemisphere];
        % lowpass filter
        samplingRate = RestData.HbT.(hemisphere).samplingRate;
        [z,p,k] = butter(4,1/(samplingRate/2),'low');
        [sos,g] = zp2sos(z,p,k);
        trialDuration_sec = RestData.HbT.LH.trialDuration_sec;
        sleepBinWidth = 5;
        lagTime = 5;
        maxLag = lagTime*samplingRate;
        %% Rest
        clear restFiltHbT restFiltNeural
        [restLogical] = FilterEvents_IOS(RestData.HbT.(hemisphere),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.HbT.(hemisphere),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.HbT.(hemisphere).fileIDs(combRestLogical,:);
        restDurations = RestData.HbT.(hemisphere).durations(combRestLogical,:);
        restEventTimes = RestData.HbT.(hemisphere).eventTimes(combRestLogical,:);
        restingHbTData = RestData.HbT.(hemisphere).data(combRestLogical,:);
        restingNeuralData = RestData.(neuralDataType).(dataType).NormData(combRestLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [restFinalRestHbTData,restFinalFileIDs,restFinalDurations,restFinalEventTimes] = RemoveInvalidData_IOS(restingHbTData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [restFinalRestNeuralData,~,~,~] = RemoveInvalidData_IOS(restingNeuralData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        % check whether the event occurs in the appropriate time frame
        for mm =  1:length(restFinalRestNeuralData)
            restFiltHbT{mm,1} = detrend(restFinalRestHbTData{mm,1},'constant');
            restFiltNeural{mm,1} = detrend(restFinalRestNeuralData{mm,1},'constant');
        end
        % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
        % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
        % set parameters for cross-correlation analysis
        % run cross-correlation analysis - average through time
        for dd = 1:length(restFiltHbT)
            restHbTarray = restFiltHbT{dd,1};
            restNeuralarray = restFiltNeural{dd,1};
            [restHbTvNeuralxcVals(dd,:),restNeural_lags] = xcorr(restHbTarray,restNeuralarray,maxLag,'coeff');
        end
        restMeanHbTvNeuralxcVals = mean(restHbTvNeuralxcVals,1);
        % save results
        Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Rest.lags = restNeural_lags;
        Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Rest.xcVals = restMeanHbTvNeuralxcVals;
        clear restProcData
        if qq == 1
            cc = 1;
            for bb = 1:length(restFinalFileIDs)
                restFileID = restFinalFileIDs{bb,1};
                % check whether the event occurs in the appropriate time frame
                restStartTime = ceil(restFinalEventTimes(bb,1)*10)/10; % *10/10 used to round to first decimal place in a floor/ceil fashion.
                restDuration = floor(restFinalDurations(bb,1)*10)/10;
                if restStartTime >= 0.5 && (restStartTime + restDuration) <= (trialDuration_sec - 0.5)
                    % remove the number of samples due to rounding up to start and rounding down to end. This is done to keep the HbT/MUA vectores aligned positionally with the upcoming
                    % spectral analysis which is at 10 Hz
                    leadSamples = round((restStartTime - restFinalEventTimes(bb,1))*samplingRate);
                    lagSamples = round((restFinalDurations(bb,1) - restDuration)*samplingRate);
                    % load in CBV_HbT from rest period
                    restHbT = restFinalRestHbTData{bb,1};
                    % remove leading/lag samples due to rounding to nearest 0.1 up/0.1 down
                    restSnipHbT = restHbT(1 + leadSamples:end - lagSamples);
                    restFiltHbT = filtfilt(sos,g,detrend(restSnipHbT,'constant'));
                    % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
                    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
                    if length(restFiltHbT) < params.minTime.Rest*samplingRate
                        restChunkSampleDiff = params.minTime.Rest*samplingRate - length(restFiltHbT);
                        restPadHbT = (ones(1,restChunkSampleDiff))*restFiltHbT(end);
                        restShortHbT = horzcat(restFiltHbT,restPadHbT);
                    else
                        restShortHbT = restFiltHbT(1:params.minTime.Rest*samplingRate);
                    end
                    restProcData.HbT{cc,1} = restShortHbT;
                    % extract LFP from spectrograms associated with the whisking indecies
                    specDataFileID = [animalID '_' restFileID '_SpecDataC.mat'];
                    clear S_data
                    for g = 1:length(AllSpecData.(neuralDataType).fileIDs)
                        if strcmp(AllSpecData.(neuralDataType).fileIDs{g,1},specDataFileID) == true
                            rest_S = AllSpecData.(neuralDataType).normS{g,1};
                            rest_T = round(AllSpecData.(neuralDataType).T{g,1},1);
                            rest_F = AllSpecData.(neuralDataType).F{g,1};
                        end
                    end
                    restStartTimeIndex = find(rest_T == restStartTime);
                    restStartTimeIndex = restStartTimeIndex(1);
                    restDurationIndex = find(rest_T == round((restStartTime + restDuration),1));
                    restDurationIndex = restDurationIndex(1);
                    restS_Vals = rest_S(:,restStartTimeIndex:restDurationIndex);
                    % only take the first min rest time in seconds
                    shortRestS_Vals = restS_Vals(:,1:params.minTime.Rest*samplingRate);
                    % mean subtract with detrend and lowpass filter each column
                    restProcData.S{cc,1} = detrend(shortRestS_Vals','constant')';
                    cc = cc + 1;
                end
                % set parameters for cross-correlation analysis
                restHbTvLFPzhold = [];
                restHbTvLFPxcVals = ones(length(rest_F),2*maxLag + 1);
                % run cross-correlation analysis - average through time
                for dd = 1:length(restProcData.HbT)
                    for ee = 1:size(restProcData.S{dd, 1}, 1)
                        restHbTarray = restProcData.HbT{dd,1};
                        restNeuralArray = restProcData.S{dd,1}(ee,:);
                        [restHbTvLFPxcVals(ee,:),restLFP_lags] = xcorr(restHbTarray,restNeuralArray,maxLag,'coeff');
                    end
                    restHbTvLFPzhold = cat(3,restHbTvLFPzhold,restHbTvLFPxcVals);
                end
                restMeanHbTvLFPxcVals = mean(restHbTvLFPzhold,3);
            end
            % save results
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Rest.lags = restLFP_lags;
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Rest.f = rest_F;
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Rest.xcVals = restMeanHbTvLFPxcVals;
        end
        %% Alert
        zz = 1; alertHbT = []; alertNeural = []; alertFileIDs = []; alertProcData = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,alertDataFileID] = GetFileInfo_IOS(procDataFileID);
            scoringLabels = [];
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(alertDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            scoringLabelsA = scoringLabels(1:60);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsA,'Not Sleep')) >= 48
                load(procDataFileID,'-mat')
                % only run on files with good pupil measurement
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    alertHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.(hemisphere)(1:300*samplingRate),'constant'));
                    alertNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(neuralDataType).(dataType)(1:300*samplingRate),'constant'));
                    alertFileIDs{zz,1} = alertDataFileID;
                    zz = zz + 1;
                end
            end
            scoringLabelsB = scoringLabels(61:120);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsB,'Not Sleep')) >= 48
                load(procDataFileID,'-mat')
                % only run on files with good pupil measurement
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    alertHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
                    alertNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(neuralDataType).(dataType)(300*samplingRate + 1:600*samplingRate),'constant'));
                    alertFileIDs{zz,1} = alertDataFileID;
                    zz = zz + 1;
                end
            end
            scoringLabelsC = scoringLabels(121:180);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsC,'Not Sleep')) >= 48
                load(procDataFileID,'-mat')
                % only run on files with good pupil measurement
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    alertHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.(hemisphere)(600*samplingRate + 1:end),'constant'));
                    alertNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(neuralDataType).(dataType)(600*samplingRate + 1:end),'constant'));
                    alertFileIDs{zz,1} = alertDataFileID;
                    zz = zz + 1;
                end
            end
        end
        if isempty(alertHbT) == false
            % run cross-correlation analysis - average through time
            for dd = 1:length(alertHbT)
                alertHbTarray = alertHbT{dd,1};
                alertNeuralarray = alertNeural{dd,1};
                [alertHbTvNeuralxcVals(dd,:),alertNeural_lags] = xcorr(alertHbTarray,alertNeuralarray,maxLag,'coeff');
            end
            alertMeanHbTvNeuralxcVals = mean(alertHbTvNeuralxcVals,1);
            % save results
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.lags = alertNeural_lags;
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.xcVals = alertMeanHbTvNeuralxcVals;
        else
            % save results
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.lags = [];
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.xcVals = [];
        end
        if qq == 1
            if isempty(alertHbT) == false
                for bb = 1:length(alertFileIDs)
                    alertFileID = alertFileIDs{bb,1};
                    procDataFileID = [animalID '_' alertFileID '_ProcData.mat'];
                    load(procDataFileID)
                    hbtData = ProcData.data.HbT.(hemisphere);
                    alertProcData.HbT{bb,1}  = filtfilt(sos,g,detrend(hbtData(22:end - 22),'constant'));
                    specDataFileID = [animalID '_' alertFileID '_SpecDataC.mat'];
                    clear S_data
                    for g = 1:length(AllSpecData.(neuralDataType).fileIDs)
                        if strcmp(AllSpecData.(neuralDataType).fileIDs{g,1},specDataFileID) == true
                            alert_S = AllSpecData.(neuralDataType).normS{g,1};
                            alert_F = AllSpecData.(neuralDataType).F{g,1};
                        end
                    end
                    % mean subtract with detrend and lowpass filter each column
                    alertProcData.S{bb,1} = detrend(alert_S','constant')';
                end
                % set parameters for cross-correlation analysis
                alertHbTvLFPzhold = [];
                alertHbTvLFPxcVals = ones(length(alert_F),2*maxLag + 1);
                % run cross-correlation analysis - average through time
                for dd = 1:length(alertProcData.HbT)
                    for ee = 1:size(alertProcData.S{dd, 1}, 1)
                        alertHbTarray = alertProcData.HbT{dd,1};
                        alertNeuralArray = alertProcData.S{dd,1}(ee,:);
                        [alertHbTvLFPxcVals(ee,:),alertLFP_lags] = xcorr(alertHbTarray,alertNeuralArray,maxLag,'coeff');
                    end
                    alertHbTvLFPzhold = cat(3,alertHbTvLFPzhold,alertHbTvLFPxcVals);
                end
                alertMeanHbTvLFPxcVals = mean(alertHbTvLFPzhold,3);
                % save results
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Alert.lags = alertLFP_lags;
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Alert.f = alert_F;
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Alert.xcVals = alertMeanHbTvLFPxcVals;
            else
                % save results
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Alert.lags = [];
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Alert.f = [];
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Alert.xcVals = [];
            end
        end
        %% Asleep
        zz = 1;
        asleepHbT = []; asleepNeural = []; asleepFileIDs = []; asleepProcData = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,asleepDataFileID] = GetFileInfo_IOS(procDataFileID);
            scoringLabels = [];
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(asleepDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            scoringLabelsA = scoringLabels(1:60);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsA,'Not Sleep')) <= 12
                load(procDataFileID,'-mat')
                % only run on files with good pupil measurement
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    asleepHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.(hemisphere)(1:300*samplingRate),'constant'));
                    asleepNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(neuralDataType).(dataType)(1:300*samplingRate),'constant'));
                    asleepFileIDs{zz,1} = asleepDataFileID;
                    zz = zz + 1;
                end
            end
            scoringLabelsB = scoringLabels(61:120);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsB,'Not Sleep')) <= 12
                load(procDataFileID,'-mat')
                % only run on files with good pupil measurement
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    asleepHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
                    asleepNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(neuralDataType).(dataType)(300*samplingRate + 1:600*samplingRate),'constant'));
                    asleepFileIDs{zz,1} = asleepDataFileID;
                    zz = zz + 1;
                end
            end
            scoringLabelsC = scoringLabels(121:180);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsC,'Not Sleep')) <= 12
                load(procDataFileID,'-mat')
                % only run on files with good pupil measurement
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    asleepHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.(hemisphere)(600*samplingRate + 1:end),'constant'));
                    asleepNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(neuralDataType).(dataType)(600*samplingRate + 1:end),'constant'));
                    asleepFileIDs{zz,1} = asleepDataFileID;
                    zz = zz + 1;
                end
            end
        end
        if isempty(asleepHbT) == false
            % run cross-correlation analysis - average through time
            for dd = 1:length(asleepHbT)
                asleepHbTarray = asleepHbT{dd,1};
                asleepNeuralarray = asleepNeural{dd,1};
                [asleepHbTvNeuralxcVals(dd,:),asleepNeural_lags] = xcorr(asleepHbTarray,asleepNeuralarray,maxLag,'coeff');
            end
            asleepMeanHbTvNeuralxcVals = mean(asleepHbTvNeuralxcVals,1);
            % save results
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.lags = asleepNeural_lags;
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.xcVals = asleepMeanHbTvNeuralxcVals;
        else
            % save results
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.lags = [];
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.xcVals = [];
        end
        if qq == 1
            if isempty(asleepHbT) == false
                for bb = 1:length(asleepFileIDs)
                    asleepFileID = asleepFileIDs{bb,1};
                    procDataFileID = [animalID '_' asleepFileID '_ProcData.mat'];
                    load(procDataFileID)
                    hbtData = ProcData.data.HbT.(hemisphere);
                    asleepProcData.HbT{bb,1}  = filtfilt(sos,g,detrend(detrend(hbtData(22:end - 22)),'constant'));
                    specDataFileID = [animalID '_' asleepFileID '_SpecDataC.mat'];
                    clear S_data
                    for g = 1:length(AllSpecData.(neuralDataType).fileIDs)
                        if strcmp(AllSpecData.(neuralDataType).fileIDs{g,1},specDataFileID) == true
                            asleep_S = AllSpecData.(neuralDataType).normS{g,1};
                            asleep_F = AllSpecData.(neuralDataType).F{g,1};
                        end
                    end
                    % mean subtract with detrend and lowpass filter each column
                    try
                        asleepProcData.S{bb,1} = detrend(asleep_S','constant')';
                    catch
                        keyboard
                    end
                end
                % set parameters for cross-correlation analysis
                asleepHbTvLFPzhold = [];
                asleepHbTvLFPxcVals = ones(length(asleep_F),2*maxLag + 1);
                % run cross-correlation analysis - average through time
                for dd = 1:length(asleepProcData.HbT)
                    for ee = 1:size(asleepProcData.S{dd, 1}, 1)
                        asleepHbTarray = asleepProcData.HbT{dd,1};
                        asleepNeuralArray = asleepProcData.S{dd,1}(ee,:);
                        [asleepHbTvLFPxcVals(ee,:),asleepLFP_lags] = xcorr(asleepHbTarray,asleepNeuralArray,maxLag,'coeff');
                    end
                    asleepHbTvLFPzhold = cat(3,asleepHbTvLFPzhold,asleepHbTvLFPxcVals);
                end
                asleepMeanHbTvLFPxcVals = mean(asleepHbTvLFPzhold,3);
                % save results
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Asleep.lags = asleepLFP_lags;
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Asleep.f = asleep_F;
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Asleep.xcVals = asleepMeanHbTvLFPxcVals;
            else
                % save results
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Asleep.lags = [];
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Asleep.f = [];
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.Asleep.xcVals = [];
            end
        end
        %% All
        zz = 1;
        allHbT = []; allNeural = []; allFileIDs = []; allProcData = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,allDataFileID] = GetFileInfo_IOS(procDataFileID);
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.(hemisphere)(1:300*samplingRate),'constant'));
                allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(neuralDataType).(dataType)(1:300*samplingRate),'constant'));
                allFileIDs{zz,1} = allDataFileID;
                zz = zz + 1;
                allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
                allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(neuralDataType).(dataType)(300*samplingRate + 1:600*samplingRate),'constant'));
                allFileIDs{zz,1} = allDataFileID;
                zz = zz + 1;
                allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.(hemisphere)(600*samplingRate + 1:end),'constant'));
                allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(neuralDataType).(dataType)(600*samplingRate + 1:end),'constant'));
                allFileIDs{zz,1} = allDataFileID;
                zz = zz + 1;
            end
        end
        if isempty(allHbT) == false
            % run cross-correlation analysis - average through time
            for dd = 1:length(allHbT)
                allHbTarray = allHbT{dd,1};
                allNeuralarray = allNeural{dd,1};
                [allHbTvNeuralxcVals(dd,:),allNeural_lags] = xcorr(allHbTarray,allNeuralarray,maxLag,'coeff');
            end
            allMeanHbTvNeuralxcVals = mean(allHbTvNeuralxcVals,1);
            % save results
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).All.lags = allNeural_lags;
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).All.xcVals = allMeanHbTvNeuralxcVals;
        else
            % save results
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).All.lags = [];
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).All.xcVals = [];
        end
        if qq == 1
            if isempty(allHbT) == false
                for bb = 1:length(allFileIDs)
                    allFileID = allFileIDs{bb,1};
                    procDataFileID = [animalID '_' allFileID '_ProcData.mat'];
                    load(procDataFileID)
                    hbtData = ProcData.data.HbT.(hemisphere);
                    allProcData.HbT{bb,1}  = filtfilt(sos,g,detrend(hbtData(22:end - 22),'constant'));
                    specDataFileID = [animalID '_' allFileID '_SpecDataC.mat'];
                    clear S_data
                    for g = 1:length(AllSpecData.(neuralDataType).fileIDs)
                        if strcmp(AllSpecData.(neuralDataType).fileIDs{g,1},specDataFileID) == true
                            all_S = AllSpecData.(neuralDataType).normS{g,1};
                            all_F = AllSpecData.(neuralDataType).F{g,1};
                        end
                    end
                    % mean subtract with detrend and lowpass filter each column
                    allProcData.S{bb,1} = detrend(all_S','constant')';
                end
                % set parameters for cross-correlation analysis
                allHbTvLFPzhold = [];
                allHbTvLFPxcVals = ones(length(all_F),2*maxLag + 1);
                % run cross-correlation analysis - average through time
                for dd = 1:length(allProcData.HbT)
                    for ee = 1:size(allProcData.S{dd, 1}, 1)
                        allHbTarray = allProcData.HbT{dd,1};
                        allNeuralArray = allProcData.S{dd,1}(ee,:);
                        [allHbTvLFPxcVals(ee,:),allLFP_lags] = xcorr(allHbTarray,allNeuralArray,maxLag,'coeff');
                    end
                    allHbTvLFPzhold = cat(3,allHbTvLFPzhold,allHbTvLFPxcVals);
                end
                allMeanHbTvLFPxcVals = mean(allHbTvLFPzhold,3);
                % save results
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.All.lags = allLFP_lags;
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.All.f = all_F;
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.All.xcVals = allMeanHbTvLFPxcVals;
            else
                % save results
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.All.lags = [];
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.All.f = [];
                Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.All.xcVals = [];
            end
        end
        %% NREM
        [NREM_finalHbT,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.HbT.(hemisphere),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [NREM_finalNeural,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(neuralDataType).(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        % adjust[HbT] and Neural events to match the edits made to the length of each spectrogram
        for mm = 1:length(NREM_finalHbT)
            NREM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(NREM_finalHbT{mm,1},'constant'));
            NREM_finalNeuralVals{mm,1} = filtfilt(sos,g,detrend(NREM_finalNeural{mm,1},'constant'));
        end
        % run cross-correlation analysis - average through time
        for nn = 1:length(NREM_finalHbTVals)
            NREM_HbT_array = NREM_finalHbTVals{nn,1};
            NREM_Neural_array = NREM_finalNeuralVals{nn,1};
            [NREM_HbTvNeuralxcVals(nn,:),NREM_Neural_lags] = xcorr(NREM_HbT_array,NREM_Neural_array,maxLag,'coeff');
        end
        NREM_meanHbTvNeuralxcVals = mean(NREM_HbTvNeuralxcVals,1);
        % save results
        Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).NREM.lags = NREM_Neural_lags;
        Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).NREM.xcVals = NREM_meanHbTvNeuralxcVals;
        clear NREM_finalHbTVals NREM_sleepNeuralVals NREM_finalSleepNeuralVals
        if qq == 1
            NREM_sleepTime = params.minTime.NREM;   % seconds
            [NREM_finalHbT,NREM_allSleepFileIDs,NREM_finalBinTimes] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.HbT.(hemisphere),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            NREM_uniqueSleepFileIDs = unique(NREM_allSleepFileIDs);
            jj = 1;
            for ff = 1:length(NREM_uniqueSleepFileIDs)
                % pull out the bin times (there may be multiple events) in each unique NREM sleep file
                NREM_uniqueSleepFileID = char(NREM_uniqueSleepFileIDs(ff));
                hh = 1;
                clear NREM_binTimes
                for gg = 1:length(NREM_allSleepFileIDs)
                    NREM_sleepFileID = char(NREM_allSleepFileIDs(gg));
                    if strcmp(NREM_uniqueSleepFileID,NREM_sleepFileID)
                        NREM_binTimes{hh,1} = NREM_finalBinTimes{gg,1};
                        hh = hh + 1;
                    end
                end
                % pull out the Spectrogram data that matches the unique NREM sleep file
                NREM_specDataFileID = [animalID '_' NREM_uniqueSleepFileID '_SpecDataC.mat'];
                load(NREM_specDataFileID,'-mat')
                NREM_S_Data = SpecData.(neuralDataType).normS;
                for ii = 1:length(NREM_binTimes)
                    NREM_Bins = NREM_binTimes{ii,1};
                    NREM_startTime = NREM_Bins(1) - sleepBinWidth;
                    NREM_endTime = NREM_Bins(end);
                    if NREM_startTime > 5 && NREM_endTime < trialDuration_sec
                        NREM_startTimeIndex = find(rest_T == NREM_startTime);
                        NREM_startTimeIndex = NREM_startTimeIndex(1);
                        NREM_durationIndex = find(rest_T == NREM_endTime);
                        NREM_durationIndex = NREM_durationIndex(1);
                        NREM_sleepNeuralVals{jj,1} = NREM_S_Data(:,NREM_startTimeIndex:NREM_durationIndex);
                        editIndex{jj,1} = {'none'};
                    elseif NREM_startTime == 5 && length(NREM_Bins) >= (params.minTime.NREM/sleepBinWidth + 1)
                        NREM_startTime = NREM_Bins(2) - sleepBinWidth;
                        NREM_endTime = NREM_Bins(end);
                        NREM_startTimeIndex = find(rest_T == NREM_startTime);
                        NREM_startTimeIndex = NREM_startTimeIndex(1);
                        NREM_durationIndex = find(rest_T == NREM_endTime);
                        NREM_durationIndex = NREM_durationIndex(1);
                        NREM_sleepNeuralVals{jj,1} = NREM_S_Data(:,NREM_startTimeIndex:NREM_durationIndex);
                        editIndex{jj,1} = {'leading'};
                    elseif NREM_endTime == 900 && length(NREM_Bins) >= (params.minTime.NREM/sleepBinWidth + 1)
                        NREM_startTime = NREM_Bins(1) - sleepBinWidth;
                        NREM_endTime = NREM_Bins(end - 1);
                        NREM_startTimeIndex = find(rest_T == NREM_startTime);
                        NREM_startTimeIndex = NREM_startTimeIndex(1);
                        NREM_durationIndex = find(rest_T == NREM_endTime);
                        NREM_durationIndex = NREM_durationIndex(1);
                        NREM_sleepNeuralVals{jj,1} = NREM_S_Data(:,NREM_startTimeIndex:NREM_durationIndex);
                        editIndex{jj,1} = {'lagging'};
                    else
                        NREM_sleepNeuralVals{jj,1} = [];
                        editIndex{jj,1} = {'delete'};
                    end
                    jj = jj + 1;
                end
            end
            % detrend spectrogram neural values
            for kk = 1:length(NREM_sleepNeuralVals)
                NREM_indSleepNeuralVals = NREM_sleepNeuralVals{kk,1};
                if isempty(NREM_indSleepNeuralVals) == false
                    NREM_indSleepNeuralVals = NREM_indSleepNeuralVals(:,1:NREM_sleepTime*samplingRate)';
                    NREM_dtSleepNeuralVals{kk,1} = detrend(NREM_indSleepNeuralVals,'constant')';
                else
                    NREM_dtSleepNeuralVals{kk,1} = [];
                end
            end
            % adjust[HbT] and MUA events to match the edits made to the length of each spectrogram
            mm = 1;
            for ll = 1:length(NREM_dtSleepNeuralVals)
                if isempty(NREM_dtSleepNeuralVals) == false
                    NREM_finalSleepNeuralVals{mm,1} = NREM_dtSleepNeuralVals{ll,1};
                    if strcmp(editIndex{ll,1},'none') == true
                        NREM_HbTVals = NREM_finalHbT{ll,1}(1:NREM_sleepTime*samplingRate);
                        NREM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(NREM_HbTVals,'constant'));
                        mm = mm + 1;
                    elseif strcmp(editIndex{ll,1},'leading') == true
                        NREM_HbTVals = NREM_finalHbT{ll,1}((samplingRate*sleepBinWidth) + 1:(NREM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                        NREM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(NREM_HbTVals,'constant'));
                        mm = mm + 1;
                    elseif strcmp(editIndex{ll,1},'lagging') == true
                        NREM_HbTVals = NREM_finalHbT{ll,1}(1:NREM_sleepTime*samplingRate);
                        NREM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(NREM_HbTVals,'constant'));
                        mm = mm + 1;
                    elseif strcmp(editIndex{ll,1},'delete') == true
                        % remove HbT/MUA from final file
                    end
                end
            end
            % run cross-correlation analysis - average through time
            NREM_F = SpecData.(neuralDataType).F;
            NREM_HbTvLFPzHold = [];
            NREM_HbTvLFPxcVals = ones(size(NREM_indSleepNeuralVals,2),2*maxLag + 1);
            for nn = 1:length(NREM_finalSleepNeuralVals)
                for oo = 1:size(NREM_finalSleepNeuralVals{nn,1},1)
                    NREM_HbT_array = NREM_finalHbTVals{nn,1};
                    NREM_Neural_array = NREM_finalSleepNeuralVals{nn,1}(oo,:);
                    [NREM_HbTvLFPxcVals(oo,:),NREM_LFP_lags] = xcorr(NREM_HbT_array,NREM_Neural_array,maxLag,'coeff');
                end
                NREM_HbTvLFPzHold = cat(3,NREM_HbTvLFPzHold,NREM_HbTvLFPxcVals);
            end
            NREM_meanHbTvLFPxcVals = mean(NREM_HbTvLFPzHold,3);
            % save results
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.NREM.lags = NREM_LFP_lags;
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.NREM.f = NREM_F;
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.NREM.xcVals = NREM_meanHbTvLFPxcVals;
        end
        %% REM
        [REM_finalHbT,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.HbT.(hemisphere),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [REM_finalNeural,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(neuralDataType).(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        % adjust[HbT] and Neural events to match the edits made to the length of each spectrogram
        for mm = 1:length(REM_finalHbT)
            REM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(REM_finalHbT{mm,1},'constant'));
            REM_finalNeuralVals{mm,1} = filtfilt(sos,g,detrend(REM_finalNeural{mm,1},'constant'));
        end
        % run cross-correlation analysis - average through time
        for nn = 1:length(REM_finalHbTVals)
            REM_HbT_array = REM_finalHbTVals{nn,1};
            REM_Neural_array = REM_finalNeuralVals{nn,1};
            [REM_HbTvNeuralxcVals(nn,:),REM_Neural_lags] = xcorr(REM_HbT_array,REM_Neural_array,maxLag,'coeff');
        end
        REM_meanHbTvNeuralxcVals = mean(REM_HbTvNeuralxcVals,1);
        % save results
        Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).REM.lags = REM_Neural_lags;
        Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).(dataType).REM.xcVals = REM_meanHbTvNeuralxcVals;
        clear REM_finalHbTVals REM_sleepNeuralVals REM_finalSleepNeuralVals
        if qq == 1
            REM_sleepTime = params.minTime.REM;   % seconds
            [REM_finalHbT,REM_allSleepFileIDs,REM_finalBinTimes] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.HbT.(hemisphere),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            REM_uniqueSleepFileIDs = unique(REM_allSleepFileIDs);
            jj = 1;
            for ff = 1:length(REM_uniqueSleepFileIDs)
                % pull out the bin times (there may be multiple events) in each unique REM sleep file
                REM_uniqueSleepFileID = char(REM_uniqueSleepFileIDs(ff));
                hh = 1;
                clear REM_binTimes
                for gg = 1:length(REM_allSleepFileIDs)
                    REM_sleepFileID = char(REM_allSleepFileIDs(gg));
                    if strcmp(REM_uniqueSleepFileID,REM_sleepFileID)
                        REM_binTimes{hh,1} = REM_finalBinTimes{gg,1};
                        hh = hh + 1;
                    end
                end
                % pull out the Spectrogram data that matches the unique REM sleep file
                REM_specDataFileID = [animalID '_' REM_uniqueSleepFileID '_SpecDataC.mat'];
                load(REM_specDataFileID,'-mat')
                REM_S_Data = SpecData.(neuralDataType).normS;
                for ii = 1:length(REM_binTimes)
                    REM_Bins = REM_binTimes{ii,1};
                    REM_startTime = REM_Bins(1) - sleepBinWidth;
                    REM_endTime = REM_Bins(end);
                    if REM_startTime > 5 && REM_endTime < trialDuration_sec
                        REM_startTimeIndex = find(rest_T == REM_startTime);
                        REM_startTimeIndex = REM_startTimeIndex(1);
                        REM_durationIndex = find(rest_T == REM_endTime);
                        REM_durationIndex = REM_durationIndex(1);
                        REM_sleepNeuralVals{jj,1} = REM_S_Data(:,REM_startTimeIndex:REM_durationIndex);
                        editIndex{jj,1} = {'none'};
                    elseif REM_startTime == 5 && length(REM_Bins) >= (params.minTime.REM/sleepBinWidth + 1)
                        REM_startTime = REM_Bins(2) - sleepBinWidth;
                        REM_endTime = REM_Bins(end);
                        REM_startTimeIndex = find(rest_T == REM_startTime);
                        REM_startTimeIndex = REM_startTimeIndex(1);
                        REM_durationIndex = find(rest_T == REM_endTime);
                        REM_durationIndex = REM_durationIndex(1);
                        REM_sleepNeuralVals{jj,1} = REM_S_Data(:,REM_startTimeIndex:REM_durationIndex);
                        editIndex{jj,1} = {'leading'};
                    elseif REM_endTime == 900 && length(REM_Bins) >= (params.minTime.REM/sleepBinWidth + 1)
                        REM_startTime = REM_Bins(1) - sleepBinWidth;
                        REM_endTime = REM_Bins(end - 1);
                        REM_startTimeIndex = find(rest_T == REM_startTime);
                        REM_startTimeIndex = REM_startTimeIndex(1);
                        REM_durationIndex = find(rest_T == REM_endTime);
                        REM_durationIndex = REM_durationIndex(1);
                        REM_sleepNeuralVals{jj,1} = REM_S_Data(:,REM_startTimeIndex:REM_durationIndex);
                        editIndex{jj,1} = {'lagging'};
                    else
                        REM_sleepNeuralVals{jj,1} = [];
                        editIndex{jj,1} = {'delete'};
                    end
                    jj = jj + 1;
                end
            end
            % detrend spectrogram neural values
            for kk = 1:length(REM_sleepNeuralVals)
                REM_indSleepNeuralVals = REM_sleepNeuralVals{kk,1};
                if isempty(REM_indSleepNeuralVals) == false
                    REM_indSleepNeuralVals = REM_indSleepNeuralVals(:,1:REM_sleepTime*samplingRate)';
                    REM_dtSleepNeuralVals{kk,1} = detrend(REM_indSleepNeuralVals,'constant')';
                else
                    REM_dtSleepNeuralVals{kk,1} = [];
                end
            end
            % adjust[HbT] and MUA events to match the edits made to the length of each spectrogram
            mm = 1;
            for ll = 1:length(REM_dtSleepNeuralVals)
                if isempty(REM_dtSleepNeuralVals) == false
                    REM_finalSleepNeuralVals{mm,1} = REM_dtSleepNeuralVals{ll,1};
                    if strcmp(editIndex{ll,1},'none') == true
                        REM_HbTVals = REM_finalHbT{ll,1}(1:REM_sleepTime*samplingRate);
                        REM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(REM_HbTVals,'constant'));
                        mm = mm + 1;
                    elseif strcmp(editIndex{ll,1},'leading') == true
                        REM_HbTVals = REM_finalHbT{ll,1}((samplingRate*sleepBinWidth) + 1:(REM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                        REM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(REM_HbTVals,'constant'));
                        mm = mm + 1;
                    elseif strcmp(editIndex{ll,1},'lagging') == true
                        REM_HbTVals = REM_finalHbT{ll,1}(1:REM_sleepTime*samplingRate);
                        REM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(REM_HbTVals,'constant'));
                        mm = mm + 1;
                    elseif strcmp(editIndex{ll,1},'delete') == true
                        % remove HbT/MUA from final file
                    end
                end
            end
            % run cross-correlation analysis - average through time
            REM_F = SpecData.(neuralDataType).F;
            REM_HbTvLFPzHold = [];
            REM_HbTvLFPxcVals = ones(size(REM_indSleepNeuralVals,2),2*maxLag + 1);
            for nn = 1:length(REM_finalSleepNeuralVals)
                for oo = 1:size(REM_finalSleepNeuralVals{nn,1},1)
                    REM_HbT_array = REM_finalHbTVals{nn,1};
                    REM_Neural_array = REM_finalSleepNeuralVals{nn,1}(oo,:);
                    [REM_HbTvLFPxcVals(oo,:),REM_LFP_lags] = xcorr(REM_HbT_array,REM_Neural_array,maxLag,'coeff');
                end
                REM_HbTvLFPzHold = cat(3,REM_HbTvLFPzHold,REM_HbTvLFPxcVals);
            end
            REM_meanHbTvLFPxcVals = mean(REM_HbTvLFPzHold,3);
            % save results
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.REM.lags = REM_LFP_lags;
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.REM.f = REM_F;
            Results_CrossCorr_Ephys.(group).(animalID).(hemisphere).LFP.REM.xcVals = REM_meanHbTvLFPxcVals;
        end
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_CrossCorr_Ephys.mat','Results_CrossCorr_Ephys')
cd([rootFolder delim 'Data'])