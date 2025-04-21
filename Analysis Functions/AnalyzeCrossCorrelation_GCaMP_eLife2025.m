function [Results_CrossCorr_GCaMP] = AnalyzeCrossCorrelation_GCaMP_eLife2025(animalID,group,set,rootFolder,delim,Results_CrossCorr_GCaMP)
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
hemispheres = {'LH','RH','fLH','fRH'};
dataTypes = {'HbT','HbO','HbR'};
% go through each valid data type for arousal-based cross-correlation analysis
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for qq = 1:length(dataTypes)
        dataType = dataTypes{1,qq};
        % lowpass filter
        samplingRate = RestData.(dataType).(hemisphere).samplingRate;
        [z,p,k] = butter(4,1/(samplingRate/2),'low');
        [sos,g] = zp2sos(z,p,k);
        lagTime = 5;
        maxLag = lagTime*samplingRate;
        %% Rest
        [restLogical] = FilterEvents_IOS_eLife2025(RestData.(dataType).(hemisphere),RestCriteria);
        [puffLogical] = FilterEvents_IOS_eLife2025(RestData.(dataType).(hemisphere),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.(dataType).(hemisphere).fileIDs(combRestLogical,:);
        restDurations = RestData.(dataType).(hemisphere).durations(combRestLogical,:);
        restEventTimes = RestData.(dataType).(hemisphere).eventTimes(combRestLogical,:);
        restingHbTData = RestData.(dataType).(hemisphere).data(combRestLogical,:);
        restingNeuralData = RestData.GCaMP.(hemisphere).data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [restFinalRestHbTData,~,~,~] = RemoveInvalidData_IOS_eLife2025(restingHbTData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [restFinalRestNeuralData,~,~,~] = RemoveInvalidData_IOS_eLife2025(restingNeuralData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        % check whether the event occurs in the appropriate time frame
        for mm =  1:length(restFinalRestNeuralData)
            restFiltHbT{mm,1} = filtfilt(sos,g,detrend(restFinalRestHbTData{mm,1},'constant'));
            restFiltNeural{mm,1} = filtfilt(sos,g,detrend(restFinalRestNeuralData{mm,1},'constant'));
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
        Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Rest.lags = restNeural_lags;
        Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Rest.xcVals = restMeanHbTvNeuralxcVals;
        %% Alert
        zz = 1; alertHbT = []; alertNeural = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,alertDataFileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
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
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    alertHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(dataType).(hemisphere)(1:300*samplingRate),'constant'));
                    alertNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP.(hemisphere)(1:300*samplingRate),'constant'));
                    zz = zz + 1;
                end
            end
            scoringLabelsB = scoringLabels(61:120);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsB,'Not Sleep')) >= 48
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    alertHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(dataType).(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
                    alertNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP.(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
                    zz = zz + 1;
                end
            end
            scoringLabelsC = scoringLabels(121:180);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsC,'Not Sleep')) >= 48
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    alertHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(dataType).(hemisphere)(600*samplingRate + 1:end),'constant'));
                    alertNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP.(hemisphere)(600*samplingRate + 1:end),'constant'));
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
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.lags = alertNeural_lags;
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.xcVals = alertMeanHbTvNeuralxcVals;
        else
            % save results
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.lags = [];
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.xcVals = [];
        end
        %% Asleep
        zz = 1;
        asleepHbT = []; asleepNeural = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,asleepDataFileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
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
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    asleepHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(dataType).(hemisphere)(1:300*samplingRate),'constant'));
                    asleepNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP.(hemisphere)(1:300*samplingRate),'constant'));
                    zz = zz + 1;
                end
            end
            scoringLabelsB = scoringLabels(61:120);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsB,'Not Sleep')) <= 12
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    asleepHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(dataType).(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
                    asleepNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP.(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
                    zz = zz + 1;
                end
            end
            scoringLabelsC = scoringLabels(121:180);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsC,'Not Sleep')) <= 12
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    asleepHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(dataType).(hemisphere)(600*samplingRate + 1:end),'constant'));
                    asleepNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP.(hemisphere)(600*samplingRate + 1:end),'constant'));
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
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.lags = asleepNeural_lags;
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.xcVals = asleepMeanHbTvNeuralxcVals;
        else
            % save results
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.lags = [];
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.xcVals = [];
        end
        %% All
        zz = 1;
        allHbT = []; allNeural = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(dataType).(hemisphere)(1:300*samplingRate),'constant'));
                allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP.(hemisphere)(1:300*samplingRate),'constant'));
                zz = zz + 1;
                allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(dataType).(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
                allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP.(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
                zz = zz + 1;
                allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.(dataType).(hemisphere)(600*samplingRate + 1:end),'constant'));
                allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP.(hemisphere)(600*samplingRate + 1:end),'constant'));
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
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).All.lags = allNeural_lags;
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).All.xcVals = allMeanHbTvNeuralxcVals;
        else
            % save results
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).All.lags = [];
            Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).All.xcVals = [];
        end
        %% NREM
        [NREM_finalHbT,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).(hemisphere),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [NREM_finalNeural,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.GCaMP.(hemisphere),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
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
        Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).NREM.lags = NREM_Neural_lags;
        Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).NREM.xcVals = NREM_meanHbTvNeuralxcVals;
        %% REM
        [REM_finalHbT,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).(hemisphere),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [REM_finalNeural,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.GCaMP.(hemisphere),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
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
        Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).REM.lags = REM_Neural_lags;
        Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).(dataType).REM.xcVals = REM_meanHbTvNeuralxcVals;
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_CrossCorr_GCaMP.mat','Results_CrossCorr_GCaMP')
cd([rootFolder delim 'Data'])