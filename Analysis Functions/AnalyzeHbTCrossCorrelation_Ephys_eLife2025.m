function [Results_HbTCrossCorr_Ephys] = AnalyzeHbTCrossCorrelation_Ephys_eLife2025(animalID,group,set,rootFolder,delim,Results_HbTCrossCorr_Ephys)
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
% lowpass filter
samplingRate = RestData.HbT.LH.samplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
trialDuration_sec = RestData.HbT.LH.trialDuration_sec;
sleepBinWidth = 5;
lagTime = 5;
maxLag = lagTime*samplingRate;
%% Rest
clear restFiltHbT restFiltNeural
[restLogical] = FilterEvents_IOS_eLife2025(RestData.HbT.LH,RestCriteria);
[puffLogical] = FilterEvents_IOS_eLife2025(RestData.HbT.LH,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFileIDs = RestData.HbT.LH.fileIDs(combRestLogical,:);
restDurations = RestData.HbT.LH.durations(combRestLogical,:);
restEventTimes = RestData.HbT.LH.eventTimes(combRestLogical,:);
restingHbTData = RestData.HbT.LH.data(combRestLogical,:);
restingNeuralData = RestData.HbT.RH.data(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[restFinalRestHbTData,restFinalFileIDs,restFinalDurations,restFinalEventTimes] = RemoveInvalidData_IOS_eLife2025(restingHbTData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
[restFinalRestNeuralData,~,~,~] = RemoveInvalidData_IOS_eLife2025(restingNeuralData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
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
Results_HbTCrossCorr_Ephys.(group).(animalID).Rest.lags = restNeural_lags;
Results_HbTCrossCorr_Ephys.(group).(animalID).Rest.xcVals = restMeanHbTvNeuralxcVals;
clear restProcData
%% Alert
zz = 1; alertHbT = []; alertNeural = []; alertFileIDs = []; alertProcData = [];
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
        % only run on files with good pupil measurement
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            alertHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.LH(1:300*samplingRate),'constant'));
            alertNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.RH(1:300*samplingRate),'constant'));
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
            alertHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.LH(300*samplingRate + 1:600*samplingRate),'constant'));
            alertNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.RH(300*samplingRate + 1:600*samplingRate),'constant'));
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
            alertHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.LH(600*samplingRate + 1:end),'constant'));
            alertNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.RH(600*samplingRate + 1:end),'constant'));
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
    Results_HbTCrossCorr_Ephys.(group).(animalID).Alert.lags = alertNeural_lags;
    Results_HbTCrossCorr_Ephys.(group).(animalID).Alert.xcVals = alertMeanHbTvNeuralxcVals;
else
    % save results
    Results_HbTCrossCorr_Ephys.(group).(animalID).Alert.lags = [];
    Results_HbTCrossCorr_Ephys.(group).(animalID).Alert.xcVals = [];
end
%% Asleep
zz = 1;
asleepHbT = []; asleepNeural = []; asleepFileIDs = []; asleepProcData = [];
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
        % only run on files with good pupil measurement
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            asleepHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.LH(1:300*samplingRate),'constant'));
            asleepNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.RH(1:300*samplingRate),'constant'));
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
            asleepHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.LH(300*samplingRate + 1:600*samplingRate),'constant'));
            asleepNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.RH(300*samplingRate + 1:600*samplingRate),'constant'));
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
            asleepHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.LH(600*samplingRate + 1:end),'constant'));
            asleepNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.RH(600*samplingRate + 1:end),'constant'));
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
    Results_HbTCrossCorr_Ephys.(group).(animalID).Asleep.lags = asleepNeural_lags;
    Results_HbTCrossCorr_Ephys.(group).(animalID).Asleep.xcVals = asleepMeanHbTvNeuralxcVals;
else
    % save results
    Results_HbTCrossCorr_Ephys.(group).(animalID).Asleep.lags = [];
    Results_HbTCrossCorr_Ephys.(group).(animalID).Asleep.xcVals = [];
end
%% All
zz = 1;
allHbT = []; allNeural = []; allFileIDs = []; allProcData = [];
for cc = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(cc,:);
    [~,~,allDataFileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
    load(procDataFileID,'-mat')
    puffs = ProcData.data.stimulations.LPadSol;
    % don't include trials with stimulation
    if isempty(puffs) == true
        allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.LH(1:300*samplingRate),'constant'));
        allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.RH(1:300*samplingRate),'constant'));
        allFileIDs{zz,1} = allDataFileID;
        zz = zz + 1;
        allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.LH(300*samplingRate + 1:600*samplingRate),'constant'));
        allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.RH(300*samplingRate + 1:600*samplingRate),'constant'));
        allFileIDs{zz,1} = allDataFileID;
        zz = zz + 1;
        allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.LH(600*samplingRate + 1:end),'constant'));
        allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.HbT.RH(600*samplingRate + 1:end),'constant'));
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
    Results_HbTCrossCorr_Ephys.(group).(animalID).All.lags = allNeural_lags;
    Results_HbTCrossCorr_Ephys.(group).(animalID).All.xcVals = allMeanHbTvNeuralxcVals;
else
    % save results
    Results_HbTCrossCorr_Ephys.(group).(animalID).All.lags = [];
    Results_HbTCrossCorr_Ephys.(group).(animalID).All.xcVals = [];
end
%% NREM
[NREM_finalHbT,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.HbT.LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
[NREM_finalNeural,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.HbT.RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
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
Results_HbTCrossCorr_Ephys.(group).(animalID).NREM.lags = NREM_Neural_lags;
Results_HbTCrossCorr_Ephys.(group).(animalID).NREM.xcVals = NREM_meanHbTvNeuralxcVals;
%% REM
[REM_finalHbT,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.HbT.LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
[REM_finalNeural,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.HbT.RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
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
Results_HbTCrossCorr_Ephys.(group).(animalID).REM.lags = REM_Neural_lags;
Results_HbTCrossCorr_Ephys.(group).(animalID).REM.xcVals = REM_meanHbTvNeuralxcVals;
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_HbTCrossCorr_Ephys.mat','Results_HbTCrossCorr_Ephys')
cd([rootFolder delim 'Data'])