function [Results_BilatCoher_GCaMP] = AnalyzeBilateralCoherence_GCaMP_eLife2025(animalID,group,set,rootFolder,delim,Results_BilatCoher_GCaMP)
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
RestStimCriteria.Fieldname = {'stimDistances'};
RestStimCriteria.Comparison = {'gt'};
RestStimCriteria.Value = {5};
% loop variables
dataTypes = {'HbT','HbO','HbR','GCaMP'};
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% Rest
    clear LH_finalRestData RH_finalRestData LH_ProcRestData RH_ProcRestData LH_restData RH_restData
    clear fLH_finalRestData fRH_finalRestData fLH_ProcRestData fRH_ProcRestData fLH_restData fRH_restData
    samplingRate = RestData.(dataType).LH.samplingRate;
    [restLogical] = FilterEvents_IOS_eLife2025(RestData.(dataType).LH,RestCriteria);
    [stimLogical] = FilterEvents_IOS_eLife2025(RestData.(dataType).LH,RestStimCriteria);
    combRestLogical = logical(restLogical.*stimLogical);
    restFileIDs = RestData.(dataType).LH.fileIDs(combRestLogical,:);
    restEventTimes = RestData.(dataType).LH.eventTimes(combRestLogical,:);
    restDurations = RestData.(dataType).LH.durations(combRestLogical,:);
    LH_RestingData = RestData.(dataType).LH.data(combRestLogical,:);
    RH_RestingData = RestData.(dataType).RH.data(combRestLogical,:);
    fLH_RestingData = RestData.(dataType).fLH.data(combRestLogical,:);
    fRH_RestingData = RestData.(dataType).fRH.data(combRestLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2025(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2025(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [fLH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2025(fLH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [fRH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2025(fRH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    % detrend and truncate data to minimum length to match events
    for bb = 1:length(LH_finalRestData)
        if length(LH_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
            restChunkSampleDiff = params.minTime.Rest*samplingRate - length(LH_finalRestData{bb,1});
            LH_restPad = (ones(1,restChunkSampleDiff))*LH_finalRestData{bb,1}(end);
            RH_restPad = (ones(1,restChunkSampleDiff))*RH_finalRestData{bb,1}(end);
            LH_ProcRestData{bb,1} = horzcat(LH_finalRestData{bb,1},LH_restPad);
            RH_ProcRestData{bb,1} = horzcat(RH_finalRestData{bb,1},RH_restPad);
            LH_ProcRestData{bb,1} = detrend(LH_ProcRestData{bb,1},'constant');
            RH_ProcRestData{bb,1} = detrend(RH_ProcRestData{bb,1},'constant');
            fLH_ProcRestData{bb,1} = horzcat(fLH_finalRestData{bb,1},LH_restPad);
            fRH_ProcRestData{bb,1} = horzcat(fRH_finalRestData{bb,1},RH_restPad);
            fLH_ProcRestData{bb,1} = detrend(fLH_ProcRestData{bb,1},'constant');
            fRH_ProcRestData{bb,1} = detrend(fRH_ProcRestData{bb,1},'constant');
        else
            LH_ProcRestData{bb,1} = detrend(LH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            RH_ProcRestData{bb,1} = detrend(RH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            fLH_ProcRestData{bb,1} = detrend(fLH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            fRH_ProcRestData{bb,1} = detrend(fRH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
        end
    end
    % pre-allocate coherence matrix
    LH_restData = zeros(length(LH_ProcRestData{1,1}),length(LH_ProcRestData));
    RH_restData = zeros(length(RH_ProcRestData{1,1}),length(RH_ProcRestData));
    fLH_restData = zeros(length(fLH_ProcRestData{1,1}),length(fLH_ProcRestData));
    fRH_restData = zeros(length(fRH_ProcRestData{1,1}),length(fRH_ProcRestData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for cc = 1:length(LH_ProcRestData)
        LH_restData(:,cc) = LH_ProcRestData{cc,1};
        RH_restData(:,cc) = RH_ProcRestData{cc,1};
        fLH_restData(:,cc) = fLH_ProcRestData{cc,1};
        fRH_restData(:,cc) = fRH_ProcRestData{cc,1};
    end
    % parameters for coherencyc - information available in function
    params.tapers = [1,1]; % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,1]; % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2,0.05];
    % calculate the coherence between desired signals
    [C_RestData,phi_RestData,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc(LH_restData,RH_restData,params);
    % save results
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.C = C_RestData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.phi = phi_RestData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.f = f_RestData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.confC = confC_RestData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.cErr = cErr_RestData;
    % calculate the coherence between desired signals
    [fC_RestData,fphi_RestData,~,~,~,ff_RestData,fconfC_RestData,~,fcErr_RestData] = coherencyc(fLH_restData,fRH_restData,params);
    % save results
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.fC = fC_RestData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.fphi = fphi_RestData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.ff = ff_RestData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.fconfC = fconfC_RestData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Rest.fcErr = fcErr_RestData;
    %% Alert
    clear LH_AlertData LH_ProcAlertData LH_alertData RH_AlertData RH_ProcAlertData RH_alertData scoringLabels
    clear fLH_AlertData fLH_ProcAlertData fLH_alertData fRH_AlertData fRH_ProcAlertData fRH_alertData
    LH_AlertData = []; % for loop pre-allocation
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,~,allDataFileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
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
    % detrend data
    if isempty(LH_AlertData) == false
        for bb = 1:length(LH_AlertData)
            LH_ProcAlertData{bb,1} = detrend(LH_AlertData{bb,1},'constant');
            RH_ProcAlertData{bb,1} = detrend(RH_AlertData{bb,1},'constant');
            fLH_ProcAlertData{bb,1} = detrend(fLH_AlertData{bb,1},'constant');
            fRH_ProcAlertData{bb,1} = detrend(fRH_AlertData{bb,1},'constant');
        end
        % preallocate coherence matrix
        LH_alertData = zeros(length(LH_ProcAlertData{1,1}),length(LH_ProcAlertData));
        RH_alertData = zeros(length(RH_ProcAlertData{1,1}),length(RH_ProcAlertData));
        fLH_alertData = zeros(length(fLH_ProcAlertData{1,1}),length(fLH_ProcAlertData));
        fRH_alertData = zeros(length(fRH_ProcAlertData{1,1}),length(fRH_ProcAlertData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(LH_ProcAlertData)
            LH_alertData(:,cc) = LH_ProcAlertData{cc,1};
            RH_alertData(:,cc) = RH_ProcAlertData{cc,1};
            fLH_alertData(:,cc) = fLH_ProcAlertData{cc,1};
            fRH_alertData(:,cc) = fRH_ProcAlertData{cc,1};
        end
        % calculate the coherence between desired signals
        params.tapers = [7,13]; % Tapers [n, 2n - 1]
        [C_AlertData,phi_AlertData,~,~,~,f_AlertData,confC_AlertData,~,cErr_AlertData] = coherencyc(LH_alertData,RH_alertData,params);
        % save results
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.C = C_AlertData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.phi = phi_AlertData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.f = f_AlertData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.confC = confC_AlertData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.cErr = cErr_AlertData;
        % calculate the coherence between desired signals
        [fC_AlertData,fphi_AlertData,~,~,~,ff_AlertData,fconfC_AlertData,~,fcErr_AlertData] = coherencyc(fLH_alertData,fRH_alertData,params);
        % save results
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.fC = fC_AlertData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.fphi = fphi_AlertData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.ff = ff_AlertData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.fconfC = fconfC_AlertData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.fcErr = fcErr_AlertData;
    else
        % save results
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.C = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.f = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.confC = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.cErr = [];
        % save results
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.fC = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.ff = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.fconfC = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Alert.fcErr = [];
    end
    %% Asleep
    clear LH_AsleepData LH_ProcAsleepData LH_asleepData RH_AsleepData RH_ProcAsleepData RH_asleepData scoringLabels
    clear fLH_AsleepData fLH_ProcAsleepData fLH_asleepData fRH_AsleepData fRH_ProcAsleepData fRH_asleepData
    LH_AsleepData = []; % for loop pre-allocation
    zz = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,~,allDataFileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
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
    % detrend data
    if isempty(LH_AsleepData) == false
        for bb = 1:length(LH_AsleepData)
            LH_ProcAsleepData{bb,1} = detrend(LH_AsleepData{bb,1},'constant');
            RH_ProcAsleepData{bb,1} = detrend(RH_AsleepData{bb,1},'constant');
            fLH_ProcAsleepData{bb,1} = detrend(fLH_AsleepData{bb,1},'constant');
            fRH_ProcAsleepData{bb,1} = detrend(fRH_AsleepData{bb,1},'constant');
        end
        % preallocate coherence matrix
        LH_asleepData = zeros(length(LH_ProcAsleepData{1,1}),length(LH_ProcAsleepData));
        RH_asleepData = zeros(length(RH_ProcAsleepData{1,1}),length(RH_ProcAsleepData));
        fLH_asleepData = zeros(length(fLH_ProcAsleepData{1,1}),length(fLH_ProcAsleepData));
        fRH_asleepData = zeros(length(fRH_ProcAsleepData{1,1}),length(fRH_ProcAsleepData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(LH_ProcAsleepData)
            LH_asleepData(:,cc) = LH_ProcAsleepData{cc,1};
            RH_asleepData(:,cc) = RH_ProcAsleepData{cc,1};
            fLH_asleepData(:,cc) = fLH_ProcAsleepData{cc,1};
            fRH_asleepData(:,cc) = fRH_ProcAsleepData{cc,1};
        end
        % calculate the coherence between desired signals
        params.tapers = [7,13]; % Tapers [n, 2n - 1]
        [C_AsleepData,phi_AsleepData,~,~,~,f_AsleepData,confC_AsleepData,~,cErr_AsleepData] = coherencyc(LH_asleepData,RH_asleepData,params);
        % save results
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.C = C_AsleepData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.phi = phi_AsleepData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.f = f_AsleepData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.confC = confC_AsleepData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.cErr = cErr_AsleepData;
        % calculate the coherence between desired signals
        [fC_AsleepData,fphi_AsleepData,~,~,~,ff_AsleepData,fconfC_AsleepData,~,fcErr_AsleepData] = coherencyc(fLH_asleepData,fRH_asleepData,params);
        % save results
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.fC = fC_AsleepData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.fphi = fphi_AsleepData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.ff = ff_AsleepData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.fconfC = fconfC_AsleepData;
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.fcErr = fcErr_AsleepData;
    else
        % save results
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.C = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.f = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.confC = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.cErr = [];
        % save results
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.fC = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.ff = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.fconfC = [];
        Results_BilatCoher_GCaMP.(group).(animalID).(dataType).Asleep.fcErr = [];
    end
    %% All
    clear LH_AllData LH_ProcAllData LH_allData RH_AllData RH_ProcAllData RH_allData
    clear fLH_AllData fLH_ProcAllData fLH_allData fRH_AllData fRH_ProcAllData fRH_allData
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
    % detrend data
    for cc = 1:length(LH_AllData)
        LH_ProcAllData{cc,1} = detrend(LH_AllData{cc,1},'constant');
        RH_ProcAllData{cc,1} = detrend(RH_AllData{cc,1},'constant');
        fLH_ProcAllData{cc,1} = detrend(fLH_AllData{cc,1},'constant');
        fRH_ProcAllData{cc,1} = detrend(fRH_AllData{cc,1},'constant');
    end
    % preallocate coherence matrix
    LH_allData = zeros(length(LH_ProcAllData{1,1}),length(LH_ProcAllData));
    RH_allData = zeros(length(RH_ProcAllData{1,1}),length(RH_ProcAllData));
    fLH_allData = zeros(length(fLH_ProcAllData{1,1}),length(fLH_ProcAllData));
    fRH_allData = zeros(length(fRH_ProcAllData{1,1}),length(fRH_ProcAllData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for cc = 1:length(LH_ProcAllData)
        LH_allData(:,cc) = LH_ProcAllData{cc,1};
        RH_allData(:,cc) = RH_ProcAllData{cc,1};
        fLH_allData(:,cc) = fLH_ProcAllData{cc,1};
        fRH_allData(:,cc) = fRH_ProcAllData{cc,1};
    end
    % calculate the coherence between desired signals
    params.tapers = [7,13]; % Tapers [n, 2n - 1]
    [C_AllData,phi_AllData,~,~,~,f_AllData,confC_AllData,~,cErr_AllData] = coherencyc(LH_allData,RH_allData,params);
    % save results
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.C = C_AllData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.phi = phi_AllData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.f = f_AllData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.confC = confC_AllData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.cErr = cErr_AllData;
    % calculate the coherence between desired signals
    [fC_AllData,fphi_AllData,~,~,~,ff_AllData,fconfC_AllData,~,fcErr_AllData] = coherencyc(fLH_allData,fRH_allData,params);
    % save results
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.fC = fC_AllData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.fphi = fphi_AllData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.ff = ff_AllData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.fconfC = fconfC_AllData;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).All.fcErr = fcErr_AllData;
    %% NREM
    clear LH_nremData LH_nrem RH_nremData RH_nrem
    clear fLH_nremData fLH_nrem fRH_nremData fRH_nrem
    % pull data from SleepData.mat structure
    [LH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [RH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [fLH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).fLH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [fRH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).fRH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % detrend and truncate data to minimum length to match events
    for ee = 1:length(LH_nremData)
        LH_nremData{ee,1} = detrend(LH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{ee,1} = detrend(RH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        fLH_nremData{ee,1} = detrend(fLH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        fRH_nremData{ee,1} = detrend(fRH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
    end
    % pre-allocate coherence matrix
    LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
    RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
    fLH_nrem = zeros(length(fLH_nremData{1,1}),length(fLH_nremData));
    fRH_nrem = zeros(length(fRH_nremData{1,1}),length(fRH_nremData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for ff = 1:length(LH_nremData)
        LH_nrem(:,ff) = LH_nremData{ff,1};
        RH_nrem(:,ff) = RH_nremData{ff,1};
        fLH_nrem(:,ff) = fLH_nremData{ff,1};
        fRH_nrem(:,ff) = fRH_nremData{ff,1};
    end
    % calculate the coherence between desired signals
    params.tapers = [3,5]; % Tapers [n, 2n - 1]
    [C_nrem,phi_nrem,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc(LH_nrem,RH_nrem,params);
    % save results
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.C = C_nrem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.phi = phi_nrem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.f = f_nrem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.confC = confC_nrem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.cErr = cErr_nrem;
    % calculate the coherence between desired signals
    [fC_nrem,fphi_nrem,~,~,~,ff_nrem,fconfC_nrem,~,fcErr_nrem] = coherencyc(fLH_nrem,fRH_nrem,params);
    % save results
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.fC = fC_nrem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.fphi = fphi_nrem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.ff = ff_nrem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.fconfC = fconfC_nrem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).NREM.fcErr = fcErr_nrem;
    %% REM
    clear LH_remData LH_rem RH_remData RH_rem
    clear fLH_remData fLH_rem fRH_remData fRH_rem
    % pull data from SleepData.mat structure
    [LH_remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [RH_remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [fLH_remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).fLH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [fRH_remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).fRH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % detrend and truncate data to minimum length to match events
    for ee = 1:length(LH_remData)
        LH_remData{ee,1} = detrend(LH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
        RH_remData{ee,1} = detrend(RH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
        fLH_remData{ee,1} = detrend(fLH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
        fRH_remData{ee,1} = detrend(fRH_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
    end
    % pre-allocate coherence matrix
    LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
    RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
    fLH_rem = zeros(length(fLH_remData{1,1}),length(fLH_remData));
    fRH_rem = zeros(length(fRH_remData{1,1}),length(fRH_remData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for ff = 1:length(LH_remData)
        LH_rem(:,ff) = LH_remData{ff,1};
        RH_rem(:,ff) = RH_remData{ff,1};
        fLH_rem(:,ff) = fLH_remData{ff,1};
        fRH_rem(:,ff) = fRH_remData{ff,1};
    end
    % calculate the coherence between desired signals
    params.tapers = [5,9]; % Tapers [n, 2n - 1]
    [C_rem,phi_rem,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc(LH_rem,RH_rem,params);
    % save results
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.C = C_rem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.phi = phi_rem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.f = f_rem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.confC = confC_rem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.cErr = cErr_rem;
    % calculate the coherence between desired signals
    [fC_rem,fphi_rem,~,~,~,ff_rem,fconfC_rem,~,fcErr_rem] = coherencyc(fLH_rem,fRH_rem,params);
    % save results
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.fC = fC_rem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.fphi = fphi_rem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.ff = ff_rem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.fconfC = fconfC_rem;
    Results_BilatCoher_GCaMP.(group).(animalID).(dataType).REM.fcErr = fcErr_rem;
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_BilatCoher_GCaMP.mat','Results_BilatCoher_GCaMP')
cd([rootFolder delim 'Data'])