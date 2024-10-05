function [Results_NeuralHemoCoher_Ephys] = AnalyzeNeuralHemoCoherence_Ephys_nNOS(animalID,group,set,rootFolder,delim,Results_NeuralHemoCoher_Ephys)
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
% loop variables
hemispheres = {'LH','RH'};
dataTypes = {'deltaBandPower','gammaBandPower'};
for zzz = 1:length(hemispheres)
    hemisphere = hemispheres{1,zzz};
    for xx = 1:length(dataTypes)
        dataType = dataTypes{1,xx};
        %% Rest
        samplingRate = RestData.HbT.LH.samplingRate;
        [restLogical] = FilterEvents_IOS_nNOS(RestData.HbT.(hemisphere),RestCriteria);
        [puffLogical] = FilterEvents_IOS_nNOS(RestData.HbT.(hemisphere),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.HbT.(hemisphere).fileIDs(combRestLogical,:);
        restEventTimes = RestData.HbT.(hemisphere).eventTimes(combRestLogical,:);
        restDurations = RestData.HbT.(hemisphere).durations(combRestLogical,:);
        HbT_RestingData = RestData.HbT.(hemisphere).data(combRestLogical,:);
        Neural_RestingData = RestData.(['cortical_' hemisphere]).(dataType).NormData(combRestLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [HbT_finalRestData,~,~,~] = RemoveInvalidData_IOS_nNOS(HbT_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [Neural_finalRestData,~,~,~] = RemoveInvalidData_IOS_nNOS(Neural_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        clear HbT_ProcRestData Neural_ProcRestData
        % filter, detrend, and truncate data to minimum length to match events
        for bb = 1:length(HbT_finalRestData)
            if length(HbT_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(HbT_finalRestData{bb,1});
                HbT_restPad = (ones(1,restChunkSampleDiff))*HbT_finalRestData{bb,1}(end);
                Neural_restPad = (ones(1,restChunkSampleDiff))*Neural_finalRestData{bb,1}(end);
                HbT_ProcRestData{bb,1} = horzcat(HbT_finalRestData{bb,1},HbT_restPad);
                Neural_ProcRestData{bb,1} = horzcat(Neural_finalRestData{bb,1},Neural_restPad);
                HbT_ProcRestData{bb,1} = detrend(HbT_ProcRestData{bb,1},'constant');
                Neural_ProcRestData{bb,1} = detrend(Neural_ProcRestData{bb,1},'constant');
            else
                HbT_ProcRestData{bb,1} = detrend(HbT_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
                Neural_ProcRestData{bb,1} = detrend(Neural_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            end
        end
        % pre-allocate coherence matrix
        HbT_restData = zeros(length(HbT_ProcRestData{1,1}),length(HbT_ProcRestData));
        Neural_restData = zeros(length(Neural_ProcRestData{1,1}),length(Neural_ProcRestData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(HbT_ProcRestData)
            HbT_restData(:,cc) = HbT_ProcRestData{cc,1};
            Neural_restData(:,cc) = Neural_ProcRestData{cc,1};
        end
        % parameters for coherencyc - information available in function
        params.tapers = [1,1]; % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1]; % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc(HbT_restData,Neural_restData,params);
        % save results
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Rest.C = C_RestData;
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Rest.f = f_RestData;
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Rest.confC = confC_RestData;
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Rest.cErr = cErr_RestData;
        %% Alert
        zz = 1;
        clear HbT_AlertData Neural_AlertData HbT_ProcAlertData Neural_ProcAlertData
        HbT_AlertData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS_nNOS(procDataFileID);
            strDay = ConvertDate_IOS_nNOS(allDataFileDate);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            scoringLabelsA = scoringLabels(1:60);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsA,'Not Sleep')) >= 48
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        HbT_AlertData{zz,1} = ProcData.data.HbT.(hemisphere)(1:300*samplingRate);
                        Neural_AlertData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType)(1:300*samplingRate) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
            scoringLabelsB = scoringLabels(61:120);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsB,'Not Sleep')) >= 48
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        HbT_AlertData{zz,1} = ProcData.data.HbT.(hemisphere)(300*samplingRate + 1:600*samplingRate);
                        Neural_AlertData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType)(300*samplingRate + 1:600*samplingRate) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
            scoringLabelsC = scoringLabels(121:180);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsC,'Not Sleep')) >= 48
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        HbT_AlertData{zz,1} = ProcData.data.HbT.(hemisphere)(600*samplingRate + 1:end);
                        Neural_AlertData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType)(600*samplingRate + 1:end) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(HbT_AlertData) == false
            for bb = 1:length(HbT_AlertData)
                HbT_ProcAlertData{bb,1} = detrend(HbT_AlertData{bb,1},'constant');
                Neural_ProcAlertData{bb,1} = detrend(Neural_AlertData{bb,1},'constant');
            end
            % pre-allocate coherence matrix
            HbT_alertData = zeros(length(HbT_ProcAlertData{1,1}),length(HbT_ProcAlertData));
            Neural_alertData = zeros(length(Neural_ProcAlertData{1,1}),length(Neural_ProcAlertData));
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
            for cc = 1:length(HbT_ProcAlertData)
                HbT_alertData(:,cc) = HbT_ProcAlertData{cc,1};
                Neural_alertData(:,cc) = Neural_ProcAlertData{cc,1};
            end
            % calculate the coherence between desired signals
            params.tapers = [7,13]; % Tapers [n, 2n - 1]
            [C_AlertData,~,~,~,~,f_AlertData,confC_AlertData,~,cErr_AlertData] = coherencyc(HbT_alertData,Neural_alertData,params);
            % save results
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.C = C_AlertData;
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.f = f_AlertData;
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.confC = confC_AlertData;
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.cErr = cErr_AlertData;
        else
            % save results
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.C = [];
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.f = [];
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.confC = [];
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.cErr = [];
        end
        %% Asleep
        zz = 1;
        clear HbT_AsleepData Neural_AsleepData HbT_ProcAsleepData Neural_ProcAsleepData
        HbT_AsleepData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS_nNOS(procDataFileID);
            strDay = ConvertDate_IOS_nNOS(allDataFileDate);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            scoringLabelsA = scoringLabels(1:60);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsA,'Not Sleep')) <= 12
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        HbT_AsleepData{zz,1} = ProcData.data.HbT.(hemisphere)(1:300*samplingRate);
                        Neural_AsleepData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType)(1:300*samplingRate) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
            scoringLabelsB = scoringLabels(61:120);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsB,'Not Sleep')) <= 12
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        HbT_AsleepData{zz,1} = ProcData.data.HbT.(hemisphere)(300*samplingRate + 1:600*samplingRate);
                        Neural_AsleepData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType)(300*samplingRate + 1:600*samplingRate) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
            scoringLabelsC = scoringLabels(121:180);
            % check labels to match arousal state
            if sum(strcmp(scoringLabelsC,'Not Sleep')) <= 12
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        HbT_AsleepData{zz,1} = ProcData.data.HbT.(hemisphere)(600*samplingRate + 1:end);
                        Neural_AsleepData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType)(600*samplingRate + 1:end) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(HbT_AsleepData) == false
            for bb = 1:length(HbT_AsleepData)
                HbT_ProcAsleepData{bb,1} = detrend(HbT_AsleepData{bb,1},'constant');
                Neural_ProcAsleepData{bb,1} = detrend(Neural_AsleepData{bb,1},'constant');
            end
            % pre-allocate coherence matrix
            HbT_asleepData = zeros(length(HbT_ProcAsleepData{1,1}),length(HbT_ProcAsleepData));
            Neural_asleepData = zeros(length(Neural_ProcAsleepData{1,1}),length(Neural_ProcAsleepData));
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
            for cc = 1:length(HbT_ProcAsleepData)
                HbT_asleepData(:,cc) = HbT_ProcAsleepData{cc,1};
                Neural_asleepData(:,cc) = Neural_ProcAsleepData{cc,1};
            end
            % calculate the coherence between desired signals
            params.tapers = [7,13]; % Tapers [n, 2n - 1]
            [C_AsleepData,~,~,~,~,f_AsleepData,confC_AsleepData,~,cErr_AsleepData] = coherencyc(HbT_asleepData,Neural_asleepData,params);
            % save results
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.C = C_AsleepData;
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.f = f_AsleepData;
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.confC = confC_AsleepData;
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.cErr = cErr_AsleepData;
        else
            % save results
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.C = [];
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.f = [];
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.confC = [];
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.cErr = [];
        end
        %% All
        zz = 1;
        clear HbT_AllData Neural_AllData HbT_ProcAllData Neural_ProcAllData
        HbT_AllData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allDataFileDate,~] = GetFileInfo_IOS_nNOS(procDataFileID);
            strDay = ConvertDate_IOS_nNOS(allDataFileDate);
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                motionArtifact = ProcData.notes.motionArtifact;
                if motionArtifact == false
                    HbT_AllData{zz,1} = ProcData.data.HbT.(hemisphere)(1:300*samplingRate);
                    Neural_AllData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType)(1:300*samplingRate) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                    zz = zz + 1;
                    HbT_AllData{zz,1} = ProcData.data.HbT.(hemisphere)(300*samplingRate + 1:600*samplingRate);
                    Neural_AllData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType)(300*samplingRate + 1:600*samplingRate) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                    zz = zz + 1;
                    HbT_AllData{zz,1} = ProcData.data.HbT.(hemisphere)(600*samplingRate + 1:end);
                    Neural_AllData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType)(600*samplingRate + 1:end) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                    zz = zz + 1;
                end
            end
        end
        % filter and detrend data
        if isempty(HbT_AllData) == false
            for cc = 1:length(HbT_AllData)
                HbT_ProcAllUnstimData{cc,1} = detrend(HbT_AllData{cc,1},'constant');
                Neural_ProcAllUnstimData{cc,1} = detrend(Neural_AllData{cc,1},'constant');
            end
            % pre-allocate coherence matrix
            HbT_allData = zeros(length(HbT_ProcAllUnstimData{1,1}),length(HbT_ProcAllUnstimData));
            Neural_allData = zeros(length(Neural_ProcAllUnstimData{1,1}),length(Neural_ProcAllUnstimData));
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
            for cc = 1:length(HbT_ProcAllUnstimData)
                HbT_allData(:,cc) = HbT_ProcAllUnstimData{cc,1};
                Neural_allData(:,cc) = Neural_ProcAllUnstimData{cc,1};
            end
            % calculate the coherence between desired signals
            params.tapers = [7,13]; % Tapers [n, 2n - 1]
            [C_AllUnstimData,~,~,~,~,f_AllUnstimData,confC_AllUnstimData,~,cErr_AllUnstimData] = coherencyc(HbT_allData,Neural_allData,params);
            % save results
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).All.C = C_AllUnstimData;
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).All.f = f_AllUnstimData;
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).All.confC = confC_AllUnstimData;
            Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).All.cErr = cErr_AllUnstimData;
        end
        %% NREM
        [HbT_nremData,~,~] = RemoveStimSleepData_IOS_nNOS(animalID,SleepData.(modelType).NREM.data.HbT.(hemisphere),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [Neural_nremData,~,~] = RemoveStimSleepData_IOS_nNOS(animalID,SleepData.(modelType).NREM.data.(['cortical_' hemisphere]).(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        % filter, detrend, and truncate data to minimum length to match events
        for ee = 1:length(HbT_nremData)
            HbT_nremData{ee,1} = detrend(HbT_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
            Neural_nremData{ee,1} = detrend(Neural_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        end
        % pre-allocate coherence matrix
        HbT_nrem = zeros(length(HbT_nremData{1,1}),length(HbT_nremData));
        Neural_nrem = zeros(length(Neural_nremData{1,1}),length(Neural_nremData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for ff = 1:length(HbT_nremData)
            HbT_nrem(:,ff) = HbT_nremData{ff,1};
            Neural_nrem(:,ff) = Neural_nremData{ff,1};
        end
        % calculate the coherence between desired signals
        params.tapers = [3,5]; % Tapers [n, 2n - 1]
        [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc(HbT_nrem,Neural_nrem,params);
        % save results
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).NREM.C = C_nrem;
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).NREM.f = f_nrem;
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).NREM.confC = confC_nrem;
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).NREM.cErr = cErr_nrem;
        %% REM
        [HbT_remData,~,~] = RemoveStimSleepData_IOS_nNOS(animalID,SleepData.(modelType).REM.data.HbT.(hemisphere),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [Neural_remData,~,~] = RemoveStimSleepData_IOS_nNOS(animalID,SleepData.(modelType).REM.data.(['cortical_' hemisphere]).(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        % filter, detrend, and truncate data to minimum length to match events
        for gg = 1:length(HbT_remData)
            HbT_remData{gg,1} = detrend(HbT_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant');
            Neural_remData{gg,1} = detrend(Neural_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant');
        end
        % pre-allocate coherence matrix
        HbT_rem = zeros(length(HbT_remData{1,1}),length(HbT_remData));
        Neural_rem = zeros(length(Neural_remData{1,1}),length(Neural_remData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for hh = 1:length(HbT_remData)
            HbT_rem(:,hh) = HbT_remData{hh,1};
            Neural_rem(:,hh) = Neural_remData{hh,1};
        end
        % calculate the coherence between desired signals
        params.tapers = [5,9]; % Tapers [n, 2n - 1]
        [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc(HbT_rem,Neural_rem,params);
        % save results
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).REM.C = C_rem;
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).REM.f = f_rem;
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).REM.confC = confC_rem;
        Results_NeuralHemoCoher_Ephys.(group).(animalID).(hemisphere).(dataType).REM.cErr = cErr_rem;
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_NeuralHemoCoher_Ephys.mat','Results_NeuralHemoCoher_Ephys')
cd([rootFolder delim 'Data'])