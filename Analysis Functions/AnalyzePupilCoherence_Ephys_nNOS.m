function [Results_PupilCoher_Ephys] = AnalyzePupilCoherence_Ephys(animalID,group,set,rootFolder,delim,Results_PupilCoher_Ephys)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% go to animal's data location
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% find and load manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% find and load RestingBaselines.mat struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData.mat struct
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% find and load ScoringResults.mat struct
scoringResultsFileStruct = dir('*Forest_ScoringResults.mat');
scoringResultsFile = {scoringResultsFileStruct.name}';
scoringResultsFileID = char(scoringResultsFile);
load(scoringResultsFileID,'-mat')
% lowpass filter
samplingRate = RestData.CBV.LH.samplingRate;
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
hemDataTypes = {'LH_HbT','RH_HbT','LH_gammaBandPower','RH_gammaBandPower'};
% go through each valid data type for arousal-based coherence analysis
for zzz = 1:length(hemDataTypes)
    hemDataType = hemDataTypes{1,zzz};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        %% analyze neural-hemo coherence during periods of rest
        % pull data from RestData.mat structure
        [restLogical] = FilterEvents_IOS(RestData.Pupil.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.Pupil.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.Pupil.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.Pupil.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.Pupil.(dataType).durations(combRestLogical,:);
        HbT_restData = RestData.Pupil.(hemDataType).data(combRestLogical,:);
        Pupil_restData = RestData.Pupil.(dataType).data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [HbT_finalRestData,~,~,~] = RemoveInvalidData_IOS(HbT_restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [Pupil_FinalRestData,~,~,~] = RemoveInvalidData_IOS(Pupil_restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        clear HbT_procRestData Pupil_procRestData
        % filter, detrend, and truncate data to minimum length to match events
        for bb = 1:length(HbT_finalRestData)
            if length(HbT_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(HbT_finalRestData{bb,1});
                HbT_restPad = (ones(1,restChunkSampleDiff))*HbT_finalRestData{bb,1}(end);
                Pupil_restPad = (ones(1,restChunkSampleDiff))*Pupil_FinalRestData{bb,1}(end);
                HbT_procRestData{bb,1} = horzcat(HbT_finalRestData{bb,1},HbT_restPad);
                Pupil_procRestData{bb,1} = horzcat(Pupil_FinalRestData{bb,1},Pupil_restPad);
                HbT_procRestData{bb,1} = detrend(HbT_procRestData{bb,1},'constant');
                Pupil_procRestData{bb,1} = detrend(Pupil_procRestData{bb,1},'constant');
            else
                HbT_procRestData{bb,1} = detrend(HbT_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
                Pupil_procRestData{bb,1} = detrend(Pupil_FinalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            end
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        zz = 1;
        for cc = 1:length(HbT_procRestData)
            if sum(isnan(Pupil_procRestData{cc,1})) == 0
                HbT_restDataMat(:,zz) = HbT_procRestData{cc,1};
                Pupil_restDataMat(:,zz) = Pupil_procRestData{cc,1};
                zz = zz + 1;
            end
        end
        % parameters for coherencyc - information available in function
        params.tapers = [1,1]; % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1]; % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc(HbT_restDataMat,Pupil_restDataMat,params);
        % save results
        Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Rest.C = C_RestData;
        Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Rest.f = f_RestData;
        Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Rest.confC = confC_RestData;
        Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Rest.cErr = cErr_RestData;
        %% analyze neural-hemo coherence during periods of alert
        zz = 1;
        clear HbT_alertData Pupil_alertData HbT_procAlertData Pupil_procAlertData
        HbT_alertData = []; Pupil_alertData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,fileDate,alertDataFileID] = GetFileInfo_IOS(procDataFileID);
            strDay = ConvertDate_IOS(fileDate);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(alertDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % 36 bins (180 total) or 3 minutes of asleep
                load(procDataFileID,'-mat')
                % don't include trials with stimulation
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                    if isempty(puffs) == true
                        if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                            if strcmp(hemDataType,'LH_HbT') == true
                                HbT_alertData{zz,1} = ProcData.data.HbT.LH;
                            elseif strcmp(hemDataType,'RH_HbT') == true
                                HbT_alertData{zz,1} = ProcData.data.HbT.RH;
                            elseif strcmp(hemDataType,'LH_gammaBandPower') == true
                                HbT_alertData{zz,1} = (ProcData.data.cortical_LH.gammaBandPower - RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay).mean;
                            elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                                HbT_alertData{zz,1} = (ProcData.data.cortical_RH.gammaBandPower - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                            end
                            Pupil_alertData{zz,1} = ProcData.data.Pupil.(dataType);
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(HbT_alertData) == false
            for bb = 1:length(HbT_alertData)
                HbT_procAlertData{bb,1} = detrend(HbT_alertData{bb,1},'constant');
                Pupil_procAlertData{bb,1} = detrend(Pupil_alertData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            HbT_alertDataMat = zeros(length(HbT_procAlertData{1,1}),length(HbT_procAlertData));
            Pupil_alertDataMat = zeros(length(Pupil_procAlertData{1,1}),length(Pupil_procAlertData));
            for cc = 1:length(HbT_procAlertData)
                HbT_alertDataMat(:,cc) = HbT_procAlertData{cc,1};
                Pupil_alertDataMat(:,cc) = Pupil_procAlertData{cc,1};
            end
            % calculate the coherence between desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            [C_AlertData,~,~,~,~,f_AlertData,confC_AlertData,~,cErr_AlertData] = coherencyc(HbT_alertDataMat,Pupil_alertDataMat,params);
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Alert.C = C_AlertData;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Alert.f = f_AlertData;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Alert.confC = confC_AlertData;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Alert.cErr = cErr_AlertData;
        else
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Alert.C = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Alert.f = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Alert.confC = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Alert.cErr = [];
        end
        % analyze neural-hemo coherence during periods of aasleep
        zz = 1;
        clear HbT_asleepData Pupil_asleepData HbT_procAsleepData Pupil_procAsleepData
        HbT_asleepData = []; Pupil_asleepData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,fileDate,asleepDataFileID] = GetFileInfo_IOS(procDataFileID);
            strDay = ConvertDate_IOS(fileDate);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(asleepDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) < 36 % 36 bins (180 total) or 3 minutes of alert
                load(procDataFileID,'-mat')
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                    % don't include trials with stimulation
                    if isempty(puffs) == true
                        if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                            if strcmp(hemDataType,'LH_HbT') == true
                                HbT_asleepData{zz,1} = ProcData.data.HbT.LH;
                            elseif strcmp(hemDataType,'RH_HbT') == true
                                HbT_asleepData{zz,1} = ProcData.data.HbT.RH;
                            elseif strcmp(hemDataType,'LH_gammaBandPower') == true
                                HbT_asleepData{zz,1} = (ProcData.data.cortical_LH.gammaBandPower - RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay).mean;
                            elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                                HbT_asleepData{zz,1} = (ProcData.data.cortical_RH.gammaBandPower - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                            end
                            Pupil_asleepData{zz,1} = ProcData.data.Pupil.(dataType);
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(HbT_asleepData) == false
            for bb = 1:length(HbT_asleepData)
                HbT_procAsleepData{bb,1} = detrend(HbT_asleepData{bb,1},'constant');
                Pupil_procAsleepData{bb,1} = detrend(Pupil_asleepData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            HbT_asleepDataMat = zeros(length(HbT_procAsleepData{1,1}),length(HbT_procAsleepData));
            Pupil_asleepDataMat = zeros(length(Pupil_procAsleepData{1,1}),length(Pupil_procAsleepData));
            for cc = 1:length(HbT_procAsleepData)
                HbT_asleepDataMat(:,cc) = HbT_procAsleepData{cc,1};
                Pupil_asleepDataMat(:,cc) = Pupil_procAsleepData{cc,1};
            end
            % calculate the coherence between desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            [C_AsleepData,~,~,~,~,f_AsleepData,confC_AsleepData,~,cErr_AsleepData] = coherencyc(HbT_asleepDataMat,Pupil_asleepDataMat,params);
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Asleep.C = C_AsleepData;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Asleep.f = f_AsleepData;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Asleep.confC = confC_AsleepData;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Asleep.cErr = cErr_AsleepData;
        else
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Asleep.C = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Asleep.f = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Asleep.confC = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).Asleep.cErr = [];
        end
        % analyze neural-hemo coherence during periods of all data
        zz = 1;
        clear HbT_allData Pupil_allData HbT_allData Pupil_allData
        HbT_allData = []; Pupil_allData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
            strDay = ConvertDate_IOS(fileDate);
            load(procDataFileID,'-mat')
            if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                try
                    puffs = ProcData.data.stimulations.LPadSol;
                catch
                    puffs = ProcData.data.solenoids.LPadSol;
                end
                % don't include trials with stimulation
                if isempty(puffs) == true
                    if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                        if strcmp(hemDataType,'LH_HbT') == true
                            HbT_allData{zz,1} = ProcData.data.HbT.LH;
                        elseif strcmp(hemDataType,'RH_HbT') == true
                            HbT_allData{zz,1} = ProcData.data.HbT.RH;
                        elseif strcmp(hemDataType,'LH_gammaBandPower') == true
                            HbT_allData{zz,1} = (ProcData.data.cortical_LH.gammaBandPower - RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay).mean;
                        elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                            HbT_allData{zz,1} = (ProcData.data.cortical_RH.gammaBandPower - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                        end
                        Pupil_allData{zz,1} = ProcData.data.Pupil.(dataType);
                        zz = zz + 1;
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(HbT_allData) == false
            for bb = 1:length(HbT_allData)
                HbT_procAllData{bb,1} = detrend(HbT_allData{bb,1},'constant');
                Pupil_procAllData{bb,1} = detrend(Pupil_allData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            HbT_allDataMat = zeros(length(HbT_procAllData{1,1}),length(HbT_procAllData));
            Pupil_allDataMat = zeros(length(Pupil_procAllData{1,1}),length(Pupil_procAllData));
            for cc = 1:length(HbT_procAllData)
                HbT_allDataMat(:,cc) = HbT_procAllData{cc,1};
                Pupil_allDataMat(:,cc) = Pupil_procAllData{cc,1};
            end
            % calculate the coherence between desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            [C_AllData,~,~,~,~,f_AllData,confC_AllData,~,cErr_AllData] = coherencyc(HbT_allDataMat,Pupil_allDataMat,params);
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).All.C = C_AllData;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).All.f = f_AllData;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).All.confC = confC_AllData;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).All.cErr = cErr_AllData;
        else
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).All.C = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).All.f = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).All.confC = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).All.cErr = [];
        end
        % analyze neural-hemo coherence during periods of NREM
        % pull data from AsleepData.mat structure
        if isempty(SleepData.(modelType).NREM.data.Pupil) == false
            [HbT_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.(hemDataType).data,SleepData.(modelType).NREM.data.Pupil.fileIDs,SleepData.(modelType).NREM.data.Pupil.binTimes);
            [Pupil_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.(dataType).data,SleepData.(modelType).NREM.data.Pupil.fileIDs,SleepData.(modelType).NREM.data.Pupil.binTimes);
        else
            HbT_nremData = [];
            Pupil_nremData = [];
        end
        if isempty(HbT_nremData) == false
            clear HbT_procNremData Pupil_procNremData
            % filter, detrend, and truncate data to minimum length to match events
            for ee = 1:length(HbT_nremData)
                HbT_procNremData{ee,1} = detrend(HbT_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
                Pupil_procNremData{ee,1} = detrend(Pupil_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            zz = 1;
            for ff = 1:length(HbT_procNremData)
                if sum(isnan(Pupil_procNremData{ff,1})) == 0
                    HbT_nremMat(:,zz) = HbT_procNremData{ff,1};
                    Pupil_nremMat(:,zz) = Pupil_procNremData{ff,1};
                    zz = zz + 1;
                end
            end
            % calculate the coherence between desired signals
            params.tapers = [3,5]; % Tapers [n, 2n - 1]
            [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc(HbT_nremMat,Pupil_nremMat,params);
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).NREM.C = C_nrem;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).NREM.f = f_nrem;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).NREM.confC = confC_nrem;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).NREM.cErr = cErr_nrem;
        else
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).NREM.C = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).NREM.f = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).NREM.confC = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).NREM.cErr = [];
        end
        % analyze neural-hemo coherence during periods of REM
        % pull data from AsleepData.mat structure
        if isempty(SleepData.(modelType).REM.data.Pupil) == false
            [HbT_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.(hemDataType).data,SleepData.(modelType).REM.data.Pupil.fileIDs,SleepData.(modelType).REM.data.Pupil.binTimes);
            [Pupil_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.(dataType).data,SleepData.(modelType).REM.data.Pupil.fileIDs,SleepData.(modelType).REM.data.Pupil.binTimes);
        else
            HbT_remData = [];
            Pupil_remData = [];
        end
        if isempty(HbT_remData) == false
            clear HbT_procRemData Pupil_procRemData
            % filter, detrend, and truncate data to minimum length to match events
            for ee = 1:length(HbT_remData)
                HbT_procRemData{ee,1} = detrend(HbT_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
                Pupil_procRemData{ee,1} = detrend(Pupil_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            zz = 1;
            for ff = 1:length(HbT_procRemData)
                if sum(isnan(Pupil_procRemData{ff,1})) == 0
                    HbT_remMat(:,zz) = HbT_procRemData{ff,1};
                    Pupil_remMat(:,zz) = Pupil_procRemData{ff,1};
                    zz = zz + 1;
                end
            end
            % calculate the coherence between desired signals
            params.tapers = [5,9]; % Tapers [n, 2n - 1]
            [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc(HbT_remMat,Pupil_remMat,params);
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).REM.C = C_rem;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).REM.f = f_rem;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).REM.confC = confC_rem;
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).REM.cErr = cErr_rem;
        else
            % save results
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).REM.C = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).REM.f = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).REM.confC = [];
            Results_PupilCoher_Ephys.(group).(animalID).(dataType).(hemDataType).REM.cErr = [];
        end
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PupilCoher_Ephys.mat','Results_PupilCoher_Ephys')
cd([rootFolder delim 'Data'])