function [Results_PreWhitenedPowerSpec_Ephys] = AnalyzePreWhitenedPowerSpectrum_Ephys_eLife2025(animalID,group,set,rootFolder,delim,Results_PreWhitenedPowerSpec_Ephys)
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
hemispheres = {'LH','RH'};
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        %% Rest
        clear restingData procRestData restData restData
        if strcmp(dataType,'HbT') == true
            samplingRate = RestData.(dataType).(hemisphere).samplingRate;
            [restLogical] = FilterEvents_IOS_eLife2025(RestData.(dataType).(hemisphere),RestCriteria);
            [stimLogical] = FilterEvents_IOS_eLife2025(RestData.(dataType).(hemisphere),RestStimCriteria);
            combRestLogical = logical(restLogical.*stimLogical);
            restFileIDs = RestData.(dataType).(hemisphere).fileIDs(combRestLogical,:);
            restEventTimes = RestData.(dataType).(hemisphere).eventTimes(combRestLogical,:);
            restDurations = RestData.(dataType).(hemisphere).durations(combRestLogical,:);
            restingData = RestData.(dataType).(hemisphere).data(combRestLogical,:);
            % keep only the data that occurs within the manually-approved alert regions
            [restData,~,~,~] = RemoveInvalidData_IOS_eLife2025(restingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        else
            samplingRate = RestData.(['cortical_' hemisphere]).(dataType).samplingRate;
            [restLogical] = FilterEvents_IOS_eLife2025(RestData.(['cortical_' hemisphere]).(dataType),RestCriteria);
            [stimLogical] = FilterEvents_IOS_eLife2025(RestData.(['cortical_' hemisphere]).(dataType),RestStimCriteria);
            combRestLogical = logical(restLogical.*stimLogical);
            restFileIDs = RestData.(['cortical_' hemisphere]).(dataType).fileIDs(combRestLogical,:);
            restEventTimes = RestData.(['cortical_' hemisphere]).(dataType).eventTimes(combRestLogical,:);
            restDurations = RestData.(['cortical_' hemisphere]).(dataType).durations(combRestLogical,:);
            restingData = RestData.(['cortical_' hemisphere]).(dataType).NormData(combRestLogical,:);
            % keep only the data that occurs within the manually-approved awake regions
            [restData,~,~,~] = RemoveInvalidData_IOS_eLife2025(restingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        end
        % detrend and truncate data to minimum length to match events
        for cc = 1:length(restData)
            if length(restData{cc,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(restData{cc,1});
                restPad = (ones(1,restChunkSampleDiff))*restData{cc,1}(end);
                procRestData{cc,1} = horzcat(restData{cc,1},restPad);
                procRestData{cc,1} = diff(detrend(procRestData{cc,1},'constant'),1);
                finalRestData(:,cc) = procRestData{cc,1};
            else
                procRestData{cc,1} = diff(detrend(restData{cc,1}(1:(params.minTime.Rest*samplingRate)),'constant'),1);
                finalRestData(:,cc) = procRestData{cc,1};
            end
        end
        % parameters for mtspectrumc - information available in function
        params.tapers = [1,1]; % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1]; % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        [restS,restf,restsErr] = mtspectrumc(finalRestData,params);
        % save results
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Rest.S = restS;
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Rest.f = restf;
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Rest.sErr = restsErr;
        %% Alert
        clear alertData procAlertData finalAlertData scoringLabels
        alertData = []; % for loop pre-allocation
        zz = 1;
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,fileDate,fileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
            strDay = ConvertDate_IOS_eLife2025(fileDate);
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(fileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % less than 3 minutes of asleep
                load(procDataFileID,'-mat')
                stims = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(stims) == true
                    if strcmp(dataType,'HbT') == true
                        alertData{zz,1} = ProcData.data.(dataType).(hemisphere);
                        zz = zz + 1;
                    else
                        motionArtifact = ProcData.notes.motionArtifact;
                        if motionArtifact == false
                            alertData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % detrend data
        if isempty(alertData) == false
            for cc = 1:length(alertData)
                procAlertData{cc,1} = diff(detrend(alertData{cc,1},'constant'),1);
                finalAlertData(:,cc) = procAlertData{cc,1};
            end
            % calculate the power spectra of the desired signals
            params.tapers = [5,9]; % Tapers [n, 2n - 1]
            [alertS,alertf,alertsErr] = mtspectrumc(finalAlertData,params);
            % save results
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.S = alertS;
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.f = alertf;
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.sErr = alertsErr;
        else
            % save results
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.S = [];
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.f = [];
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.sErr = [];
        end
        %% Asleep
        clear asleepData procAsleepData finalAsleepData scoringLabels
        asleepData = []; % for loop pre-allocation
        zz = 1;
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,fileDate,fileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
            strDay = ConvertDate_IOS_eLife2025(fileDate);
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(fileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) < 36 % less than 3 minutes of alert
                load(procDataFileID,'-mat')
                stims = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(stims) == true
                    if strcmp(dataType,'HbT') == true
                        asleepData{zz,1} = ProcData.data.(dataType).(hemisphere);
                        zz = zz + 1;
                    else
                        motionArtifact = ProcData.notes.motionArtifact;
                        if motionArtifact == false
                            asleepData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % detrend data
        if isempty(asleepData) == false
            for cc = 1:length(asleepData)
                procAsleepData{cc,1} = diff(detrend(asleepData{cc,1},'constant'),1);
                finalAsleepData(:,cc) = procAsleepData{cc,1};
            end
            % calculate the power spectra of the desired signals
            params.tapers = [5,9]; % Tapers [n, 2n - 1]
            [asleepS,asleepf,asleepsErr] = mtspectrumc(finalAsleepData,params);
            % save results
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.S = asleepS;
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.f = asleepf;
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.sErr = asleepsErr;
        else
            % save results
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.S = [];
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.f = [];
            Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.sErr = [];
        end
        %% All
        clear allData procAllData finalAllData
        zz = 1;
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,fileDate,~] = GetFileInfo_IOS_eLife2025(procDataFileID);
            strDay = ConvertDate_IOS_eLife2025(fileDate);
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                if strcmp(dataType,'HbT') == true
                    allData{zz,1} = ProcData.data.(dataType).(hemisphere);
                    zz = zz + 1;
                else
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        allData{zz,1} = (ProcData.data.(['cortical_' hemisphere]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' hemisphere]).(dataType).(strDay).mean;
                        zz = zz + 1;
                    end
                end
            end
        end
        % detrend data
        for cc = 1:length(allData)
            procAllData{cc,1} = diff(detrend(allData{cc,1},'constant'),1);
            finalAllData(:,cc) = procAllData{cc,1};
        end
        % calculate the power spectra of the desired signals
        params.tapers = [5,9]; % Tapers [n, 2n - 1]
        [allS,allf,allsErr] = mtspectrumc(finalAllData,params);
        % save results
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).All.S = allS;
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).All.f = allf;
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).All.sErr = allsErr;
        %% NREM
        clear nremData procNREMData finalNREMData
        if strcmp(dataType,'HbT') == true
            [nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(dataType).(hemisphere),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        else
            [nremData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).NREM.data.(['cortical_' hemisphere]).(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        end
        % detrend and truncate data to minimum length to match events
        for ee = 1:length(nremData)
            procNREMData{ee,1} = diff(detrend(nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant'),1);
            finalNREM(:,ee) = procNREMData{ee,1};
        end
        % calculate the power spectra of the desired signals
        params.tapers = [3,5]; % Tapers [n, 2n - 1]
        [nremS,nremf,nremsErr] = mtspectrumc(finalNREM,params);
        % save results
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).NREM.S = nremS;
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).NREM.f = nremf;
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).NREM.sErr = nremsErr;
        %% REM
        clear remData procREMData finalREMData
        if strcmp(dataType,'HbT') == true
            [remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(dataType).(hemisphere),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        else
            [remData,~,~] = RemoveStimSleepData_IOS_eLife2025(animalID,SleepData.(modelType).REM.data.(['cortical_' hemisphere]).(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        end
        % detrend and truncate data to minimum length to match events
        for ee = 1:length(remData)
            procREMData{ee,1} = diff(detrend(remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant'),1);
            finalREM(:,ee) = procREMData{ee,1};
        end
        % calculate the power spectra of the desired signals
        params.tapers = [3,5]; % Tapers [n, 2n - 1]
        [remS,remf,remsErr] = mtspectrumc(finalREM,params);
        % save results
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).REM.S = remS;
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).REM.f = remf;
        Results_PreWhitenedPowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).REM.sErr = remsErr;
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PreWhitenedPowerSpec_Ephys.mat','Results_PreWhitenedPowerSpec_Ephys')
cd([rootFolder delim 'Data'])