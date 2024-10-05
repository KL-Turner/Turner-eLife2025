function [Results_PupilCrossCorr_Ephys] = AnalyzePupilCrossCorrelation_Ephys_nNOS(animalID,group,set,rootFolder,delim,Results_PupilCrossCorr_Ephys)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% identify and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% identify and load ManualDecisions struct
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% identify and load RestingBaselines struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% identify and load SleepData struct
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% identify and load ScoringResults struct
scoringResultsFileStruct = dir('*Forest_ScoringResults.mat');
scoringResultsFile = {scoringResultsFileStruct.name}';
scoringResultsFileID = char(scoringResultsFile);
load(scoringResultsFileID,'-mat')
% parameters
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% filter characteristics & resting criteria
samplingRate = RestData.HbT.LH.samplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
dataTypes = {'LH_HbT','RH_HbT','LH_gammaBandPower','RH_gammaBandPower'};
pupilDataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
% go through each valid data type for arousal-dependent cross-correlation analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(pupilDataTypes)
        pupilDataType = pupilDataTypes{1,bb};
        %% cross-correlation analysis for resting data
        % pull data from RestData.mat structure
        [restLogical] = FilterEvents_IOS_nNOS(RestData.Pupil.(pupilDataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS_nNOS(RestData.Pupil.(pupilDataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.Pupil.(pupilDataType).fileIDs(combRestLogical,:);
        restDurations = RestData.Pupil.(pupilDataType).durations(combRestLogical,:);
        restEventTimes = RestData.Pupil.(pupilDataType).eventTimes(combRestLogical,:);
        restData = RestData.Pupil.(dataType).data(combRestLogical,:);
        pupilRestData = RestData.Pupil.(pupilDataType).data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [finalRestData,~,~,~] = RemoveInvalidData_IOS_nNOS(restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [finalPupilRestData,~,~,~] = RemoveInvalidData_IOS_nNOS(pupilRestData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        % process, filter + detrend each array
        catRestData = [];
        cc = 1;
        for dd = 1:length(finalRestData)
            if sum(isnan(finalPupilRestData{dd,1})) == 0
                if length(finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                    restChunkSampleDiff = params.minTime.Rest*samplingRate - length(finalRestData{bb,1});
                    restPad = (ones(1,restChunkSampleDiff))*finalRestData{bb,1}(end);
                    restPupilPad = (ones(1,restChunkSampleDiff))*finalPupilRestData{bb,1}(end);
                    procRestData = horzcat(finalRestData{bb,1},restPad);
                    procPupilRestData = horzcat(finalPupilRestData{bb,1},restPupilPad);
                    procRestData = filtfilt(sos,g,detrend(procRestData{bb,1},'constant'));
                    procPupilRestData = filtfilt(sos,g,detrend(procPupilRestData{bb,1},'constant'));
                else
                    procRestData = filtfilt(sos,g,detrend(finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                    procPupilRestData = filtfilt(sos,g,detrend(finalPupilRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                end
                catRestData.data{cc,1} = procRestData;
                catRestData.pupil{cc,1} = procPupilRestData;
                cc = cc + 1;
            end
        end
        % run cross correlation between data types
        restXcVals = [];
        if isempty(catRestData) == false
            restLagTime = 5; % seconds
            restFrequency = samplingRate; % Hz
            restMaxLag = restLagTime*restFrequency;
            % run cross-correlation analysis - average through time
            for ee = 1:length(catRestData.data)
                restArray = catRestData.data{ee,1};
                restPupilarray = catRestData.pupil{ee,1};
                [restXcVals(ee,:),restPupilLags] = xcorr(restArray,restPupilarray,restMaxLag,'coeff');
            end
            restMeanXcVals = mean(restXcVals,1);
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Rest.lags = restPupilLags;
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Rest.xcVals = restMeanXcVals;
        else
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Rest.lags = [];
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Rest.xcVals = [];
        end
        %% cross-correlation analysis for alert data
        zz = 1;
        alertData = []; alertPupilData = []; alertProcData = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,alertDataFileID] = GetFileInfo_IOS_nNOS(procDataFileID);
            scoringLabels = [];
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(alertDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % 12 minutes of awake
                load(procDataFileID,'-mat')
                % only run on files with good pupil measurement
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                    % don't include trials with stimulation
                    if isempty(puffs) == true
                        if sum(isnan(ProcData.data.Pupil.(pupilDataType))) == 0
                            % pull data based on data type
                            if strcmp(dataType,'LH_HbT') == true
                                alertData{zz,1} = ProcData.data.HbT.LH;
                            elseif strcmp(dataType,'RH_HbT') == true
                                alertData{zz,1} = ProcData.data.HbT.RH;
                            elseif strcmp(dataType,'LH_gammaBandPower') == true
                                alertData{zz,1} = ProcData.data.cortical_LH.gammaBandPower;
                            elseif strcmp(dataType,'RH_gammaBandPower') == true
                                alertData{zz,1} = ProcData.data.cortical_RH.gammaBandPower;
                            end
                            alertPupilData{zz,1} = ProcData.data.Pupil.(pupilDataType);
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % process, filter + detrend each array
        alertXcVals = [];
        if isempty(alertData) == false
            for ee = 1:length(alertData)
                alertProcData.data{ee,1} = filtfilt(sos,g,detrend(alertData{ee,1},'constant'));
                alertProcData.pupil{ee,1} = filtfilt(sos,g,detrend(alertPupilData{ee,1},'constant'));
            end
            % set parameters for cross-correlation analysis
            alertLagTime = 30; % seconds
            alertFrequency = samplingRate; % Hz
            alertMaxLag = alertLagTime*alertFrequency;
            % run cross-correlation analysis - average through time
            for ff = 1:length(alertProcData.data)
                alertArray = alertProcData.data{ff,1};
                alertPupilarray = alertProcData.pupil{ff,1};
                [alertXcVals(ff,:),alertPupilLags] = xcorr(alertArray,alertPupilarray,alertMaxLag,'coeff');
            end
            alertMeanXcVals = mean(alertXcVals,1);
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Alert.lags = alertPupilLags;
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Alert.xcVals = alertMeanXcVals;
        else
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Alert.lags = [];
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Alert.xcVals = [];
        end
        %% cross-correlation analysis for asleep data
        zz = 1;
        asleepData = []; asleepPupilData = []; asleepProcData = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,asleepDataFileID] = GetFileInfo_IOS_nNOS(procDataFileID);
            scoringLabels = [];
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(asleepDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) < 36 % 12 minutes of asleep
                load(procDataFileID,'-mat')
                % only run on files with good pupil measurement
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                    % don't include trials with stimulation
                    if isempty(puffs) == true
                        if sum(isnan(ProcData.data.Pupil.(pupilDataType))) == 0
                            if strcmp(dataType,'LH_HbT') == true
                                asleepData{zz,1} = ProcData.data.HbT.LH;
                            elseif strcmp(dataType,'RH_HbT') == true
                                asleepData{zz,1} = ProcData.data.HbT.RH;
                            elseif strcmp(dataType,'LH_gammaBandPower') == true
                                asleepData{zz,1} = ProcData.data.cortical_LH.gammaBandPower;
                            elseif strcmp(dataType,'RH_gammaBandPower') == true
                                asleepData{zz,1} = ProcData.data.cortical_RH.gammaBandPower;
                            end
                            asleepPupilData{zz,1} = ProcData.data.Pupil.(pupilDataType);
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % process, filter + detrend each array
        asleepXcVals = [];
        if isempty(asleepData) == false
            for ee = 1:length(asleepData)
                asleepProcData.data{ee,1} = filtfilt(sos,g,detrend(asleepData{ee,1},'constant'));
                asleepProcData.pupil{ee,1} = filtfilt(sos,g,detrend(asleepPupilData{ee,1},'constant'));
            end
            % set parameters for cross-correlation analysis
            asleepLagTime = 30; % seconds
            asleepFrequency = samplingRate; % Hz
            asleepMaxLag = asleepLagTime*asleepFrequency;
            % run cross-correlation analysis - average through time
            for ff = 1:length(asleepProcData.data)
                asleepArray = asleepProcData.data{ff,1};
                asleepPupilarray = asleepProcData.pupil{ff,1};
                [asleepXcVals(ff,:),asleepPupilLags] = xcorr(asleepArray,asleepPupilarray,asleepMaxLag,'coeff');
            end
            asleepMeanXcVals = mean(asleepXcVals,1);
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Asleep.lags = asleepPupilLags;
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Asleep.xcVals = asleepMeanXcVals;
        else
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Asleep.lags = [];
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).Asleep.xcVals = [];
        end
        %% cross-correlation analysis for all data
        zz = 1;
        allData = []; allPupilData = []; allProcData = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,~] = GetFileInfo_IOS_nNOS(procDataFileID);
            load(procDataFileID,'-mat')
            % only run on files with good pupil measurement
            if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                try
                    puffs = ProcData.data.stimulations.LPadSol;
                catch
                    puffs = ProcData.data.solenoids.LPadSol;
                end
                % don't include trials with stimulation
                if isempty(puffs) == true
                    if sum(isnan(ProcData.data.Pupil.(pupilDataType))) == 0
                        if strcmp(dataType,'LH_HbT') == true
                            allData{zz,1} = ProcData.data.HbT.LH;
                        elseif strcmp(dataType,'RH_HbT') == true
                            allData{zz,1} = ProcData.data.HbT.RH;
                        elseif strcmp(dataType,'LH_gammaBandPower') == true
                            allData{zz,1} = ProcData.data.cortical_LH.gammaBandPower;
                        elseif strcmp(dataType,'RH_gammaBandPower') == true
                            allData{zz,1} = ProcData.data.cortical_RH.gammaBandPower;
                        end
                        allPupilData{zz,1} = ProcData.data.Pupil.(pupilDataType);
                        if sum(isinf(ProcData.data.Pupil.(pupilDataType))) > 1
                            keyboard
                        end
                        zz = zz + 1;
                    end
                end
            end
        end
        % process, filter + detrend each array
        allXcVals = [];
        if isempty(allData) == false
            for dd = 1:length(allData)
                allProcData.data{dd,1} = filtfilt(sos,g,detrend(allData{dd,1},'constant'));
                allProcData.pupil{dd,1} = filtfilt(sos,g,detrend(allPupilData{dd,1},'constant'));
            end
            % set parameters for cross-correlation analysis
            allLagTime = 30; % seconds
            allFrequency = samplingRate; % Hz
            allMaxLag = allLagTime*allFrequency;
            % run cross-correlation analysis - average through time
            for ee = 1:length(allProcData.data)
                allHbTarray = allProcData.data{ee,1};
                allPupilarray = allProcData.pupil{ee,1};
                [allXcVals(ee,:),allPupilLags] = xcorr(allHbTarray,allPupilarray,allMaxLag,'coeff');
            end
            allMeanXcVals = mean(allXcVals,1);
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).All.lags = allPupilLags;
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).All.xcVals = allMeanXcVals;
        else
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).All.lags = [];
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).All.xcVals = [];
        end
        %% NREM
        if isempty(SleepData.(modelType).NREM.data.Pupil) == false
            NREM_sleepTime = params.minTime.NREM; % seconds
            [NREM_finalData,~,~] = RemoveStimSleepData_IOS_nNOS(animalID,SleepData.(modelType).NREM.data.Pupil.(dataType).data,SleepData.(modelType).NREM.data.Pupil.fileIDs,SleepData.(modelType).NREM.data.Pupil.binTimes);
            [NREM_finalPupilData,~,~] = RemoveStimSleepData_IOS_nNOS(animalID,SleepData.(modelType).NREM.data.Pupil.(pupilDataType).data,SleepData.(modelType).NREM.data.Pupil.fileIDs,SleepData.(modelType).NREM.data.Pupil.binTimes);
            NREM_finalVals = []; NREM_finalPupilVals = [];
            if isempty(NREM_finalData) == false
                % ust events to match the edits made to the length of each spectrogram
                dd = 1;
                for cc = 1:length(NREM_finalData)
                    if sum(isnan(NREM_finalPupilData{cc,1})) == 0
                        NREM_vals = NREM_finalData{cc,1}(1:NREM_sleepTime*samplingRate);
                        NREM_pupilVals = NREM_finalPupilData{cc,1}(1:NREM_sleepTime*samplingRate);
                        NREM_finalVals{dd,1} = filtfilt(sos,g,detrend(NREM_vals,'constant'));
                        NREM_finalPupilVals{dd,1} = filtfilt(sos,g,detrend(NREM_pupilVals,'constant'));
                        dd = dd + 1;
                    end
                end
                % process, filter + detrend each array
                NREM_xcVals = [];
                if isempty(NREM_finalVals) == false
                    % run cross-correlation analysis - average through time
                    NREM_lagTime = 15; % seconds
                    NREM_frequency = samplingRate; % Hz
                    NREM_maxLag = NREM_lagTime*NREM_frequency;
                    for ee = 1:length(NREM_finalVals)
                        NREM_array = NREM_finalVals{ee,1};
                        NREM_pupilArray = NREM_finalPupilVals{ee,1};
                        [NREM_xcVals(ee,:),NREM_PupilLags] = xcorr(NREM_array,NREM_pupilArray,NREM_maxLag,'coeff');
                    end
                    NREM_meanXcVals = mean(NREM_xcVals,1);
                    % save results
                    Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).NREM.lags = NREM_PupilLags;
                    Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).NREM.xcVals = NREM_meanXcVals;
                else
                    % save results
                    Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).NREM.lags = [];
                    Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).NREM.xcVals = [];
                end
            end
        else
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).NREM.lags = [];
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).NREM.xcVals = [];
        end
        %% REM
        if isempty(SleepData.(modelType).REM.data.Pupil) == false
            REM_sleepTime = params.minTime.REM; % seconds
            [REM_finalData,~,~] = RemoveStimSleepData_IOS_nNOS(animalID,SleepData.(modelType).REM.data.Pupil.(dataType).data,SleepData.(modelType).REM.data.Pupil.fileIDs,SleepData.(modelType).REM.data.Pupil.binTimes);
            [REM_finalPupilData,~,~] = RemoveStimSleepData_IOS_nNOS(animalID,SleepData.(modelType).REM.data.Pupil.(pupilDataType).data,SleepData.(modelType).REM.data.Pupil.fileIDs,SleepData.(modelType).REM.data.Pupil.binTimes);
            REM_finalVals = []; REM_finalPupilVals = [];
            if isempty(REM_finalData) == false
                % ust events to match the edits made to the length of each spectrogram
                dd = 1;
                for cc = 1:length(REM_finalData)
                    if sum(isnan(REM_finalPupilData{cc,1})) == 0
                        REM_vals = REM_finalData{cc,1}(1:REM_sleepTime*samplingRate);
                        REM_pupilVals = REM_finalPupilData{cc,1}(1:REM_sleepTime*samplingRate);
                        REM_finalVals{dd,1} = filtfilt(sos,g,detrend(REM_vals,'constant'));
                        REM_finalPupilVals{dd,1} = filtfilt(sos,g,detrend(REM_pupilVals,'constant'));
                        dd = dd + 1;
                    end
                end
                % process, filter + detrend each array
                REM_xcVals = [];
                if isempty(REM_finalVals) == false
                    % run cross-correlation analysis - average through time
                    REM_lagTime = 15; % seconds
                    REM_frequency = samplingRate; % Hz
                    REM_maxLag = REM_lagTime*REM_frequency;
                    for ee = 1:length(REM_finalVals)
                        REM_array = REM_finalVals{ee,1};
                        REM_pupilArray = REM_finalPupilVals{ee,1};
                        [REM_xcVals(ee,:),REM_PupilLags] = xcorr(REM_array,REM_pupilArray,REM_maxLag,'coeff');
                    end
                    REM_meanXcVals = mean(REM_xcVals,1);
                    % save results
                    Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).REM.lags = REM_PupilLags;
                    Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).REM.xcVals = REM_meanXcVals;
                else
                    % save results
                    Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).REM.lags = [];
                    Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).REM.xcVals = [];
                end
            end
        else
            % save results
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).REM.lags = [];
            Results_PupilCrossCorr_Ephys.(group).(animalID).(dataType).(pupilDataType).REM.xcVals = [];
        end
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PupilCrossCorr_Ephys.mat','Results_PupilCrossCorr_Ephys')
cd([rootFolder delim 'Data'])