function [Results_Derivative_Ephys] = AnalyzeArousalDerivative_Ephys(animalID,group,set,rootFolder,delim,Results_Derivative_Ephys)
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
% find and load Forest_ScoringResults struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% loop variables
hemispheres = {'LH','RH'};
dataTypes = {'HbT'};
% lowpass filter
samplingRate = 30;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        %% Alert
        AlertData = []; % for loop pre-allocation
        zz = 1;
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,allDataFileID] = GetFileInfo_IOS(procDataFileID);
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(allDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % less than 3 minutes of asleep
                load(procDataFileID,'-mat')
                stims = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(stims) == true
                    AlertData{zz,1} = diff(filtfilt(sos,g,ProcData.data.(dataType).(hemisphere)));
                    zz = zz + 1;
                end
            end
        end
        Results_Derivative_Ephys.(group).(animalID).(hemisphere).(dataType).Alert.data = AlertData;
        %% Asleep
        AsleepData = []; % for loop pre-allocation
        zz = 1;
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,asleepDataFileID] = GetFileInfo_IOS(procDataFileID);
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(asleepDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) < 36 % less than 3 minutes of alert
                load(procDataFileID,'-mat')
                stims = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(stims) == true
                    AsleepData{zz,1} = diff(filtfilt(sos,g,ProcData.data.(dataType).(hemisphere)));
                    zz = zz + 1;
                end
            end
        end
        Results_Derivative_Ephys.(group).(animalID).(hemisphere).(dataType).Asleep.data = AsleepData;
        %% All
        AllData = []; % for loop pre-allocation
        zz = 1;
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            load(procDataFileID,'-mat')
            stims = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(stims) == true
                AllData{zz,1} = diff(filtfilt(sos,g,ProcData.data.(dataType).(hemisphere)));
                zz = zz + 1;
            end
        end
        Results_Derivative_Ephys.(group).(animalID).(hemisphere).(dataType).All.data = AllData;
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_Derivative_Ephys.mat','Results_Derivative_Ephys')
cd([rootFolder delim 'Data'])