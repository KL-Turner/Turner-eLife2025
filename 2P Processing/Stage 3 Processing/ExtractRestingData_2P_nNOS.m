function [RestData] = ExtractRestingData_2P(mergedDataFileIDs,dataTypes)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
% Purpose: Extracts all resting data periods from the data using behavioral flags
%________________________________________________________________________________________________________________________

% control for singular input
if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end
% go through each datatype and extract the corresponding data
for dT = 1:length(dataTypes)
    dataType = dataTypes(dT);
    if strcmp(dataType,'vesselDiameter') == true
        subDataTypes = {'data'};
    elseif strcmp(dataType,'EMG') == true
        subDataTypes = {'data'};
    else
        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
    end
    % go through subdata types (hemisphere,etc)
    for b = 1:length(subDataTypes)
        subDataType = subDataTypes{1,b};
        % Initialize cell arrays for resting data and other information.
        restVals = cell(size(mergedDataFileIDs,1),1);
        eventTimes = cell(size(mergedDataFileIDs,1),1);
        durations = cell(size(mergedDataFileIDs,1),1);
        puffDistances = cell(size(mergedDataFileIDs,1),1);
        fileIDs = cell(size(mergedDataFileIDs,1),1);
        fileDates = cell(size(mergedDataFileIDs,1),1);
        vesselIDs = cell(size(mergedDataFileIDs,1),1);
        for f = 1:size(mergedDataFileIDs,1)
            disp(['Gathering rest ' char(dataType) ' data from file ' num2str(f) ' of ' num2str(size(mergedDataFileIDs,1)) '...']); disp(' ')
            mergedDataFileID = mergedDataFileIDs(f,:);
            load(mergedDataFileID);
            % Get the date and file identifier for the data to be saved with each resting event
            [animalID,~,fileDate,fileID,~,vesselID] = GetFileInfo2_2P(mergedDataFileID);
            % Sampling frequency for element of dataTypes
            if strcmp(dataType,'vesselDiameter')
                Fs = floor(MergedData.notes.p2Fs);
            else
                Fs = floor(MergedData.notes.dsFs);
            end            
            % Expected number of samples for element of dataType
            expectedLength = MergedData.notes.trialDuration_Sec*Fs;
            % Get information about periods of rest from the loaded file
            trialEventTimes = MergedData.flags.rest.eventTime';
            trialPuffDistances = MergedData.flags.rest.puffDistance;
            trialDurations = MergedData.flags.rest.duration';
            % Initialize cell array for all periods of rest from the loaded file
            trialRestVals = cell(size(trialEventTimes'));
            for tET = 1:length(trialEventTimes)
                % Extract the whole duration of the resting event. Coerce the
                % start index to values above 1 to preclude rounding to 0.
                startInd = max(floor(trialEventTimes(tET)*Fs),1);                
                % Convert the duration from seconds to samples.
                dur = round(trialDurations(tET)*Fs);                
                % Get ending index for data chunk. If event occurs at the end of
                % the trial, assume animal whisks as soon as the trial ends and
                % give a 200ms buffer.
                stopInd = min(startInd + dur,expectedLength - round(0.2*Fs));                
                % Extract data from the trial and add to the cell array for the current loaded file
                trialRestVals{tET} = MergedData.data.(dataTypes{dT}).(subDataType)(:,startInd:stopInd);
            end
            % Add all periods of rest to a cell array for all files
            restVals{f} = trialRestVals';            
            % Transfer information about resting periods to the new structure
            eventTimes{f} = trialEventTimes';
            durations{f} = trialDurations';
            puffDistances{f} = trialPuffDistances';
            fileIDs{f} = repmat({fileID},1,length(trialEventTimes));
            fileDates{f} = repmat({fileDate},1,length(trialEventTimes));
            vesselIDs{f} = repmat({vesselID},1,length(trialEventTimes));
        end
        % Combine the cells from separate files into a single cell array of all resting periods
        RestData.(dataTypes{dT}).(subDataType).data = [restVals{:}]';
        RestData.(dataTypes{dT}).(subDataType).eventTimes = cell2mat(eventTimes);
        RestData.(dataTypes{dT}).(subDataType).durations = cell2mat(durations);
        RestData.(dataTypes{dT}).(subDataType).puffDistances = [puffDistances{:}]';
        RestData.(dataTypes{dT}).(subDataType).fileIDs = [fileIDs{:}]';
        RestData.(dataTypes{dT}).(subDataType).fileDates = [fileDates{:}]';
        RestData.(dataTypes{dT}).(subDataType).vesselIDs = [vesselIDs{:}]';
        RestData.(dataTypes{dT}).(subDataType).samplingRate = Fs;
    end
end
save([animalID '_RestData.mat'],'RestData');

end
