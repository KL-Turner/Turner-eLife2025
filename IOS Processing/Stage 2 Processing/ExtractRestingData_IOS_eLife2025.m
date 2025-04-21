function [RestData] = ExtractRestingData_IOS_eLife2025(procdataFiles,iteration)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
load(procdataFiles(1,:));
imagingWavelengths = ProcData.notes.imagingWavelengths;
[dataTypes] = DetermineWavelengthDatatypes_IOS_eLife2025(imagingWavelengths,iteration);
% go through each datatype and extract the corresponding data
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    subDataTypes = fieldnames(ProcData.data.(dataType));
    if any(strcmpi(dataType,{'green','lime','blue','red','HbT','HbO','HbR','GCaMP'})) == true
        samplingRate = ProcData.notes.wavelengthSamplingRate;
    else
        samplingRate = ProcData.notes.dsFs;
    end
    % go through subdata types (hemisphere,etc)
    for bb = 1:length(subDataTypes)
        % initialize cell arrays for resting data and other information.
        restVals = cell(size(procdataFiles,1),1);
        eventTimes = cell(size(procdataFiles,1),1);
        durations = cell(size(procdataFiles,1),1);
        stimDistances = cell(size(procdataFiles,1),1);
        fileIDs = cell(size(procdataFiles,1),1);
        fileDates = cell(size(procdataFiles,1),1);
        for cc = 1:size(procdataFiles,1)
            disp(['Extracting ' char(dataType) ' ' char(subDataTypes(bb)) ' rest data from file ' num2str(cc) ' of ' num2str(size(procdataFiles,1)) '...']); disp(' ')
            procdataFile = procdataFiles(cc,:);
            load(procdataFile);
            % get the date and file identifier for the data to be saved with each resting event
            [animalID,fileDate,fileID] = GetFileInfo_IOS_eLife2025(procdataFile);
            % expected number of samples for element of dataType
            trialDuration_sec = ProcData.notes.trialDuration_sec;
            expectedLength = trialDuration_sec*samplingRate;
            % get information about periods of rest from the loaded file
            trialEventTimes = ProcData.flags.rest.eventTime';
            trialStimDistances = ProcData.flags.rest.stimDistance;
            trialDurations = ProcData.flags.rest.duration';
            % initialize cell array for all periods of rest from the loaded file
            trialRestVals = cell(size(trialEventTimes'));
            for d = 1:length(trialEventTimes)
                % extract the whole duration of the resting event. Coerce the
                % start index to values above 1 to preclude rounding to 0.
                startInd = max(floor(trialEventTimes(d)*samplingRate),1);
                % convert the duration from seconds to samples.
                dur = round(trialDurations(d)*samplingRate);
                % get ending index for data chunk. If event occurs at the end of
                % the trial, assume animal whisks as soon as the trial ends and give a 200ms buffer.
                stopInd = min(startInd + dur,expectedLength - round(0.2*samplingRate));
                try
                    % extract data from the trial and add to the cell array for the current loaded file
                    trialRestVals{d} = ProcData.data.(dataTypes{aa}).(subDataTypes{bb})(:,startInd:stopInd);
                catch % some files don't have certain fields. Skip those
                    trialRestVals{d} = [];
                end
            end
            % add all periods of rest to a cell array for all files
            restVals{cc} = trialRestVals';
            % transfer information about resting periods to the new structure
            eventTimes{cc} = trialEventTimes';
            durations{cc} = trialDurations';
            stimDistances{cc} = trialStimDistances';
            fileIDs{cc} = repmat({fileID},1,length(trialEventTimes));
            fileDates{cc} = repmat({fileDate},1,length(trialEventTimes));
        end
        % combine the cells from separate files into a single cell array of all resting periods
        RestData.(dataTypes{aa}).(subDataTypes{bb}).data = [restVals{:}]';
        RestData.(dataTypes{aa}).(subDataTypes{bb}).eventTimes = cell2mat(eventTimes);
        RestData.(dataTypes{aa}).(subDataTypes{bb}).durations = cell2mat(durations);
        RestData.(dataTypes{aa}).(subDataTypes{bb}).stimDistances = [stimDistances{:}]';
        RestData.(dataTypes{aa}).(subDataTypes{bb}).fileIDs = [fileIDs{:}]';
        RestData.(dataTypes{aa}).(subDataTypes{bb}).fileDates = [fileDates{:}]';
        RestData.(dataTypes{aa}).(subDataTypes{bb}).samplingRate = samplingRate;
        RestData.(dataTypes{aa}).(subDataTypes{bb}).trialDuration_sec = trialDuration_sec;
        RestData.(dataTypes{aa}).(subDataTypes{bb}).animalID = animalID;
    end
end
save([animalID '_RestData.mat'],'RestData','-v7.3');