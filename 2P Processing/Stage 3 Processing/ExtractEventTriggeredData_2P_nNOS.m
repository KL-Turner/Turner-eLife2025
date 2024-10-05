function [EventData] = ExtractEventTriggeredData_2P_nNOS(mergedDataFileIDs,dataTypes)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%________________________________________________________________________________________________________________________
%
% Purpose: Extracts all event-triggered data from the data using behavioral flags
%________________________________________________________________________________________________________________________

EventData = [];
epoch.duration = 12;
epoch.offset = 2;
% Control for dataTypes as string
if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end
for dT = 1:length(dataTypes)
    temp = struct();
    dataType = dataTypes{dT};
    if strcmp(dataType,'vesselDiameter') == true
        subDataTypes = {'data'};
    elseif strcmp(dataType,'EMG') == true
        subDataTypes = {'data'};
    else
        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
    end
    for f = 1:size(mergedDataFileIDs,1)
        % Load MergedData File
        mergedDataFileID = mergedDataFileIDs(f,:);
        load(mergedDataFileID);
        % Get the date and file ID to include in the EventData structure
        [animalID,~,fileDate,fileID,imageID,vesselID] = GetFileInfo2_2P_nNOS(mergedDataFileIDs(f,:));
        % Get the types of behaviors present in the file (stim,whisk,rest)
        holdData = fieldnames(MergedData.flags);
        behaviorFields = holdData([1,2],1);      
        % Sampling frequency for element of dataTypes
        if strcmp(dataType,'vesselDiameter')
            samplingRate = floor(MergedData.notes.p2Fs);
        else
            samplingRate = floor(MergedData.notes.dsFs);
        end
        for c = 1:length(subDataTypes)
            sDT = char(subDataTypes(c));
            % Loop over the behaviors present in the file
            for bF = 1:length(behaviorFields)
                if not(isfield(temp,sDT))
                    temp.(sDT) = [];
                end
                % Create behavioral subfields for the temp structure, if needed
                if not(isfield(temp.(sDT),behaviorFields{bF}))
                    subFields = fieldnames(MergedData.flags.(behaviorFields{bF}));
                    blankCell = cell(1,size(mergedDataFileIDs,1));
                    structVals = cell(size(subFields));
                    structVals(:) = {blankCell};
                    temp.(sDT).(behaviorFields{bF}) = cell2struct(structVals,subFields,1)';
                    temp.(sDT).(behaviorFields{bF}).fileIDs = blankCell;
                    temp.(sDT).(behaviorFields{bF}).imageIDs = blankCell;
                    temp.(sDT).(behaviorFields{bF}).fileDates = blankCell;
                    temp.(sDT).(behaviorFields{bF}).data = blankCell;
                    temp.(sDT).(behaviorFields{bF}).vesselIDs = blankCell;
                end
                % Assemble a structure to send to the sub-functions
                fieldName2 = dataType;
                data = MergedData.data.(fieldName2);
                data.flags = MergedData.flags;
                data.notes = MergedData.notes;
                % Extract the data from the epoch surrounding the event
                disp(['Extracting event-triggered ' dataType ' ' sDT ' ' behaviorFields{bF} ' data from file ' num2str(f) ' of ' num2str(size(mergedDataFileIDs,1)) '...']); disp(' ');
                [chunkData,evFilter] = ExtractBehavioralData_2P_nNOS(data,epoch,sDT,samplingRate,behaviorFields{bF});
                % Add epoch details to temp struct
                [temp] = AddEpochInfo_2P_nNOS(data,sDT,behaviorFields{bF},temp,fileID,fileDate,imageID,vesselID,evFilter,f);
                temp.(sDT).(behaviorFields{bF}).data{f} = chunkData;
            end
        end
    end
    % Convert the temporary stuct into a final structure
    [EventData] = ProcessTempStruct_2P_nNOS(EventData,dataType,temp,epoch);
end
save([animalID '_EventData.mat'],'EventData');

end

function [chunkData,evFilter] = ExtractBehavioralData_2P_nNOS(data,epoch,dataType,samplingRate,behavior)
% Setup variables
eventTimes = data.flags.(behavior).eventTime;
trialDuration = data.notes.trialDuration_Sec;
% Get the content from data.(dataType)
data = getfield(data,{},dataType,{});
% Calculate start/stop times (seconds) for the events
allEpochStarts = eventTimes - epoch.offset*ones(size(eventTimes));
allEpochEnds = allEpochStarts + epoch.duration*ones(size(eventTimes));
% Filter out events which are too close to the beginning or end of trials
startFilter = allEpochStarts > 0;
stopFilter = round(allEpochEnds) < trialDuration; % Apply "round" to give an extra half second buffer and prevent indexing errors
evFilter = logical(startFilter.*stopFilter);
% Convert the starts from seconds to samples, round down to the nearest
% sample, coerce the value above 1.
epochStarts = max(floor(allEpochStarts(evFilter)*samplingRate),1);
% Calculate stops indices using the duration of the epoch, this avoids
% potential matrix dimension erros caused by differences in rounding when
% converting from seconds to samples.
sampleDur = round(epoch.duration*samplingRate);
epochStops = epochStarts + sampleDur*ones(size(epochStarts));
% Extract the chunk of data from the trial
chunkData = zeros(sampleDur + 1,length(epochStarts),size(data,1));
for eS = 1:length(epochStarts)
    chunkInds = epochStarts(eS):epochStops(eS);
    chunkData(:,eS,:) = data(:,chunkInds)';
end

end

function [temp] = AddEpochInfo_2P_nNOS(data,dataType,behavior,temp,fileID,fileDate,imageID,vesselID,evFilter,f)
% Get the field names for each behavior
fields = fieldnames(data.flags.(behavior));
% Filter out the events which are too close to the trial edge
for flds = 1:length(fields)
    field = fields{flds};
    temp.(dataType).(behavior).(field){f} = data.flags.(behavior).(field)(evFilter,:)';
end
% Tag each event with the file ID, arrange cell array horizontally for
% later processing.
temp.(dataType).(behavior).fileIDs{f} = repmat({fileID},1,sum(evFilter));
temp.(dataType).(behavior).fileDates{f} = repmat({fileDate},1,sum(evFilter));
temp.(dataType).(behavior).imageIDs{f} = repmat({imageID},1,sum(evFilter));
temp.(dataType).(behavior).vesselIDs{f} = repmat({vesselID},1,sum(evFilter));

end

function [EventData] = ProcessTempStruct_2P_nNOS(EventData,dataType,temp,epoch)
% Get the dataTypes from temp
dTs = fieldnames(temp);
for a = 1:length(dTs)
    dT = dTs{a};
    % Get dataType names
    behaviorFields = fieldnames(temp.(dT));
    % Intialize Behavior fields of the dataType sub-structure
    structArray2 = cell(size(behaviorFields));
    EventData.(dataType).(dT) = cell2struct(structArray2,behaviorFields,1);
    for bF = 1:length(behaviorFields)
        behavior = behaviorFields{bF};
        % Get Behavior names
        eventFields = fieldnames(temp.(dT).(behavior));
        % Initialize Event fields for the Behavior sub-structure
        structArray3 = cell(size(eventFields));
        EventData.(dataType).(dT).(behavior) = cell2struct(structArray3,eventFields,1);
        for eF = 1:length(eventFields)
            evField = eventFields{eF};
            transferArray = [temp.(dT).(behavior).(evField){:}];
            EventData.(dataType).(dT).(behavior).(evField) = permute(transferArray,unique([2,1,ndims(transferArray)],'stable'));
        end
        EventData.(dataType).(dT).(behavior).epoch = epoch;
    end
end

end
