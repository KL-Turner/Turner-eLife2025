function [EventData] = ExtractEventTriggeredData_IOS(procDataFileIDs)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
EventData = [];
epoch.duration = 12;
epoch.offset = 2;
load(procDataFileIDs(1,:));
imagingWavelengths = ProcData.notes.imagingWavelengths;
[dataTypes] = DetermineWavelengthDatatypes_IOS(imagingWavelengths,3);
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    subDataTypes = fieldnames(ProcData.data.(dataType));
    if any(strcmpi(dataType,{'green','lime','blue','red','HbT','HbO','HbR','GCaMP'})) == true
        samplingRate = ProcData.notes.wavelengthSamplingRate;
    else
        samplingRate = ProcData.notes.dsFs;
    end
    temp = struct();
    for bb = 1:size(procDataFileIDs,1)
        % load ProcData File
        filename = procDataFileIDs(bb,:);
        load(filename);
        % get the date and file ID to include in the EventData structure
        [~,fileDate,fileID] = GetFileInfo_IOS(procDataFileIDs(bb,:));
        % get the types of behaviors present in the file (stim, whisk, rest)
        holddata = fieldnames(ProcData.flags);
        behaviorFields = holddata([1,2],1);
        for cc = 1:length(subDataTypes)
            sDT = char(subDataTypes(cc));
            % set the sampling frequency for the dataType
            trialDuration_sec = ProcData.notes.trialDuration_sec;
            % loop over the behaviors present in the file
            for dd = 1:length(behaviorFields)
                % Pre-allocate space for unknown number of events using a
                % 'temporary' structure of cells
                if not(isfield(temp,sDT))
                    temp.(sDT) = [];
                end
                % create behavioral subfields for the temp structure, if needed
                if not(isfield(temp.(sDT),behaviorFields{dd}))
                    subFields = fieldnames(ProcData.flags.(behaviorFields{dd}));
                    blankCell = cell(1,size(procDataFileIDs,1));
                    structVals = cell(size(subFields));
                    structVals(:) = {blankCell};
                    temp.(sDT).(behaviorFields{dd}) = cell2struct(structVals,subFields,1)';
                    temp.(sDT).(behaviorFields{dd}).fileIDs = blankCell;
                    temp.(sDT).(behaviorFields{dd}).fileDates = blankCell;
                    temp.(sDT).(behaviorFields{dd}).data = blankCell;
                end
                % assemble a structure to send to the sub-functions
                fieldName2 = dataType;
                try
                    data = ProcData.data.(fieldName2);
                catch % some files don't have certain fields. Skip those
                    data = [];
                end
                data.Flags = ProcData.flags;
                data.notes = ProcData.notes;
                % extract the data from the epoch surrounding the event
                disp(['Extracting ' dataType ' ' sDT ' event-triggered ' behaviorFields{dd} ' data from file ' num2str(bb) ' of ' num2str(size(procDataFileIDs,1)) '...']); disp(' ');
                try
                    [chunkdata,evFilter] = ExtractBehavioralData_IOS(data,epoch,sDT,behaviorFields{dd},samplingRate);
                catch
                    chunkdata = [];
                    evFilter = [];
                end
                % add epoch details to temp struct
                [temp] = AddEpochInfo_IOS(data,sDT,behaviorFields{dd},temp,fileID,fileDate,evFilter,bb);
                temp.(sDT).(behaviorFields{dd}).data{bb} = chunkdata;
                % add the sampling frequency, assume all Fs are the same for given datatype
                temp.(sDT).(behaviorFields{dd}).samplingRate = {samplingRate};
                temp.(sDT).(behaviorFields{dd}).trialDuration_sec = {trialDuration_sec};
            end
        end
    end
    % convert the temporary stuct into a final structure
    [EventData] = ProcessTempStruct_IOS(EventData,dataType,temp,epoch);
end