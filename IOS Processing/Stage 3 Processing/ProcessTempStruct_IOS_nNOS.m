function [EventData] = ProcessTempStruct_IOS_nNOS(EventData,dataType,temp,epoch)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% get the dataTypes from temp
dTs = fieldnames(temp);
for a = 1:length(dTs)
    dT = dTs{a};
    % get dataType names
    behaviorFields = fieldnames(temp.(dT));
    % intialize Behavior fields of the dataType sub-structure
    structArray2 = cell(size(behaviorFields));
    EventData.(dataType).(dT) = cell2struct(structArray2,behaviorFields,1);
    for b = 1:length(behaviorFields)
        behavior = behaviorFields{b};
        % get Behavior names
        eventFields = fieldnames(temp.(dT).(behavior));
        % initialize Event fields for the Behavior sub-structure
        structArray3 = cell(size(eventFields));
        EventData.(dataType).(dT).(behavior) = cell2struct(structArray3,eventFields,1);
        for c = 1:length(eventFields)
            evField = eventFields{c};
            transferArray = [temp.(dT).(behavior).(evField){:}];
            EventData.(dataType).(dT).(behavior).(evField) = permute(transferArray,unique([2,1,ndims(transferArray)],'stable'));
        end
        EventData.(dataType).(dT).(behavior).epoch = epoch;
    end
end