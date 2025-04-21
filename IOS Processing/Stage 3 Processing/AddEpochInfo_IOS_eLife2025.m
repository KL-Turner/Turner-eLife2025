function [temp] = AddEpochInfo_IOS_eLife2025(data,dataType,behavior,temp,fileID,fileDate,evFilter,f)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% get the field names for each behavior
fields = fieldnames(data.Flags.(behavior));
% filter out the events which are too close to the trial edge
for a = 1:length(fields)
    field = fields{a};
    temp.(dataType).(behavior).(field){f} = data.Flags.(behavior).(field)(evFilter,:)';
end
% tag each event with the file ID, arrange cell array horizontally for later processing.
temp.(dataType).(behavior).fileIDs{f} = repmat({fileID},1,sum(evFilter));
temp.(dataType).(behavior).fileDates{f} = repmat({fileDate},1,sum(evFilter));
