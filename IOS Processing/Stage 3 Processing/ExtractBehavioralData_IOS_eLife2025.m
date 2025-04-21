function [chunkdata,evFilter] = ExtractBehavioralData_IOS_eLife2025(data,epoch,dataType,behavior,samplingRate)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% setup variables
eventTime = data.Flags.(behavior).eventTime;
trialDuration = data.notes.trialDuration_sec;
% get the content from data.(dataType)
data = getfield(data,{},dataType,{});
% calculate start/stop times (seconds) for the events
allEpochStarts = eventTime - epoch.offset*ones(size(eventTime));
allEpochEnds = allEpochStarts + epoch.duration*ones(size(eventTime));
% filter out events which are too close to the beginning or end of trials
startFilter = allEpochStarts > 0;
stopFilter = round(allEpochEnds) < trialDuration; % apply "round" to give an extra half second buffer and prevent indexing errors
evFilter = logical(startFilter.*stopFilter);
% convert the starts from seconds to samples, round down to the nearest sample, coerce the value above 1.
epochStarts = max(floor(allEpochStarts(evFilter)*samplingRate),1);
% calculate stops indices using the duration of the epoch, this avoids potential matrix dimension erros caused by differences in rounding when converting from seconds to samples.
sampleDur = round(epoch.duration*samplingRate);
epochStops = epochStarts + sampleDur*ones(size(epochStarts));
% extract the chunk of data from the trial
chunkdata = zeros(sampleDur + 1,length(epochStarts),size(data,1));
for a = 1:length(epochStarts)
    chunkInds = epochStarts(a):epochStops(a);
    chunkdata(:,a,:) = data(:,chunkInds)';
end