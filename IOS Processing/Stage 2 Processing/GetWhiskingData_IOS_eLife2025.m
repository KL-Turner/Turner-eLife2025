function [Whisk] = GetWhiskingData_IOS_eLife2025(ProcData,binWhiskerAngle)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% setup
whiskerSamplingRate = ProcData.notes.dsFs;
forceSensorSamplingRate = ProcData.notes.dsFs;
% get Stim Times
[stimTimes] = GetStimTimes_IOS_eLife2025(ProcData);
% find the starts of whisking
whiskEdge = diff(binWhiskerAngle);
whiskSamples = find(whiskEdge > 0);
whiskStarts = whiskSamples/whiskerSamplingRate;
% classify each whisking event by duration, whisking intensity, rest durations
sampleVec = 1:length(binWhiskerAngle);
% identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
highSamples = unique([1,sampleVec(binWhiskerAngle),sampleVec(end)]);
lowSamples = unique([1,sampleVec(not(binWhiskerAngle)),sampleVec(end)]);
% calculate the number of samples between consecutive high/low samples.
dHigh = diff(highSamples);
dLow = diff(lowSamples);
% identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
restLength = dHigh(dHigh > 1);
whiskLength = dLow(dLow > 1);
restDur = restLength/whiskerSamplingRate;
whiskDur = whiskLength/whiskerSamplingRate;
% control for the beginning/end of the trial to correctly map rests/whisks
% onto the whisk_starts.
if binWhiskerAngle(1)
    whiskDur(1) = [];
    whiskLength(1) = [];
end
if not(binWhiskerAngle(end))
    restDur(end) = [];
end
% calculate the whisking intensity -> sum(ProcData.Bin_wwf)/sum(Bin_wwf)
% over the duration of the whisk. Calculate the movement intensity over the same interval.
whiskInt = zeros(size(whiskStarts));
movementInt = zeros(size(whiskStarts));
for wS = 1:length(whiskSamples)
    % whisking intensity
    whiskInds = whiskSamples(wS):whiskSamples(wS) + whiskLength(wS);
    whiskInt(wS) = sum(ProcData.data.whiskerAngle.binarization(whiskInds))/numel(whiskInds);
    % movement intensity
    movementStart = round(whiskStarts(wS)*forceSensorSamplingRate);
    movementDur = round(whiskDur(wS)*forceSensorSamplingRate);
    movementInds = max(movementStart,1):min(movementStart + movementDur,length(ProcData.data.forceSensor.binarization));
    movementInt(wS) = sum(ProcData.data.forceSensor.binarization(movementInds))/numel(movementInds);
end
% calculate the time to the closest stim
% if no stim occurred during the trial, store 0 as a place holder.
if isempty(stimTimes)
    stimTimes = 0;
end
stimMat = ones(length(whiskSamples),1)*stimTimes;
whiskMat = whiskSamples'*ones(1,length(stimTimes))/whiskerSamplingRate;
stimTimeElapsed = abs(whiskMat - stimMat);
% convert to cell
stimTimeCell = mat2cell(stimTimeElapsed,ones(length(whiskStarts),1));
% error handle
if length(restDur) ~= length(whiskDur)
    disp('Error in GetWhiskdata! The number of whisks does not equal the number of rests...'); disp(' ')
    keyboard;
end
% compile into final structure
Whisk.eventTime = whiskStarts';
Whisk.duration = whiskDur';
Whisk.restTime = restDur';
Whisk.whiskScore = whiskInt';
Whisk.movementScore = movementInt';
Whisk.stimDistance = stimTimeCell;