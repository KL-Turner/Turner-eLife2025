function [Rest] = GetRestData_IOS_nNOS(ProcData)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% setup
whiskerSamplingRate = ProcData.notes.dsFs;
forceSensorSamplingRate = ProcData.notes.dsFs;
% get stimulation times
[stimTimes] = GetStimTimes_IOS_nNOS(ProcData);
% recalculate linked binarized wwf without omitting any possible whisks
modBinarizedWhiskers = ProcData.data.whiskerAngle.binarization;
modBinarizedWhiskers([1,end]) = 1;
modBinarizedForceSensor = ProcData.data.forceSensor.binarization;
modBinarizedForceSensor([1,end]) = 1;
linkThresh = 0.5; % seconds
breakThresh = 0; % seconds
binWhiskerAngle = LinkBinaryEvents_IOS_nNOS(gt(modBinarizedWhiskers,0),[linkThresh breakThresh]*whiskerSamplingRate);
binForceSensor = LinkBinaryEvents_IOS_nNOS(modBinarizedForceSensor,[linkThresh breakThresh]*forceSensorSamplingRate);
% combine binWhiskerAngle, binForceSensor, and stimTimes, to find periods of rest
sampleVec = 1:length(binWhiskerAngle);
whiskHigh = sampleVec(binWhiskerAngle)/whiskerSamplingRate;
dsBinarizedWhiskers = zeros(size(binForceSensor));
% find Bin_wwf == 1. Convert indexes into pswf time. Coerce converted indexes
% between 1 and length(Bin_pswf). Take only unique values.
dsInds = min(max(round(whiskHigh*forceSensorSamplingRate),1),length(binForceSensor));
dsBinarizedWhiskers(unique(dsInds)) = 1;
% combine binarized whisking and body movement
wfBin = logical(min(dsBinarizedWhiskers + binForceSensor,1));
samplingRate = forceSensorSamplingRate;
% add stim times into the Bin_wf
stimInds = round(stimTimes*samplingRate);
wfBin(stimInds) = 1;
% find index for end of whisking event
edge = diff(wfBin);
samples = find([not(wfBin(1)),edge < 0]);
stops = samples/samplingRate;
% identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
sampleVec = 1:length(logical(wfBin));
highSamples = unique([1,sampleVec(wfBin),sampleVec(end)]);
lowSamples = unique([1,sampleVec(not(wfBin)),sampleVec(end)]);
% calculate the number of samples between consecutive high/low samples.
dHigh = diff(highSamples);
dLow = diff(lowSamples);
% identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
restLength = dHigh(dHigh > 1);
restDur = restLength/samplingRate;
whiskLength = dLow(dLow > 1);
whiskDur = whiskLength/samplingRate;
% control for the beginning/end of the trial to correctly map rests/whisks
if not(wfBin(2))
    whiskDur = [NaN,whiskDur];
end
% control for the beginning/end of the trial to correctly map rests/whisks
if wfBin(end - 1)
    whiskDur(end) = [];
end
% calculate the time to the closest stim
if isempty(stimTimes)
    stimTimes = 0;
end
stimMat = ones(length(samples),1)*stimTimes;
restMat = samples'*ones(1,length(stimTimes))/samplingRate;
stimTimeElapsed = abs(restMat - stimMat);
% convert to cell
stimTimeCell = mat2cell(stimTimeElapsed,ones(length(samples),1));
% compile into a structure
Rest.eventTime = stops';
Rest.duration = restDur';
Rest.stimDistance = stimTimeCell;
Rest.whiskDuration = whiskDur';