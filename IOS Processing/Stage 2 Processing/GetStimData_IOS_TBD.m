function [Stim] = GetStimData_IOS(ProcData)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% get stimulation times
whiskerSamplingRate = ProcData.notes.dsFs;
forceSensorSamplingRate = ProcData.notes.dsFs;
stimTimes = GetStimTimes_IOS(ProcData);
trialDuration = ProcData.notes.trialDuration_sec;
% set time intervals for calculation of the whisk scores
preTime = 1;
postTime = 1;
% get stimer IDs
stimNames = fieldnames(ProcData.data.stimulations);
Stim.solenoidName = cell(length(stimTimes),1);
Stim.eventTime = zeros(length(stimTimes),1);
Stim.whiskScore_Pre = zeros(length(stimTimes),1);
Stim.whiskScore_Post = zeros(length(stimTimes),1);
Stim.movementScore_Pre = zeros(length(stimTimes),1);
Stim.movementScore_Post = zeros(length(stimTimes),1);
j = 1;
for sN = 1:length(stimNames)
    solStimTimes = ProcData.data.stimulations.(stimNames{sN});
    for spT = 1:length(solStimTimes)
        if trialDuration - solStimTimes(spT) <= postTime
            disp(['Stim at time: ' solStimTimes(spT) ' is too close to trial end'])
            continue;
        end
        % set indexes for pre and post periods
        wStimInd = round(solStimTimes(spT)*whiskerSamplingRate);
        mStimInd = round(solStimTimes(spT)*forceSensorSamplingRate);
        wPreStart = max(round((solStimTimes(spT) - preTime)*whiskerSamplingRate),1);
        mPreStart = max(round((solStimTimes(spT) - preTime)*forceSensorSamplingRate),1);
        wPostEnd = round((solStimTimes(spT) + postTime)*whiskerSamplingRate);
        mPostEnd = round((solStimTimes(spT) + postTime)*forceSensorSamplingRate);
        % calculate the percent of the pre-stim time that the animal moved or whisked
        whiskScorePre = sum(ProcData.data.whiskerAngle.binarization(wPreStart:wStimInd))/(preTime*whiskerSamplingRate);
        whiskScorePost = sum(ProcData.data.whiskerAngle.binarization(wStimInd:wPostEnd))/(postTime*whiskerSamplingRate);
        moveScorePre = sum(ProcData.data.forceSensor.binarization(mPreStart:mStimInd))/(preTime*forceSensorSamplingRate);
        moveScorePost = sum(ProcData.data.forceSensor.binarization(mStimInd:mPostEnd))/(postTime*forceSensorSamplingRate);
        % add to Stim structure
        Stim.solenoidName{j} = stimNames{sN};
        Stim.eventTime(j) = solStimTimes(spT)';
        Stim.whiskScore_Pre(j) = whiskScorePre';
        Stim.whiskScore_Post(j) = whiskScorePost';
        Stim.movementScore_Pre(j) = moveScorePre';
        Stim.movementScore_Post(j) = moveScorePost';
        j = j + 1;
    end
end
% calculate the time to the closest stim, omit comparison of stim to itself
stimMat = ones(length(stimTimes),1)*stimTimes;
timeElapsed = abs(nonzeros(stimMat - stimMat'));
% if no other stim occurred during the trial, store 0 as a place holder.
if isempty(timeElapsed)
    stimTimeElapsed = 0;
else
    % if not empty, Reshape the array to compensate for nonzeros command
    stimTimeElapsed = reshape(timeElapsed,numel(stimTimes) - 1,numel(stimTimes));
end
% convert to cell and add to struct, if length of Stim_Times = 0, coerce to
% 1 to accommodate the NaN entry.
stimTimeCell = mat2cell(stimTimeElapsed',ones(max(length(stimTimes),1),1));
Stim.StimDistance = stimTimeCell;