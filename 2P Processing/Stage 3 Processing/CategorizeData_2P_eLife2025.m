function [] = CategorizeData_2P_eLife2025(mergedDataFileID,stimulationType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Originally written by Aaron T. Winder
%________________________________________________________________________________________________________________________
%
% Purpose: Catagorizes data based on behavioral flags from whisking/movement events
%________________________________________________________________________________________________________________________

% Load and Setup
disp(['Categorizing data for: ' mergedDataFileID]); disp(' ')
load(mergedDataFileID)
% condense pulsed stimulations
if strcmp(stimulationType,'pulse') == true
    stimulationFields = fieldnames(MergedData.data.solenoids);
    for aa = 1:length(stimulationFields)
        stimulationTimes = MergedData.data.solenoids.(stimulationFields{aa,1});
        MergedData.data.solenoidsOriginal.(stimulationFields{aa,1}) = stimulationTimes;
        condensedStimulationTimes = [];
        cc = 1;
        if isempty(stimulationTimes) == false
            for bb = 1:length(stimulationTimes)
                if bb == 1
                    condensedStimulationTimes(1,bb) = stimulationTimes(1,bb);
                    cc = cc + 1;
                else
                    timeDifference = stimulationTimes(1,bb) - stimulationTimes(1,bb - 1);
                    if timeDifference > 1 % remove stimulations that are closer than 1 second to the previous
                        condensedStimulationTimes(1,cc) = stimulationTimes(1,bb);
                        cc = cc + 1;
                    end
                end
            end
            MergedData.data.solenoids.(stimulationFields{aa,1}) = condensedStimulationTimes;
        end
    end
end
whiskerSamplingRate = MergedData.notes.dsFs;
% Process binary whisking waveform to detect whisking events
% Setup parameters for link_binary_events
linkThresh = 0.5;   % seconds, Link events < 0.5 seconds apart
breakThresh = 0;   % seconds changed by atw on 2/6/18 from 0.07
% Assume that whisks at the beginning/end of trial continue outside of the
% trial time. This will link any event occurring within "link_thresh"
% seconds to the beginning/end of the trial rather than assuming that it is
% a new/isolated event.
modBinWhiskers = MergedData.data.binWhiskerAngle;
% Link the binarized whisking for use in GetWhiskingData function
binWhiskers = LinkBinaryEvents_2P_eLife2025(gt(modBinWhiskers,0),[linkThresh breakThresh]*whiskerSamplingRate);
% Added 2/6/18 with atw. Code throws errors if binWhiskers(1)=1 and binWhiskers(2) = 0, or if 
% binWhiskers(1) = 0 and binWhiskers(2) = 1. This happens in GetWhiskingData because starts of 
% whisks are detected by taking the derivative of binWhiskers. Purpose of following lines is to 
% handle trials where the above conditions occur and avoid difficult dimension errors.
if binWhiskers(1) == 0 && binWhiskers(2) == 1
    binWhiskers(1) = 1;
elseif binWhiskers(1) == 1 && binWhiskers(2) == 0
    binWhiskers(1) = 0;
end
if binWhiskers(end) == 0 && binWhiskers(end - 1) == 1
    binWhiskers(end) = 1;
elseif binWhiskers(end) == 1 && binWhiskers(end - 1) == 0
    binWhiskers(end) = 0;
end
% Categorize data by behavior
% Retrieve details on whisking events
[MergedData.flags.whisk] = GetWhiskingData_2P_eLife2025(MergedData,binWhiskers);
% Identify and separate resting data
[MergedData.flags.stim] = GetStimData_2P_eLife2025(MergedData);
% Identify and separate resting data
[MergedData.flags.rest] = GetRestData_2P_eLife2025(MergedData);
% Save MergedData structure
save(mergedDataFileID,'MergedData');
end

function [puffTimes] = GetPuffTimes_2P_eLife2025(MergedData)
solNames = fieldnames(MergedData.data.solenoids);
puffList = cell(1,length(solNames));
for sN = 1:length(solNames)
    puffList{sN} = MergedData.data.solenoids.(solNames{sN});
end
puffTimes = cell2mat(puffList);
end

function [Stim] = GetStimData_2P_eLife2025(MergedData)
% Setup
whiskerSamplingRate = MergedData.notes.dsFs;
forceSensorSamplingRate = MergedData.notes.dsFs;
puffTimes = GetPuffTimes_2P_eLife2025(MergedData);
trialDuration = MergedData.notes.trialDuration_Sec;
% Set time intervals for calculation of the whisk scores
preTime = 1;
postTime = 1;
% Get puffer IDs
solNames = fieldnames(MergedData.data.solenoids);
Stim.solenoidName = cell(length(puffTimes),1);
Stim.eventTime = zeros(length(puffTimes),1);
Stim.whiskScore_Pre = zeros(length(puffTimes),1);
Stim.whiskScore_Post = zeros(length(puffTimes),1);
Stim.movementScore_Pre = zeros(length(puffTimes),1);
Stim.movementScore_Post = zeros(length(puffTimes),1);
ii = 1;
for sN = 1:length(solNames)
    solPuffTimes = MergedData.data.solenoids.(solNames{sN});
    for spT = 1:length(solPuffTimes) 
        if trialDuration - solPuffTimes(spT) <= postTime
            disp(['Puff at time: ' solPuffTimes(spT) ' is too close to trial end'])
            continue;
        end
        % Set indexes for pre and post periods
        wPuffInd = round(solPuffTimes(spT)*whiskerSamplingRate);
        mPuffInd = round(solPuffTimes(spT)*forceSensorSamplingRate);
        wPreStart = max(round((solPuffTimes(spT) - preTime)*whiskerSamplingRate),1);
        mPreStart = max(round((solPuffTimes(spT) - preTime)*forceSensorSamplingRate),1);
        wPostEnd = round((solPuffTimes(spT) + postTime)*whiskerSamplingRate);
        mPostEnd = round((solPuffTimes(spT) + postTime)*forceSensorSamplingRate);
        % Calculate the percent of the pre-stim time that the animal moved or whisked
        whiskScorePre = sum(MergedData.data.binWhiskerAngle(wPreStart:wPuffInd))/(preTime*whiskerSamplingRate);
        whiskScorePost = sum(MergedData.data.binWhiskerAngle(wPuffInd:wPostEnd))/(postTime*whiskerSamplingRate);
        moveScorePre = sum(MergedData.data.binForceSensorL(mPreStart:mPuffInd))/(preTime*forceSensorSamplingRate);
        moveScorePost = sum(MergedData.data.binForceSensorL(mPuffInd:mPostEnd))/(postTime*forceSensorSamplingRate);     
        % Add to Stim structure
        Stim.solenoidName{ii} = solNames{sN};
        Stim.eventTime(ii) = solPuffTimes(spT)';
        Stim.whiskScore_Pre(ii) = whiskScorePre';
        Stim.whiskScore_Post(ii) = whiskScorePost';
        Stim.movementScore_Pre(ii) = moveScorePre'; 
        Stim.movementScore_Post(ii) = moveScorePost';
        ii = ii + 1;
    end
end
% Calculate the time to the closest puff, omit comparison of puff to itself
% (see nonzeros)
puffMat = ones(length(puffTimes),1)*puffTimes;
timeElapsed = abs(nonzeros(puffMat - puffMat'));
% If no other puff occurred during the trial, store 0 as a place holder.
if isempty(timeElapsed)
    puffTimeElapsed = 0;
else
% if not empty, Reshape the array to compensate for nonzeros command
    puffTimeElapsed = reshape(timeElapsed,numel(puffTimes) - 1,numel(puffTimes));
end
% Convert to cell and add to struct, if length of Puff_Times = 0, coerce to
% 1 to accommodate the NaN entry.
puffTimeCell = mat2cell(puffTimeElapsed',ones(max(length(puffTimes),1),1));
Stim.PuffDistance = puffTimeCell;
end

function [Whisk] = GetWhiskingData_2P_eLife2025(MergedData, binarizedWhiskers)
% Setup
whiskerSamplingRate = MergedData.notes.dsFs;
forceSensorSamplingRate = MergedData.notes.dsFs;
% Get Puff Times
[puffTimes] = GetPuffTimes_2P_eLife2025(MergedData);
% Find the starts of whisking
whiskEdge = diff(binarizedWhiskers);
whiskSamples = find(whiskEdge > 0);
whiskStarts = whiskSamples/whiskerSamplingRate;
% Classify each whisking event by duration, whisking intensity, rest durations
sampleVec = 1:length(binarizedWhiskers); 
% Identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
highSamples = unique([1,sampleVec(binarizedWhiskers),sampleVec(end)]); 
lowSamples = unique([1,sampleVec(not(binarizedWhiskers)),sampleVec(end)]);
% Calculate the number of samples between consecutive high/low samples.
dHigh = diff(highSamples);
dLow = diff(lowSamples);
% Identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
restLength = dHigh(dHigh > 1);
whiskLength = dLow(dLow > 1);
restDur = restLength/whiskerSamplingRate;
whiskDur = whiskLength/whiskerSamplingRate;
% Control for the beginning/end of the trial to correctly map rests/whisks
% onto the whisk_starts.
if binarizedWhiskers(1)
    whiskDur(1) = [];
    whiskLength(1) = [];
end
if not(binarizedWhiskers(end))
    restDur(end) = [];
end
% Calculate the whisking intensity -> sum(MergedData.Bin_wwf)/sum(Bin_wwf)
% over the duration of the whisk. Calculate the movement intensity over the same interval.
whiskInt = zeros(size(whiskStarts));
movementInt = zeros(size(whiskStarts));
for wS = 1:length(whiskSamples)
    % Whisking intensity
    whiskInds = whiskSamples(wS):whiskSamples(wS) + whiskLength(wS);
    whiskInt(wS) = sum(MergedData.data.binWhiskerAngle(whiskInds))/numel(whiskInds);
    % Movement intensity
    movementStart = round(whiskStarts(wS)*forceSensorSamplingRate);
    movementDur = round(whiskDur(wS)*forceSensorSamplingRate);
    movementInds = max(movementStart,1):min(movementStart + movementDur,length(MergedData.data.binForceSensorM));
    movementInt(wS) = sum(MergedData.data.binForceSensorM(movementInds))/numel(movementInds);
end
% Calculate the time to the closest puff
% If no puff occurred during the trial, store 0 as a place holder.
if isempty(puffTimes)
    puffTimes = 0;
end
puffMat = ones(length(whiskSamples),1)*puffTimes;
whiskMat = whiskSamples'*ones(1,length(puffTimes))/whiskerSamplingRate;
puffTimeElapsed = abs(whiskMat - puffMat);
% Convert to cell
puffTimeCell = mat2cell(puffTimeElapsed,ones(length(whiskStarts),1));
% Error handle
if length(restDur) ~= length(whiskDur)
    disp('Error in GetWhiskData! The number of whisks does not equal the number of rests...'); disp(' ')
    keyboard;
end
% Compile into final structure
Whisk.eventTime = whiskStarts';
Whisk.duration = whiskDur';
Whisk.restTime = restDur';
Whisk.whiskScore = whiskInt';
Whisk.movementScore = movementInt';
Whisk.puffDistance = puffTimeCell;
end

function [Rest] = GetRestData_2P_eLife2025(MergedData)
% Setup
whiskerSamplingRate = MergedData.notes.dsFs;
forceSensorSamplingRate = MergedData.notes.dsFs;
% Get stimulation times
[puffTimes] = GetPuffTimes_2P_eLife2025(MergedData);
% Recalculate linked binarized wwf without omitting any possible whisks,
% this avoids inclusion of brief whisker movements in periods of rest.
% Assume that whisks at the beginning/end of trial continue outside of the
% trial time. This will link any event occurring within "link_thresh"
% seconds to the beginning/end of the trial rather than assuming that it is a new/isolated event.
modBinarizedWhiskers = MergedData.data.binWhiskerAngle;
modBinarizedWhiskers([1, end]) = 1;
modBinarizedForceSensor = MergedData.data.binForceSensorM;
modBinarizedForceSensor([1,end]) = 1;
linkThresh = 0.5;   % seconds
breakThresh = 0;   % seconds
binarizedWhiskers = LinkBinaryEvents_2P_eLife2025(gt(modBinarizedWhiskers,0),[linkThresh breakThresh]*whiskerSamplingRate);
binarizedForceSensor = LinkBinaryEvents_2P_eLife2025(modBinarizedForceSensor,[linkThresh breakThresh]*forceSensorSamplingRate);
% Combine binarizedWhiskers, binarizedForceSensor, and puffTimes, to find periods of rest. 
% Downsample bin_wwf to match length of bin_pswf
sampleVec = 1:length(binarizedWhiskers); 
whiskHigh = sampleVec(binarizedWhiskers)/whiskerSamplingRate;
dsBinarizedWhiskers = zeros(size(binarizedForceSensor));
% Find Bin_wwf == 1. Convert indexes into pswf time. Coerce converted indexes
% between 1 and length(Bin_pswf). Take only unique values.
dsInds = min(max(round(whiskHigh*forceSensorSamplingRate),1),length(binarizedForceSensor));
dsBinarizedWhiskers(unique(dsInds)) = 1;
% Combine binarized whisking and body movement
wfBin = logical(min(dsBinarizedWhiskers + binarizedForceSensor,1));
Fs = forceSensorSamplingRate;
% Add puff times into the Bin_wf
puffInds = round(puffTimes*Fs);
wfBin(puffInds) = 1;
% Find index for end of whisking event
edge = diff(wfBin);
samples = find([not(wfBin(1)),edge < 0]);
stops = samples/Fs;
% Identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
sampleVec = 1:length(logical(wfBin));
highSamples = unique([1,sampleVec(wfBin),sampleVec(end)]); 
lowSamples = unique([1,sampleVec(not(wfBin)),sampleVec(end)]); 
% Calculate the number of samples between consecutive high/low samples.
dHigh = diff(highSamples);
dLow = diff(lowSamples);
% Identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
restLength = dHigh(dHigh > 1);
restDur = restLength/Fs;
whiskLength = dLow(dLow > 1);
whiskDur = whiskLength/Fs;
% Control for the beginning/end of the trial to correctly map rests/whisks
% onto the whisk_starts. Use index 2 and end-1 since it is assumed that the
% first and last indexes of a trial are the end/beginning of a volitional movement.
if not(wfBin(2)) 
    whiskDur = [NaN,whiskDur];
end
if wfBin(end - 1)
    whiskDur(end) = [];
end
% Calculate the time to the closest puff
% If no puff occurred during the trial, store 0 as a place holder.
if isempty(puffTimes)
    puffTimes = 0;
end
puffMat = ones(length(samples),1)*puffTimes;
restMat = samples'*ones(1,length(puffTimes))/Fs;
puffTimeElapsed = abs(restMat - puffMat);
% Convert to cell
puffTimeCell = mat2cell(puffTimeElapsed,ones(length(samples),1));
% Compile into a structure
Rest.eventTime = stops';
Rest.duration = restDur';
Rest.puffDistance = puffTimeCell;
Rest.whiskDuration = whiskDur';
end

