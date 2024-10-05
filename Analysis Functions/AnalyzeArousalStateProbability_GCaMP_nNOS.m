function [Results_ArousalStateProb_GCaMP] = AnalyzeArousalStateProbability_GCaMP(animalID,group,set,rootFolder,delim,Results_ArousalStateProb_GCaMP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% find and load RestData struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID)
% scoring results
modelScoringResults = [animalID '_Forest_ScoringResults.mat'];
load(modelScoringResults)
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {0};
RestStimCriteria.Fieldname = {'stimDistances'};
RestStimCriteria.Comparison = {'gt'};
RestStimCriteria.Value = {5};
% pull data from RestData structure
[restLogical] = FilterEvents_IOS(RestData.HbT.LH,RestCriteria);
[stimLogical] = FilterEvents_IOS(RestData.HbT.LH,RestStimCriteria);
combRestLogical = logical(restLogical.*stimLogical);
restFileIDs = RestData.HbT.LH.fileIDs(combRestLogical,:);
restEventTimes = RestData.HbT.LH.eventTimes(combRestLogical,:);
restDurations = RestData.HbT.LH.durations(combRestLogical,:);
% 5 second events
bins = {'five','ten','fifteen','twenty','twentyfive','thirty','thirtyfive','forty','fortyfive','fifty','fiftyfive','sixty','sixtyplus'};
for aa = 1:length(bins)
    data.(bins{1,aa}) = [];
end
% group data based on resting duration
a5 = 1; a10 = 1; a15 = 1; a20 = 1; a25 = 1; a30 = 1; a35 = 1; a40 = 1; a45 = 1; a50 = 1; a55 = 1; a60 = 1; a61 = 1;
for aa = 1:length(restDurations)
    duration = restDurations(aa,1);
    if duration <= 5
        data.five.fileIDs{a5,1} = restFileIDs{aa,1};
        data.five.durations(a5,1) = restDurations(aa,1);
        data.five.eventTimes(a5,1) = restEventTimes(aa,1);
        a5 = a5 + 1;
    elseif duration > 5 && duration <= 10
        data.ten.fileIDs{a10,1} = restFileIDs{aa,1};
        data.ten.durations(a10,1) = restDurations(aa,1);
        data.ten.eventTimes(a10,1) = restEventTimes(aa,1);
        a10 = a10 + 1;
    elseif duration > 10 && duration <= 15
        data.fifteen.fileIDs{a15,1} = restFileIDs{aa,1};
        data.fifteen.durations(a15,1) = restDurations(aa,1);
        data.fifteen.eventTimes(a15,1) = restEventTimes(aa,1);
        a15 = a15 + 1;
    elseif duration > 15 && duration <= 20
        data.twenty.fileIDs{a20,1} = restFileIDs{aa,1};
        data.twenty.durations(a20,1) = restDurations(aa,1);
        data.twenty.eventTimes(a20,1) = restEventTimes(aa,1);
        a20 = a20 + 1;
    elseif duration > 20 && duration <= 25
        data.twentyfive.fileIDs{a25,1} = restFileIDs{aa,1};
        data.twentyfive.durations(a25,1) = restDurations(aa,1);
        data.twentyfive.eventTimes(a25,1) = restEventTimes(aa,1);
        a25 = a25 + 1;
    elseif duration > 25 && duration <= 30
        data.thirty.fileIDs{a30,1} = restFileIDs{aa,1};
        data.thirty.durations(a30,1) = restDurations(aa,1);
        data.thirty.eventTimes(a30,1) = restEventTimes(aa,1);
        a30 = a30 + 1;
    elseif duration > 30 && duration <= 35
        data.thirtyfive.fileIDs{a35,1} = restFileIDs{aa,1};
        data.thirtyfive.durations(a35,1) = restDurations(aa,1);
        data.thirtyfive.eventTimes(a35,1) = restEventTimes(aa,1);
        a35 = a35 + 1;
    elseif duration > 35 && duration <= 40
        data.forty.fileIDs{a40,1} = restFileIDs{aa,1};
        data.forty.durations(a40,1) = restDurations(aa,1);
        data.forty.eventTimes(a40,1) = restEventTimes(aa,1);
        a40 = a40 + 1;
    elseif duration > 40 && duration <= 45
        data.fortyfive.fileIDs{a45,1} = restFileIDs{aa,1};
        data.fortyfive.durations(a45,1) = restDurations(aa,1);
        data.fortyfive.eventTimes(a45,1) = restEventTimes(aa,1);
        a45 = a45 + 1;
    elseif duration > 45 && duration <= 50
        data.fifty.fileIDs{a50,1} = restFileIDs{aa,1};
        data.fifty.durations(a50,1) = restDurations(aa,1);
        data.fifty.eventTimes(a50,1) = restEventTimes(aa,1);
        a50 = a50 + 1;
    elseif duration > 50 && duration <= 55
        data.fiftyfive.fileIDs{a55,1} = restFileIDs{aa,1};
        data.fiftyfive.durations(a55,1) = restDurations(aa,1);
        data.fiftyfive.eventTimes(a55,1) = restEventTimes(aa,1);
        a55 = a55 + 1;
    elseif duration > 55 && duration <= 60
        data.sixty.fileIDs{a60,1} = restFileIDs{aa,1};
        data.sixty.durations(a60,1) = restDurations(aa,1);
        data.sixty.eventTimes(a60,1) = restEventTimes(aa,1);
        a60 = a60 + 1;
    elseif duration > 60
        data.sixtyplus.fileIDs{a61,1} = restFileIDs{aa,1};
        data.sixtyplus.durations(a61,1) = restDurations(aa,1);
        data.sixtyplus.eventTimes(a61,1) = restEventTimes(aa,1);
        a61 = a61 + 1;
    end
end
% go through each event and determine the probabilty of awake
for bb = 1:length(bins)
    bin = bins{1,bb};
    try
        for cc = 1:length(data.(bin).fileIDs)
            fileID = data.(bin).fileIDs{cc,1};
            eventTime = floor(data.(bin).eventTimes(cc,1));
            duration = floor(data.(bin).durations(cc,1));
            indFileScores = [];
            for ee = 1:length(ScoringResults.fileIDs)
                if strcmp(fileID,ScoringResults.fileIDs{ee,1})
                    indFileScores = ScoringResults.labels{ee,1};
                end
            end
            startRemTime = rem(eventTime,5);
            binStartTime = eventTime - startRemTime;
            binStartIndex = (binStartTime/5) + 1;
            binEndIndex = floor(duration/5);
            indexScores = indFileScores(binStartIndex:binStartIndex + binEndIndex);
            indexLogical = strcmp(indexScores,'Not Sleep');
            if sum(indexLogical) == length(indexLogical)
                data.(bin).awakeLogical(cc,1) = 1;
            else
                data.(bin).awakeLogical(cc,1) = 0;
            end
        end
    catch
        data.(bin).awakeLogical = 0;
    end
end
% save results
for ff = 1:length(bins)
    Results_ArousalStateProb_GCaMP.(group).(animalID).SleepProbability.(bins{1,ff}).awakeLogical = data.(bins{1,ff}).awakeLogical;
end
% analyze trial hypogram and awake probability based on trial time
% identify the unique file IDs, unique imaging days, and scoring labels from the file list
allScoringLabels = ScoringResults.alllabels;
allFileIDs = ScoringResults.allfileIDs;
uniqueFileIDs = unique(allFileIDs);
for aa = 1:size(uniqueFileIDs,1)
    [animalID,allFileDates{aa,1},~] = GetFileInfo_IOS(uniqueFileIDs{aa,1});
    allFileDays{aa,1} = ConvertDate_IOS(allFileDates{aa,1});
end
data.uniqueDates = unique(allFileDates);
data.uniqueDays = unique(allFileDays);
% determine how many 5 second bins there are for each individual day
for bb = 1:size(data.uniqueDays,1)
    data.dayBinLengths{bb,1} = find(~cellfun('isempty',strfind(allFileIDs,data.uniqueDates{bb,1})));
end
% extract the file IDs and scoring labels that correspond to each day's indeces
for cc = 1:size(data.dayBinLengths,1)
    dayInds = data.dayBinLengths{cc,1};
    data.dayScoreLabels{cc,1} = allScoringLabels(dayInds);
    data.dayScoreFileIDs{cc,1} = allFileIDs(dayInds);
end
trialDuration = 15; % minutes
binTime = 5; % seconds
fileBinLength = (trialDuration*60)/binTime;
% further break down each day's scores into the scores for each individual file
for dd = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{dd,1};
    uniqueDayFileIDs = unique(data.dayScoreFileIDs{dd,1});
    for ee = 1:size(uniqueDayFileIDs,1)
        if ee == 1
            data.(uniqueDay).indFileData{ee,1} = data.dayScoreLabels{dd,1}(1:fileBinLength);
        else
            data.(uniqueDay).indFileData{ee,1} = data.dayScoreLabels{dd,1}((ee - 1)*fileBinLength + 1:ee*fileBinLength);
        end
    end
end
% calculate the time difference between every file to append padding 'Time Pad' to the end of the leading file's score labels
for ff = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{ff,1};
    uniqueDayFileIDs = unique(data.dayScoreFileIDs{ff,1});
    % start with file 2 to focus on the differences between each file
    for gg = 2:size(uniqueDayFileIDs,1)
        leadFileID = uniqueDayFileIDs{gg - 1,1};
        lagFileID = uniqueDayFileIDs{gg,1};
        [~,~,leadFileInfo] = GetFileInfo_IOS(leadFileID);
        [~,~,lagFileInfo] = GetFileInfo_IOS(lagFileID);
        leadFileStr = ConvertDateTime_IOS(leadFileInfo);
        lagFileStr = ConvertDateTime_IOS(lagFileInfo);
        leadFileTime = datevec(leadFileStr);
        lagFileTime = datevec(lagFileStr);
        timeDifference = etime(lagFileTime,leadFileTime) - (trialDuration*60); % seconds
        timePadBins.(uniqueDay){gg - 1,1} = cell(floor(timeDifference/binTime),1);
        timePadBins.(uniqueDay){gg - 1,1}(:) = {'Time Pad'};
    end
end
% apply the time padding to the end of each file
for hh = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{hh,1};
    for jj = 1:size(data.(uniqueDay).indFileData,1) - 1
        data.(uniqueDay).indFileData{jj,1} = vertcat(data.(uniqueDay).indFileData{jj,1},timePadBins.(uniqueDay){jj,1});
    end
end
% concatendate the data for each day now that the padding is added at the end of each file
for kk = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{kk,1};
    data.(uniqueDay).catData = [];
    for mm = 1:size(data.(uniqueDay).indFileData,1)
        data.(uniqueDay).catData = vertcat(data.(uniqueDay).catData,data.(uniqueDay).indFileData{mm,1});
    end
end
% prepare indeces for each behavior
for nn = 1:size(data.uniqueDays,1)
    uniqueDay = data.uniqueDays{nn,1};
    data.(uniqueDay).NotSleep_inds = NaN(1,size(data.(uniqueDay).catData,1));
    data.(uniqueDay).NREM_inds = NaN(1,size(data.(uniqueDay).catData,1));
    data.(uniqueDay).REM_inds = NaN(1,size(data.(uniqueDay).catData,1));
    data.(uniqueDay).TimePad_inds = NaN(1,size(data.(uniqueDay).catData,1));
    for oo = 1:size(data.(uniqueDay).catData,1)
        if strcmp(data.(uniqueDay).catData{oo,1},'Not Sleep') == true
            data.(uniqueDay).NotSleep_inds(1,oo) = 1;
            data.(uniqueDay).NREM_inds(1,oo) = NaN;
            data.(uniqueDay).REM_inds(1,oo) = NaN;
            data.(uniqueDay).TimePad_inds(1,oo) = NaN;
            data.(uniqueDay).AwakeProb_inds(1,oo) = 1;
            data.(uniqueDay).NREMProb_inds(1,oo) = 0;
            data.(uniqueDay).REMProb_inds(1,oo) = 0;
        elseif strcmp(data.(uniqueDay).catData{oo,1},'NREM Sleep') == true
            data.(uniqueDay).NotSleep_inds(1,oo) = NaN;
            data.(uniqueDay).NREM_inds(1,oo) = 1;
            data.(uniqueDay).REM_inds(1,oo) = NaN;
            data.(uniqueDay).TimePad_inds(1,oo) = NaN;
            data.(uniqueDay).AwakeProb_inds(1,oo) = 0;
            data.(uniqueDay).NREMProb_inds(1,oo) = 1;
            data.(uniqueDay).REMProb_inds(1,oo) = 0;
        elseif strcmp(data.(uniqueDay).catData{oo,1},'REM Sleep') == true
            data.(uniqueDay).NotSleep_inds(1,oo) = NaN;
            data.(uniqueDay).NREM_inds(1,oo) = NaN;
            data.(uniqueDay).REM_inds(1,oo) = 1;
            data.(uniqueDay).TimePad_inds(1,oo) = NaN;
            data.(uniqueDay).AwakeProb_inds(1,oo) = 0;
            data.(uniqueDay).NREMProb_inds(1,oo) = 0;
            data.(uniqueDay).REMProb_inds(1,oo) = 1;
        elseif strcmp(data.(uniqueDay).catData{oo,1},'Time Pad') == true
            data.(uniqueDay).NotSleep_inds(1,oo) = NaN;
            data.(uniqueDay).NREM_inds(1,oo) = NaN;
            data.(uniqueDay).REM_inds(1,oo) = NaN;
            data.(uniqueDay).TimePad_inds(1,oo) = 1;
            data.(uniqueDay).AwakeProb_inds(1,oo) = NaN;
            data.(uniqueDay).NREMProb_inds(1,oo) = NaN;
            data.(uniqueDay).REMProb_inds(1,oo) = NaN;
        end
    end
end
% save results
for aa = 1:size(data.uniqueDays,1)
    Results_ArousalStateProb_GCaMP.(group).(animalID).Hypnogram.(data.uniqueDays{aa,1}) = data.(data.uniqueDays{aa,1});
end
% percentage of time in each state
Results_ArousalStateProb_GCaMP.(group).(animalID).numberOfScores = length(ScoringResults.alllabels);
Results_ArousalStateProb_GCaMP.(group).(animalID).awakePercent = round((sum(strcmp(ScoringResults.alllabels,'Not Sleep'))/length(ScoringResults.alllabels))*100,1);
Results_ArousalStateProb_GCaMP.(group).(animalID).nremPercent = round((sum(strcmp(ScoringResults.alllabels,'NREM Sleep'))/length(ScoringResults.alllabels))*100,1);
Results_ArousalStateProb_GCaMP.(group).(animalID).remPercent = round((sum(strcmp(ScoringResults.alllabels,'REM Sleep'))/length(ScoringResults.alllabels))*100,1);
% save results
cd([rootFolder delim 'Results_Turner'])
save('Results_ArousalStateProb_GCaMP.mat','Results_ArousalStateProb_GCaMP')
cd([rootFolder delim 'Data'])