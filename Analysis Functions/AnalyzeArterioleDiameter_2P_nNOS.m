function [Results_Diameter_2P] = AnalyzeArterioleDiameter_2P(animalID,group,set,rootFolder,delim,Results_Diameter_2P)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID)
% find and load ManualDecisions struct
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID)
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)
% find and load RestingBaselines struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% parameters 
params.Offset = 2;
params.minTime.Rest = 10;
params.minTime.Whisk = params.Offset + 5;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% criteria for whisking
WhiskCriteria.Fieldname = {'duration','duration'};
WhiskCriteria.Comparison = {'gt','lt'};
WhiskCriteria.Value = {2,5};
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
% lowpass filter
samplingRate = RestData.vesselDiameter.data.samplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
%% Rest
[restLogical] = FilterEvents_2P(RestData.vesselDiameter.data,RestCriteria);
combRestLogical = logical(restLogical);
restVesselData = RestData.vesselDiameter.data.data(combRestLogical,:);
restFileIDs = RestData.vesselDiameter.data.fileIDs(combRestLogical,:);
restVesselIDs = RestData.vesselDiameter.data.vesselIDs(combRestLogical,:);
restDurations = RestData.vesselDiameter.data.durations(combRestLogical,:);
restEventTimes = RestData.vesselDiameter.data.eventTimes(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[finalRestVesselData,finalRestFileIDs,finalRestVesselIDs,~,~] = RemoveInvalidData_2P(restVesselData,restFileIDs,restVesselIDs,restDurations,restEventTimes,ManualDecisions);
% go through the data and normalize + filter each rest epoch based on individual vessels
uniqueRestVesselIDs = unique(finalRestVesselIDs);
for aa = 1:length(uniqueRestVesselIDs)
    cc = 1;
    for bb = 1:length(finalRestVesselIDs)
        if strcmp(uniqueRestVesselIDs{aa,1},finalRestVesselIDs{bb,1})
            strDay = ConvertDate_2P(finalRestFileIDs{bb,1}(1:6));
            tempRestData.(uniqueRestVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,((finalRestVesselData{bb,1} - RestingBaselines.manualSelection.vesselDiameter.data.(uniqueRestVesselIDs{aa,1}).(strDay))/RestingBaselines.manualSelection.vesselDiameter.data.(uniqueRestVesselIDs{aa,1}).(strDay)));
            tempRestFileIDs.(uniqueRestVesselIDs{aa,1}){cc,1} = finalRestFileIDs{bb,1};
            cc = cc + 1;
        end
    end
end
% take the average of each vessel's individual resting event
for dd = 1:length(uniqueRestVesselIDs)
    for ee = 1:length(tempRestData.(uniqueRestVesselIDs{dd,1}))
        tempRestDataMeans.(uniqueRestVesselIDs{dd,1})(ee,1) = mean(tempRestData.(uniqueRestVesselIDs{dd,1}){ee,1})*100;
        tempRestDataMaxs.(uniqueRestVesselIDs{dd,1})(ee,1) = max(tempRestData.(uniqueRestVesselIDs{dd,1}){ee,1})*100;
        tempRestDataInd.(uniqueRestVesselIDs{dd,1}){ee,1} = tempRestData.(uniqueRestVesselIDs{dd,1}){ee,1}*100;
    end
end
% take the average of each vessel's total resting events
for ff = 1:length(uniqueRestVesselIDs)
    % save results
    Results_Diameter_2P.(group).(animalID).(uniqueRestVesselIDs{ff,1}).Rest.indEvents = tempRestDataInd.(uniqueRestVesselIDs{ff,1});
    Results_Diameter_2P.(group).(animalID).(uniqueRestVesselIDs{ff,1}).Rest.mean = tempRestDataMeans.(uniqueRestVesselIDs{ff,1});
    Results_Diameter_2P.(group).(animalID).(uniqueRestVesselIDs{ff,1}).Rest.max = tempRestDataMaxs.(uniqueRestVesselIDs{ff,1});
    Results_Diameter_2P.(group).(animalID).(uniqueRestVesselIDs{ff,1}).Rest.fileIDs = tempRestFileIDs.(uniqueRestVesselIDs{ff,1});
end
%% Whisk
[whiskLogical] = FilterEvents_2P(EventData.vesselDiameter.data.whisk,WhiskCriteria);
combWhiskLogical = logical(whiskLogical);
whiskVesselData = EventData.vesselDiameter.data.whisk.data(combWhiskLogical,:);
whiskFileIDs = EventData.vesselDiameter.data.whisk.fileIDs(combWhiskLogical,:);
whiskVesselIDs = EventData.vesselDiameter.data.whisk.vesselIDs(combWhiskLogical,:);
whiskDurations = EventData.vesselDiameter.data.whisk.duration(combWhiskLogical,:);
whiskEventTimes = EventData.vesselDiameter.data.whisk.eventTime(combWhiskLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[finalWhiskVesselData,finalWhiskFileIDs,finalWhiskVesselIDs,~,~] = RemoveInvalidData_2P(whiskVesselData,whiskFileIDs,whiskVesselIDs,whiskDurations,whiskEventTimes,ManualDecisions);
% go through the data and normalize + filter each whisk event based on individual vessels
uniqueWhiskVesselIDs = unique(finalWhiskVesselIDs);
for aa = 1:length(uniqueWhiskVesselIDs)
    cc = 1;
    for bb = 1:length(finalWhiskVesselIDs)
        if strcmp(uniqueWhiskVesselIDs{aa,1},finalWhiskVesselIDs{bb,1})
            strDay = ConvertDate_2P(finalWhiskFileIDs{bb,1}(1:6));
            tempWhiskData.(uniqueWhiskVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,((finalWhiskVesselData(bb,:) - RestingBaselines.manualSelection.vesselDiameter.data.(uniqueWhiskVesselIDs{aa,1}).(strDay))/RestingBaselines.manualSelection.vesselDiameter.data.(uniqueWhiskVesselIDs{aa,1}).(strDay)));
            tempWhiskFileIDs.(uniqueWhiskVesselIDs{aa,1}){cc,1} = finalWhiskFileIDs{bb,1};
            cc = cc + 1;
        end
    end
end
% mean-subtract 2 seconds prior to whisk
for dd = 1:length(uniqueWhiskVesselIDs)
    for ee = 1:length(tempWhiskData.(uniqueWhiskVesselIDs{dd,1}))
        tempWhiskDataB.(uniqueWhiskVesselIDs{dd,1}){ee,1} = tempWhiskData.(uniqueWhiskVesselIDs{dd,1}){ee,1} - mean(tempWhiskData.(uniqueWhiskVesselIDs{dd,1}){ee,1}(1:params.Offset*samplingRate));
    end
end
% take the average of each vessel's individual whisking event from onset through 5 seconds
for dd = 1:length(uniqueWhiskVesselIDs)
    for ee = 1:length(tempWhiskData.(uniqueWhiskVesselIDs{dd,1}))
        tempWhiskDataMeans.(uniqueWhiskVesselIDs{dd,1})(ee,1) = mean(tempWhiskDataB.(uniqueWhiskVesselIDs{dd,1}){ee,1}(params.Offset*samplingRate:params.minTime.Whisk*samplingRate))*100;
        tempWhiskDataMaxs.(uniqueWhiskVesselIDs{dd,1})(ee,1) = max(tempWhiskDataB.(uniqueWhiskVesselIDs{dd,1}){ee,1}(params.Offset*samplingRate:params.minTime.Whisk*samplingRate))*100;
        tempWhiskDataInd.(uniqueWhiskVesselIDs{dd,1}){ee,1} = tempWhiskDataB.(uniqueWhiskVesselIDs{dd,1}){ee,1}(params.Offset*samplingRate:params.minTime.Whisk*samplingRate)*100;
    end
end
% take the average of each vessel's total whisking events
for ff = 1:length(uniqueWhiskVesselIDs)
    % save results
    Results_Diameter_2P.(group).(animalID).(uniqueWhiskVesselIDs{ff,1}).Whisk.indEvents = tempWhiskDataInd.(uniqueWhiskVesselIDs{ff,1});
    Results_Diameter_2P.(group).(animalID).(uniqueWhiskVesselIDs{ff,1}).Whisk.mean = tempWhiskDataMeans.(uniqueWhiskVesselIDs{ff,1});
    Results_Diameter_2P.(group).(animalID).(uniqueWhiskVesselIDs{ff,1}).Whisk.max = tempWhiskDataMaxs.(uniqueWhiskVesselIDs{ff,1});
    Results_Diameter_2P.(group).(animalID).(uniqueWhiskVesselIDs{ff,1}).Whisk.fileIDs = tempWhiskFileIDs.(uniqueWhiskVesselIDs{ff,1});
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_Diameter_2P.mat','Results_Diameter_2P')
cd([rootFolder delim 'Data'])