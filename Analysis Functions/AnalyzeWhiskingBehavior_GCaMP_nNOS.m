function [Results_WhiskBehav_GCaMP] = AnalyzeWhiskingBehavior_GCaMP(animalID,group,set,rootFolder,delim,Results_WhiskBehav_GCaMP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
% analyze data during periods of whisking
whiskDurations = EventData.HbT.LH.whisk.duration;
fileIDs = unique(EventData.HbT.LH.whisk.fileIDs);
imagingDuration = 15*length(fileIDs);
% save results
Results_WhiskBehav_GCaMP.(group).(animalID).whiskDurations = whiskDurations;
Results_WhiskBehav_GCaMP.(group).(animalID).whiskDurationSec = sum(whiskDurations);
Results_WhiskBehav_GCaMP.(group).(animalID).whiskDurationPerc = ((sum(whiskDurations)/60)/imagingDuration)*100;
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_WhiskBehav_GCaMP.mat','Results_WhiskBehav_GCaMP')
cd([rootFolder delim 'Data'])