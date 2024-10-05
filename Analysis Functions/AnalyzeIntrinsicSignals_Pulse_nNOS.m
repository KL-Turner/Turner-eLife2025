function [Results_IntSig_Pulse] = AnalyzeIntrinsicSignals_Pulse(animalID,group,set,rootFolder,delim,Results_IntSig_Pulse)
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
load(restDataFileID,'-mat')
% find and load ManualDecisions struct
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
params.minTime.Rest = 10;
params.Offset = 2;
% criteria for the resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestStimCriteria.Fieldname = {'puffDistances'};
RestStimCriteria.Comparison = {'gt'};
RestStimCriteria.Value = {5};
% loop variables
% lowpass filter
samplingRate = RestData.CBV_HbT.adjBarrels.CBVCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% Rest
[restLogical] = FilterEvents_IOS(RestData.CBV_HbT.adjBarrels,RestCriteria);
[stimLogical] = FilterEvents_IOS(RestData.CBV_HbT.adjBarrels,RestStimCriteria);
combRestLogical = logical(restLogical.*stimLogical);
restFileIDs = RestData.CBV_HbT.adjBarrels.fileIDs(combRestLogical,:);
restEventTimes = RestData.CBV_HbT.adjBarrels.eventTimes(combRestLogical,:);
restDurations = RestData.CBV_HbT.adjBarrels.durations(combRestLogical,:);
restData = RestData.CBV_HbT.adjBarrels.data(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
try
    [finalRestData,~,~,~] = RemoveInvalidData_IOS(restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
catch
    finalRestData = [];
end
% filter and average
procRestData = {}; restMean = []; varInd = [];
for gg = 1:length(finalRestData)
    procRestData{gg,1} = filtfilt(sos,g,finalRestData{gg,1});
    restMean(gg,1) = mean(procRestData{gg,1}(1:end));
    varInd(gg,1) = var(filtfilt(sos,g,finalRestData{gg,1}));
end
% save results
Results_IntSig_Pulse.(group).(animalID).Rest.indHbT = procRestData;
Results_IntSig_Pulse.(group).(animalID).Rest.HbT = restMean;
cd([rootFolder delim 'Results_Turner'])
save('Results_IntSig_Pulse.mat','Results_IntSig_Pulse')
cd([rootFolder delim 'Data'])