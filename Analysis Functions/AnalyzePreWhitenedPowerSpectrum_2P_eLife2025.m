function [Results_PreWhitenedPowerSpec_2P_eLife2025] = AnalyzePreWhitenedPowerSpectrum_2P_eLife2025(animalID,group,set,rootFolder,delim,Results_PreWhitenedPowerSpec_2P_eLife2025)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% character list of MergedData files
mergedDirectory = dir('*_MergedData.mat');
mergedDataFiles = {mergedDirectory.name}';
mergedDataFileIDs = char(mergedDataFiles);
% find and load ManualDecisions struct
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID)
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% lowpass filter
samplingRate = 5;
% All
allData = [];
for aa = 1:size(mergedDataFileIDs,1)
    mergedDataFileID = mergedDataFileIDs(aa,:);
    [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P_eLife2025(mergedDataFileID);
    if strcmp(vesselID(1),'V') == false
        load(mergedDataFileID,'-mat')
        strDay = ConvertDate_2P_eLife2025(fileDate);
        if isfield(allData,vesselID) == false
            allData.(vesselID) = [];
            binWhisking.(vesselID) = [];
        end
        % filter and detrend data
        vesselDiam = MergedData.data.vesselDiameter.data;
        normVesselDiam = diff(detrend((vesselDiam - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay))./ RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay),'constant'),1);
        allData.(vesselID) = horzcat(allData.(vesselID),normVesselDiam');
        binWhisk = MergedData.data.binWhiskerAngle;
        [linkedBinarizedWhiskers] = LinkBinaryEvents_2P_eLife2025(gt(binWhisk,0),[round(30/3),0]);
        binWhiskingPercent = sum(linkedBinarizedWhiskers)/length(linkedBinarizedWhiskers)*100 ;
        binWhisking.(vesselID) = horzcat(binWhisking.(vesselID),binWhiskingPercent);
    end
end
% parameters for mtspectrumc - information available in function
params.tapers = [5,9]; % Tapers [n, 2n - 1]
params.pad = 1;
params.Fs = samplingRate;
params.fpass = [0,0.5];  % Pass band [0, nyquist]
params.trialave = 1;
params.err = [2,0.05];
if isempty(allData) == false
    allDataVesselIDs = fieldnames(allData);
    for bb = 1:length(allDataVesselIDs)
        allDataVID = allDataVesselIDs{bb,1};
        [allData_S{bb,1},allData_f{bb,1},allData_sErr{bb,1}] = mtspectrumc(allData.(allDataVID),params);
        % save results
        Results_PreWhitenedPowerSpec_2P.(group).(animalID).(allDataVID).whiskingPerc = mean(binWhisking.(allDataVID));
        Results_PreWhitenedPowerSpec_2P.(group).(animalID).(allDataVID).S = allData_S{bb,1};
        Results_PreWhitenedPowerSpec_2P.(group).(animalID).(allDataVID).f = allData_f{bb,1};
        Results_PreWhitenedPowerSpec_2P.(group).(animalID).(allDataVID).sErr = allData_sErr{bb,1};
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PreWhitenedPowerSpec_2P_eLife2025.mat','Results_PreWhitenedPowerSpec_2P')
cd([rootFolder delim 'Data'])