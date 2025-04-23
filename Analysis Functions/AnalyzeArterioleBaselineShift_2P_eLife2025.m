function [Results_BaselineShift_2P] = AnalyzeArterioleBaselineShift_2P_eLife2025(animalID,group,set,rootFolder,delim,Results_BaselineShift_2P)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
baseLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging/Combined Imaging'];
cd(baseLocation)
% find and load RestingBaselines strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% cd to data location
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging/Isoflurane Trials'];
cd(dataLocation)
% character list of all MergedData files
mergedDirectory = dir('*_MergedData.mat');
mergedDataFiles = {mergedDirectory.name}';
mergedDataFileIDs = char(mergedDataFiles);
samplingRate = 5;
for aa = 1:size(mergedDataFileIDs,1)
    mergedDataFileID = mergedDataFileIDs(aa,:);
    load(mergedDataFileID)
    [~,~,fileDate,~,~,vID] = GetFileInfo2_2P_eLife2025(mergedDataFileID);
    strDay = ConvertDate_2P_eLife2025(fileDate);
    % save results
    Results_BaselineShift_2P.(group).(animalID).(vID).diameter = mean(MergedData.data.vesselDiameter.data(1:30*samplingRate));
    Results_BaselineShift_2P.(group).(animalID).(vID).baseline = RestingBaselines.manualSelection.vesselDiameter.data.(vID).(strDay);
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_BaselineShift_2P.mat','Results_BaselineShift_2P')
cd([rootFolder delim 'Data'])