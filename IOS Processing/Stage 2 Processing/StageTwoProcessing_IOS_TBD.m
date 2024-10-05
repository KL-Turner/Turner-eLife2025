%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
zap; 
% process RawData, set behavior thresholds and create ProcData
[procDataFileIDs] = ProcessRawDataFiles_IOS();
% set ROIs for IOS camera and collect reflectance from each region of interest
ProcessIntrinsicData_IOS(procDataFileIDs)
% analyze the spectrogram for each session
CreateTrialSpectrograms_IOS();
% add Heart Rate for valid imaging wavelengths
ExtractHeartRate_IOS(procDataFileIDs)
% categorize data
[procDataFileIDs] = CategorizeData_IOS();
% create RestData data structure
[RestData] = ExtractRestingData_IOS(procDataFileIDs,1);
% create Baselines data structure
[RestingBaselines] = CalculateRestingBaselines_IOS(RestData,60);
% create Baselines for spectrogram data
[RestingBaselines] = CalculateSpectrogramBaselines_IOS(RestingBaselines,'setDuration');
% normalize spectrogram by baseline
NormalizeSpectrograms_IOS(RestingBaselines);