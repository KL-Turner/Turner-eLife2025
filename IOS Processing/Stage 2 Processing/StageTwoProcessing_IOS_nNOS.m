%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
zap; 
% process RawData, set behavior thresholds and create ProcData
[procDataFileIDs] = ProcessRawDataFiles_IOS_nNOS();
% set ROIs for IOS camera and collect reflectance from each region of interest
ProcessIntrinsicData_IOS_nNOS(procDataFileIDs)
% analyze the spectrogram for each session
CreateTrialSpectrograms_IOS_nNOS();
% add Heart Rate for valid imaging wavelengths
ExtractHeartRate_IOS_nNOS(procDataFileIDs)
% categorize data
[procDataFileIDs] = CategorizeData_IOS_nNOS();
% create RestData data structure
[RestData] = ExtractRestingData_IOS_nNOS(procDataFileIDs,1);
% create Baselines data structure
[RestingBaselines] = CalculateRestingBaselines_IOS_nNOS(RestData,60);
% create Baselines for spectrogram data
[RestingBaselines] = CalculateSpectrogramBaselines_IOS_nNOS(RestingBaselines,'setDuration');
% normalize spectrogram by baseline
NormalizeSpectrograms_IOS_nNOS(RestingBaselines);