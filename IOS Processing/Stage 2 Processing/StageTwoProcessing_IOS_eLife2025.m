%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
zap; 
% process RawData, set behavior thresholds and create ProcData
[procDataFileIDs] = ProcessRawDataFiles_IOS_eLife2025();
% set ROIs for IOS camera and collect reflectance from each region of interest
ProcessIntrinsicData_IOS_eLife2025(procDataFileIDs)
% analyze the spectrogram for each session
CreateTrialSpectrograms_IOS_eLife2025();
% add Heart Rate for valid imaging wavelengths
ExtractHeartRate_IOS_eLife2025(procDataFileIDs)
% categorize data
[procDataFileIDs] = CategorizeData_IOS_eLife2025();
% create RestData data structure
[RestData] = ExtractRestingData_IOS_eLife2025(procDataFileIDs,1);
% create Baselines data structure
[RestingBaselines] = CalculateRestingBaselines_IOS_eLife2025(RestData,60);
% create Baselines for spectrogram data
[RestingBaselines] = CalculateSpectrogramBaselines_IOS_eLife2025(RestingBaselines,'setDuration');
% normalize spectrogram by baseline
NormalizeSpectrograms_IOS_eLife2025(RestingBaselines);