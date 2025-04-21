%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
zap;
% manually select files for custom baseline calculation
[RestingBaselines,ManualDecisions,procDataFileIDs,animalID] = CalculateManualRestingBaselinesTimeIndeces_IOS_eLife2025('reflectance');
% pixel-wise resting baselines
[RestingBaselines] = CalculatePixelWiselRestingBaselines_IOS_eLife2025(procDataFileIDs,RestingBaselines);
% IOS analysis
CalculateFullSpectroscopy_IOS_eLife2025(procDataFileIDs,RestingBaselines)
% re-create the RestData structure now that HbT (and/or corrected GCaMP) is available
[RestData] = ExtractRestingData_IOS_eLife2025(procDataFileIDs,3);
% create the EventData structure for CBV and neural data
[EventData] = ExtractEventTriggeredData_IOS_eLife2025(procDataFileIDs);
% normalize RestData structures by the resting baseline
[RestData] = NormRestDataStruct_IOS_eLife2025(animalID,RestData,RestingBaselines,'manualSelection');
% normalize EventData structures by the resting baseline
[EventData] = NormEventDataStruct_IOS_eLife2025(animalID,EventData,RestingBaselines,'manualSelection');
% find spectrogram baselines for each day
[RestingBaselines] = CalculateSpectrogramBaselines_IOS_eLife2025(RestingBaselines,'manualSelection');
% normalize spectrogram by baseline
NormalizeSpectrograms_IOS_eLife2025(RestingBaselines);
% create a structure with all spectrograms for convenient analysis further downstream
CreateAllSpecDataStruct_IOS_eLife2025()
% generate single trial figures
GenerateTrialFigures_IOS_eLife2025(procDataFileIDs,RestingBaselines);