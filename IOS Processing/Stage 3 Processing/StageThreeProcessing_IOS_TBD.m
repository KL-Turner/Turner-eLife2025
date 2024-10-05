%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
zap;
% manually select files for custom baseline calculation
[RestingBaselines,ManualDecisions,procDataFileIDs,animalID] = CalculateManualRestingBaselinesTimeIndeces_IOS('reflectance');
% pixel-wise resting baselines
[RestingBaselines] = CalculatePixelWiselRestingBaselines_IOS(procDataFileIDs,RestingBaselines);
% IOS analysis
CalculateFullSpectroscopy_IOS(procDataFileIDs,RestingBaselines)
% re-create the RestData structure now that HbT (and/or corrected GCaMP) is available
[RestData] = ExtractRestingData_IOS(procDataFileIDs,3);
% create the EventData structure for CBV and neural data
[EventData] = ExtractEventTriggeredData_IOS(procDataFileIDs);
% normalize RestData structures by the resting baseline
[RestData] = NormRestDataStruct_IOS(animalID,RestData,RestingBaselines,'manualSelection');
% normalize EventData structures by the resting baseline
[EventData] = NormEventDataStruct_IOS(animalID,EventData,RestingBaselines,'manualSelection');
% find spectrogram baselines for each day
[RestingBaselines] = CalculateSpectrogramBaselines_IOS(RestingBaselines,'manualSelection');
% normalize spectrogram by baseline
NormalizeSpectrograms_IOS(RestingBaselines);
% create a structure with all spectrograms for convenient analysis further downstream
CreateAllSpecDataStruct_IOS()
% generate single trial figures
GenerateTrialFigures_IOS(procDataFileIDs,RestingBaselines);