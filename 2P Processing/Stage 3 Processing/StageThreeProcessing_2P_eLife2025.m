%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: 1) Categorize behavioral (rest,whisk,stim) data using previously processed data structures, add 'flags'  
%            2) Create a temporary RestData structure that contains periods of rest - use this for initial figures
%            3) Analyze neural data and create different spectrograms for each file's electrodes
%            4) Uses periods when animal is not being stimulated or moving to establish an initial baseline
%            5) Manually select awake files for a slightly different baseline not based on hard time vals
%            6) Use the best baseline to convert reflectance changes to total hemoglobin
%            7) Re-create the RestData structure now that we can deltaHbT
%            8) Create an EventData structure looking at the different data types after whisking or stimulation
%            9) Apply the resting baseline to each data type to create a percentage change 
%            10) Use the time indeces of the resting baseline file to apply a percentage change to the spectrograms
%            11) Use the time indeces of the resting baseline file to create a reflectance pixel-based baseline
%            12) Generate a summary figure for all of the analyzed and processed data
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
% Character list of all MergedData files
mergedDirectory = dir('*_MergedData.mat');
mergedDataFiles = {mergedDirectory.name}';
mergedDataFileIDs = char(mergedDataFiles);
[animalID,~,~,~,~,~] = GetFileInfo2_2P_eLife2025(mergedDataFileIDs(1,:));
genSampleFigs = 'y';
stimulationType = input('Input stimulation type (single or pulse): ','s'); disp(' ')
dataTypes = {'vesselDiameter','corticalNeural','hippocampalNeural','EMG'};
neuralDataTypes = {'corticalNeural','hippocampalNeural'};
specNeuralDataTypes = {'rawCorticalNeural','rawHippocampalNeural'};
%% BLOCK PURPOSE: [1] Categorize data 
disp('Analyzing Block [1] Categorizing data.'); disp(' ')
for aa = 1:size(mergedDataFileIDs,1)
    mergedDataFileID = mergedDataFileIDs(aa,:);
    disp(['Analyzing file ' num2str(aa) ' of ' num2str(size(mergedDataFileIDs,1)) '...']); disp(' ')
    CategorizeData_2P_eLife2025(mergedDataFileID,stimulationType)
end
%% BLOCK PURPOSE: [2] Create RestData data structure.
disp('Analyzing Block [2] Creating RestData struct for vessels and neural data.'); disp(' ')
[RestData] = ExtractRestingData_2P_eLife2025(mergedDataFileIDs,dataTypes);
%% BLOCK PURPOSE: [3] Create EventData data structure.
disp('Analyzing Block [3] Creating EventData struct for vessels and neural data.'); disp(' ')
[EventData] = ExtractEventTriggeredData_2P_eLife2025(mergedDataFileIDs,dataTypes);
%% BLOCK PURPOSE: [4] Analyze the spectrogram for each session.
disp('Analyzing Block [4] Analyzing the spectrogram for each file and normalizing by the resting baseline.'); disp(' ')
CreateTrialSpectrograms_2P_eLife2025(mergedDataFileIDs,specNeuralDataTypes);
%% BLOCK PURPOSE: [5] Create Baselines data structure
disp('Analyzing Block [5] Create Baselines struct for CBV and neural data.'); disp(' ')
baselineType = 'setDuration';
trialDuration_sec = 900;
targetMinutes = 60;
[RestingBaselines] = CalculateRestingBaselines_2P_eLife2025(animalID,targetMinutes,trialDuration_sec,RestData);
% Find spectrogram baselines for each day
specDirectory = dir('*_SpecData.mat');
specDataFiles = {specDirectory.name}';
specDataFileIDs = char(specDataFiles);
[RestingBaselines] = CalculateSpectrogramBaselines_2P_eLife2025(animalID,neuralDataTypes,trialDuration_sec,specDataFileIDs,RestingBaselines,baselineType);
% Normalize spectrogram by baseline
NormalizeSpectrograms_2P_eLife2025(specDataFileIDs,neuralDataTypes,RestingBaselines);
%% BLOCK PURPOSE: [6] Generate first set of figures to remove unwanted data
disp('Analyzing Block [6] Generating sample figures for inspection.'); disp(' ')
if strcmp(genSampleFigs,'y') == true
    saveFigs = 'y';
    for bb = 1:size(mergedDataFileIDs,1)
        mergedDataFileID = mergedDataFileIDs(bb,:);
        disp(['Generating single trial figure: (' num2str(bb) '/' num2str(size(mergedDataFileIDs,1)) ')']); disp(' ')
        [figHandle] = GenerateSingleFigures_2P_eLife2025(mergedDataFileID,baselineType,saveFigs,RestingBaselines);
        close(figHandle)
    end
end
%% BLOCK PURPOSE: [7] Manually select files for custom baseline calculation
disp('Analyzing Block [7] Manually select files for custom baseline calculation.'); disp(' ')
[RestingBaselines] = CalculateManualRestingBaselinesTimeIndeces_2P_eLife2025;
%% BLOCK PURPOSE: [8] Analyze the spectrogram baseline for each session.
disp('Analyzing Block [8] Analyzing the spectrogram for each file and normalizing by the resting baseline.'); disp(' ')
updatedBaselineType = 'manualSelection';
% Find spectrogram baselines for each day
specDirectory = dir('*_SpecData.mat');
specDataFiles = {specDirectory.name}';
specDataFileIDs = char(specDataFiles);
[RestingBaselines] = CalculateSpectrogramBaselines_2P_eLife2025(animalID,neuralDataTypes,trialDuration_sec,specDataFileIDs,RestingBaselines,updatedBaselineType);
% Normalize spectrogram by baseline
NormalizeSpectrograms_2P_eLife2025(specDataFileIDs,neuralDataTypes,RestingBaselines);
% Create a structure with all spectrograms for convenient analysis further downstream
CreateAllSpecDataStruct_2P_eLife2025(animalID,neuralDataTypes)
%% BLOCK PURPOSE [9] Generate single trial figures
disp('Analyzing Block [9] Generating single trial summary figures'); disp(' ')
updatedBaselineType = 'manualSelection';
saveFigs = 'y';
for bb = 1:size(mergedDataFileIDs,1)
    mergedDataFileID = mergedDataFileIDs(bb,:);
    [figHandle] = GenerateSingleFigures_2P_eLife2025(mergedDataFileID,updatedBaselineType,saveFigs,RestingBaselines);
    close(figHandle)
end
%% Set Isoflurane boundary times
SetIsofluraneBoundary_2P_eLife2025(RestingBaselines)
disp('Stage Three Processing - Complete.'); disp(' ')
