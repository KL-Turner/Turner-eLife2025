%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
%% Clear workspace/Load in file names for various analysis
clear; clc; close all
disp('Loading necessary file names...'); disp(' ')
animalIDs = {'T115','T116','T117','T118','T125','T126'};
baselineType = 'manualSelection';
startingDirectory = cd;

%% create manually-scored training data set for each animal
for aa = 1:size(animalIDs,2)
    % cd to the animal's training set folder
    trainingDirectory = [animalIDs{1,aa} '\2P Data\'];
    cd(trainingDirectory)
    % load the baseline structure
    baselinesFileStruct = dir('*_RestingBaselines.mat');
    baselinesFile = {baselinesFileStruct.name}';
    baselinesFileID = char(baselinesFile);
    load(baselinesFileID)
    % character list of all MergedData files
    mergedDataFileStruct = dir('*_MergedData.mat');
    mergedDataFiles = {mergedDataFileStruct.name}';
    mergedDataFileIDs = char(mergedDataFiles);
    % add sleep parameters (each behavior we care about during sleep)
    AddSleepParameters_2P(mergedDataFileIDs,RestingBaselines,baselineType)
    % create manual decisions for each 5 second bin
    CreateTrainingDataSet_2P(mergedDataFileIDs,RestingBaselines,baselineType)
    cd(startingDirectory)
end

%% sleep score an animal's data set and create a SleepData.mat structure for classification
modelName = 'Manual';
for bb = 1:length(animalIDs)
    SleepData = [];
    % cd to the animal's training set folder
    trainingDirectory = [animalIDs{1,bb} '\2P Data\'];
    cd(trainingDirectory)
    NREMsleepTime = 30;   % seconds
    REMsleepTime = 60;   % seconds
    % apply sleep logicals for sleep scoring
    ApplySleepLogical_2P(modelName)
    % create SleepData.mat structure
    [SleepData] = CreateSleepData_2P(NREMsleepTime,REMsleepTime,modelName,SleepData);
    save([animalIDs{1,bb} '_SleepData.mat'],'SleepData')
    cd(startingDirectory)
end
