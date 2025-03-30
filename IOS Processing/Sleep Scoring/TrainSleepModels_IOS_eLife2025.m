function [] = TrainSleepModels_IOS_nNOS()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% load resting baseline file
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFiles = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFiles);
load(baselineDataFileID)
% go to animal's data location
dataLocation = [rootFolder delim 'Data' delim animalID delim 'Training Data'];
cd(dataLocation)
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% prepare pupil training data by updating parameters
AddPupilSleepParameters_IOS_nNOS(procDataFileIDs,RestingBaselines)
CreatePupilModelDataSet_IOS_nNOS(procDataFileIDs)
UpdatePupilTrainingDataSet_IOS_nNOS(procDataFileIDs)
% prepare physio training data by updating parameters
AddSleepParameters_IOS_nNOS(procDataFileIDs,RestingBaselines,'manualSelection')
CreateModelDataSet_IOS_nNOS(procDataFileIDs)
UpdateTrainingDataSets_IOS_nNOS(procDataFileIDs)
% training data file IDs
pupilTrainingDataFileStruct = dir('*_PupilTrainingData.mat');
pupilTrainingDataFiles = {pupilTrainingDataFileStruct.name}';
pupilTrainingDataFileIDs = char(pupilTrainingDataFiles);
% only use training files that match those of the pupil
for aa = 1:size(pupilTrainingDataFileIDs)
    trainingDataFileIDs(aa,:) = strrep(pupilTrainingDataFileIDs(aa,:),'Pupil','');
end
% load each updated training set and concatenate the data into table
pupilJoinedTable = []; physioJoinedTable = []; combinedJoinedTable = [];
for bb = 1:size(pupilTrainingDataFileIDs,1)
    % pupil table
    pupilTrainingTableFileID = pupilTrainingDataFileIDs(bb,:);
    load(pupilTrainingTableFileID)
    pupilJoinedTable = vertcat(pupilJoinedTable,pupilTrainingTable);
    % physio table
    trainingTableFileID = trainingDataFileIDs(bb,:);
    load(trainingTableFileID)
    physioJoinedTable = vertcat(physioJoinedTable,trainingTable);
    % combined table
    combinedJoinedTable = vertcat(combinedJoinedTable,horzcat(trainingTable(:,1:end - 1),pupilTrainingTable));
end
shuffleSeedA = randperm(size(pupilJoinedTable,1));
P = 0.70 ;
iterations = 100;
numTrees = 128;
paroptions = statset('UseParallel',true);
% separate the manual scores into 3 groups based on arousal classification
physioJoinedAwakeTable = []; physioJoinedNREMTable = []; physioJoinedREMTable = [];
physioRandomTable = physioJoinedTable(shuffleSeedA,:);
for aa = 1:size(physioRandomTable,1)
    if strcmp(physioRandomTable.behavState{aa,1},'Not Sleep') == true
        physioJoinedAwakeTable = vertcat(physioJoinedAwakeTable,physioRandomTable(aa,:));
    elseif strcmp(physioRandomTable.behavState{aa,1},'NREM Sleep') == true
        physioJoinedNREMTable = vertcat(physioJoinedNREMTable,physioRandomTable(aa,:));
    elseif strcmp(physioRandomTable.behavState{aa,1},'REM Sleep') == true
        physioJoinedREMTable = vertcat(physioJoinedREMTable,physioRandomTable(aa,:));
    end
end
physioTrainingTable = vertcat(physioJoinedAwakeTable(1:round(P*size(physioJoinedAwakeTable,1)),:), ...
    physioJoinedNREMTable(1:round(P*size(physioJoinedNREMTable,1)),:),...
    physioJoinedREMTable(1:round(P*size(physioJoinedREMTable,1)),:));
physioTrainingTable = physioTrainingTable(shuffleSeedB,:);
physioTestingTable = vertcat(physioJoinedAwakeTable(round(P*size(physioJoinedAwakeTable,1) + 1:end),:), ...
    physioJoinedNREMTable(round(P*size(physioJoinedNREMTable,1) + 1:end),:),...
    physioJoinedREMTable(round(P*size(physioJoinedREMTable,1) + 1:end),:));
physioTestingTable = physioTestingTable(shuffleSeedC,:);
% train on odd data
physioXtraining = physioTrainingTable(:,1:end - 1);
physioYtraining = physioTrainingTable(:,end);
% test on even data
physioXtesting = physioTestingTable(:,1:end - 1);
physioYtesting = physioTestingTable(:,end);
% random forest
bestAccuracy = 0;
for aa = 1:iterations
    RF_MDL = TreeBagger(numTrees,physioXtraining,physioYtraining,'Method','Classification','Surrogate','all','OOBPrediction','on','ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'},'Options',paroptions);
    % use the model to generate a set of predictions
    [physioTestingPredictions,~] = predict(RF_MDL,physioXtesting);
    % confusion chart
    RF_confMat = figure;
    CM = confusionchart(physioYtesting.behavState,physioTestingPredictions);
    CM.ColumnSummary = 'column-normalized';
    CM.RowSummary = 'row-normalized';
    CM.Title = [animalID ' Testing Data'];
    confVals = CM.NormalizedValues;
    % totalScoresA = sum(confVals([7,8,9]));
    % modelAccuracyA = (sum(confVals(9)/totalScoresA))*100;
    % totalScoresB = sum(confVals([3,6,9]));
    % modelAccuracyB = (sum(confVals(9)/totalScoresB))*100;
    % modelAccuracy = (modelAccuracyA + modelAccuracyB)/2;
    totalScores = sum(confVals(:));
    modelAccuracy = (sum(confVals([1,5,9])/totalScores))*100;
    CM.Title = {'Physio RF',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
    close(RF_confMat)
    if modelAccuracy > bestAccuracy
        bestAccuracy = modelAccuracy;
        bestMDL = RF_MDL;
    end
end
% determine the misclassification probability (for classification trees) for out-of-bag observations in the training data
outOfBagError = oobError(bestMDL,'Mode','Ensemble');
% use the model to generate a set of predictions
[physioTestingPredictions,~] = predict(bestMDL,physioXtesting);
% save labels for later confusion matrix
Results_PhysioSleepModel.(animalID).physio.mdl = bestMDL;
Results_PhysioSleepModel.(animalID).physio.outOfBagError = outOfBagError;
Results_PhysioSleepModel.(animalID).physio.trueTestingLabels = physioYtesting.behavState;
Results_PhysioSleepModel.(animalID).physio.predictedTestingLabels = physioTestingPredictions;