function [Results_ModelAccuracy_Ephys] = AnalyzeModelAccuracy_Ephys_nNOS(animalID,group,set,rootFolder,delim,Results_ModelAccuracy_Ephys)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% training data file IDs
trainingDataFileStruct = dir('*_TrainingData.mat');
trainingDataFiles = {trainingDataFileStruct.name}';
trainingDataFileIDs = char(trainingDataFiles);
% load each updated training set and concatenate the data into table
physioJoinedTable = [];
for bb = 1:size(trainingDataFileIDs,1)
    % physio table
    trainingTableFileID = trainingDataFileIDs(bb,:);
    load(trainingTableFileID)
    physioJoinedTable = vertcat(physioJoinedTable,trainingTable);
end
shuffleSeedA = randperm(size(physioJoinedTable,1));
P = 0.70 ;
iterations = 100;
numTrees = 128;
paroptions = statset('UseParallel',true);
% physio model - separate the manual scores into 3 groups based on arousal classification
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
shuffleSeedB = randperm(size(physioTrainingTable,1));
physioTrainingTable = physioTrainingTable(shuffleSeedB,:);
physioTestingTable = vertcat(physioJoinedAwakeTable(round(P*size(physioJoinedAwakeTable,1) + 1:end),:), ...
    physioJoinedNREMTable(round(P*size(physioJoinedNREMTable,1) + 1:end),:),...
    physioJoinedREMTable(round(P*size(physioJoinedREMTable,1) + 1:end),:));
shuffleSeedC = randperm(size(physioTestingTable,1));
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
Results_ModelAccuracy_Ephys.(group).(animalID).mdl = bestMDL;
Results_ModelAccuracy_Ephys.(group).(animalID).outOfBagError = outOfBagError;
Results_ModelAccuracy_Ephys.(group).(animalID).trueTestingLabels = physioYtesting.behavState;
Results_ModelAccuracy_Ephys.(group).(animalID).predictedTestingLabels = physioTestingPredictions;
% save results
cd([rootFolder delim 'Results_Turner'])
save('Results_ModelAccuracy_Ephys.mat','Results_ModelAccuracy_Ephys')
cd([rootFolder delim 'Data'])