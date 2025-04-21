%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
zap;
[TrainingFiles,procDataFileIDs] = SelectTrainingDates_IOS_eLife2025();
% add sleep parameters (each behavior we care about during sleep)
AddSleepParameters_IOS_eLife2025(procDataFileIDs)
% create a table of values for sleep scoring model
CreateModelDataSet_IOS_eLife2025(procDataFileIDs)
% create manual decisions for each 5 second bin
CreateTrainingDataSet_IOS_eLife2025(procDataFileIDs,TrainingFiles)
% train Models - cycle through each data set and update any necessary parameters
TrainSleepModels_IOS_eLife2025();
% sleep score data set and create structure for classification
[ScoringResults] = PredictBehaviorEvents_IOS_eLife2025(modelName);
ApplySleepLogical_IOS_eLife2025(modelName,TrainingFiles,ScoringResults)
CreateSleepData_IOS_eLife2025(TrainingFiles);
