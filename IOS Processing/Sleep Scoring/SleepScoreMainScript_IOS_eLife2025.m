%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
zap;
[TrainingFiles,procDataFileIDs] = SelectTrainingDates_IOS_nNOS();
% add sleep parameters (each behavior we care about during sleep)
AddSleepParameters_IOS_nNOS(procDataFileIDs)
% create a table of values for sleep scoring model
CreateModelDataSet_IOS_nNOS(procDataFileIDs)
% create manual decisions for each 5 second bin
CreateTrainingDataSet_IOS_nNOS(procDataFileIDs,TrainingFiles)
% train Models - cycle through each data set and update any necessary parameters
TrainSleepModels_IOS_nNOS();
% sleep score data set and create structure for classification
[ScoringResults] = PredictBehaviorEvents_IOS_nNOS(modelName);
ApplySleepLogical_IOS_nNOS(modelName,TrainingFiles,ScoringResults)
CreateSleepData_IOS_nNOS(TrainingFiles);
