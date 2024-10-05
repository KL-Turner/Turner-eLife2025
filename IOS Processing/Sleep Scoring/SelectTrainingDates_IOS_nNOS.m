function [TrainingFiles,procDataFileIDs] = SelectTrainingDates_IOS_nNOS()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% check and load TrainingFileDates
trainingDatesFileStruct = dir('*_TrainingFileDates.mat');
trainingDatesFile = {trainingDatesFileStruct.name}';
trainingDatesFileID = char(trainingDatesFile);
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
if isempty(trainingDatesFileID) == true
    [animalIDs,fileDates,~] = GetFileInfo_IOS_nNOS(procDataFileIDs);
    animalID = animalIDs(1,:);
    uniqueDays = GetUniqueDays_IOS_nNOS(fileDates);
    TrainingFiles = [];
    % select day 1
    tf = 0;
    while tf == 0
        disp('Select the first day for training data'); disp(' ')
        [indxA,tf] = listdlg('PromptString',{'Select a file.','Only one file can be selected at a time.',''},'SelectionMode','single','ListString',uniqueDays);
        if tf == 0
            disp('Please select a date'); disp(' ');
        else
            disp(['Day ' num2str(uniqueDays{indxA,1}) ' selected (' ConvertDate_IOS_nNOS(uniqueDays{indxA,1}) ')']); disp(' ')
            TrainingFiles.day1 = uniqueDays{indxA,1};
        end
    end
    % select day 2
    tf = 0;
    while tf == 0
        disp('Select the second day for training data'); disp(' ')
        [indxB,tf] = listdlg('PromptString',{'Select a file.','Only one file can be selected at a time.',''},'SelectionMode','single','ListString',uniqueDays);
        if tf == 0
            disp('Please select a date'); disp(' ');
        else
            disp(['Day ' num2str(uniqueDays{indxB,1}) ' selected (' ConvertDate_IOS_nNOS(uniqueDays{indxB,1}) ')']); disp(' ')
            TrainingFiles.day2 = uniqueDays{indxB,1};
        end
    end
    % save
    save([animalID '_TrainingFileDates.mat'],'TrainingFiles')
else
    load(trainingDatesFileID,'-mat')
end