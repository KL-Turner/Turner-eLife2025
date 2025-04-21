function [] = ApplySleepLogical_IOS_eLife2025(modelName,TrainingFiles,ScoringResults)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
if strcmp(modelName,'Manual') == false
    % character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    % go through each ProcData file and apply a logical to each bin for this model
    for a = 1:size(procDataFileIDs,1)
        disp(['Adding sleep logicals to ProcData file (' num2str(a) '/' (num2str(size(procDataFileIDs,1))) ')...']); disp(' ')
        procDataFileID = procDataFileIDs(a,:);
        load(procDataFileID)
        fileID = procDataFileID(1:end - 13);
        c = 1;
        % extract the labels from the model for each file
        for b = 1:length(ScoringResults.allfileIDs)
            modelFileCheck = ScoringResults.allfileIDs{b,1};
            modelFileID = modelFileCheck(1:end - 14);
            if strcmp(fileID,modelFileID) == true
                behavState{c,1} = ScoringResults.alllabels{b,1}; %#ok<*AGROW>
                c = c + 1;
            end
        end
        % create a logical for each behavior bin
        for d = 1:length(behavState)
            if strcmp(behavState{d,1},'Not Sleep') == true
                awakeLogical(d,1) = 1;
                nremLogical(d,1) = 0;
                remLogical(d,1) = 0;
            elseif strcmp(behavState{d,1},'NREM Sleep') == true
                awakeLogical(d,1) = 0;
                nremLogical(d,1) = 1;
                remLogical(d,1) = 0;
            elseif strcmp(behavState{d,1},'REM Sleep') == true
                awakeLogical(d,1) = 0;
                nremLogical(d,1) = 0;
                remLogical(d,1) = 1;
            end
        end
        % save each logical under the sleep folder
        ProcData.sleep.logicals.(modelName).awakeLogical = logical(awakeLogical);
        ProcData.sleep.logicals.(modelName).nremLogical = logical(nremLogical);
        ProcData.sleep.logicals.(modelName).remLogical = logical(remLogical);
        save(procDataFileID,'ProcData')
    end
else
    % character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    cc = 1;
    % reduce file list to those with the training dates
    for aa = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(aa,:);
        [~,fileDate,~] = GetFileInfo_IOS_eLife2025(procDataFileID);
        if strcmp(fileDate,TrainingFiles.day1) == true || strcmp(fileDate,TrainingFiles.day2) == true
            trainingFileList(cc,:) = procDataFileID;
            cc = cc + 1;
        end
    end
    % character list of all TrainingData files
    trainingDataFileStruct = dir('*_TrainingData.mat');
    trainingDataFiles = {trainingDataFileStruct.name}';
    trainingDataFileIDs = char(trainingDataFiles);
    for e = 1:size(trainingFileList,1)
        procDataFileID = trainingFileList(e,:);
        load(procDataFileID)
        trainingDataFileID = trainingDataFileIDs(e,:);
        load(trainingDataFileID)
        behavState = trainingTable.behavState;
        % create a logical for each behavior bin
        for f = 1:length(behavState)
            if strcmp(behavState{f,1},'Not Sleep') == true
                awakeLogical(f,1) = 1;
                nremLogical(f,1) = 0;
                remLogical(f,1) = 0;
            elseif strcmp(behavState{f,1},'NREM Sleep') == true
                awakeLogical(f,1) = 0;
                nremLogical(f,1) = 1;
                remLogical(f,1) = 0;
            elseif strcmp(behavState{f,1},'REM Sleep') == true
                awakeLogical(f,1) = 0;
                nremLogical(f,1) = 0;
                remLogical(f,1) = 1;
            end
        end
        % save each logical under the sleep folder
        ProcData.sleep.logicals.(modelName).awakeLogical = logical(awakeLogical);
        ProcData.sleep.logicals.(modelName).nremLogical = logical(nremLogical);
        ProcData.sleep.logicals.(modelName).remLogical = logical(remLogical);
        save(procDataFileID,'ProcData')
    end
end