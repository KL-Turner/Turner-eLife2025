function [] = ApplySleepLogical_2P_nNOS(modelName)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% character list of all MergedData files
mergedDataFileStruct = dir('*_MergedData.mat');
mergedDataFiles = {mergedDataFileStruct.name}';
mergedDataFileIDs = char(mergedDataFiles);
% go through each MergedData file and apply a logical to each bin for this model
for aa = 1:size(mergedDataFileIDs,1)
    disp(['Adding sleep logicals to MergedData file (' num2str(aa) '/' (num2str(size(mergedDataFileIDs,1))) ')...']); disp(' ')
    mergedDataFileID = mergedDataFileIDs(aa,:);
    load(mergedDataFileID)
    trainingFileID = [mergedDataFileID(1:end - 14) 'TrainingData.mat'];
    load(trainingFileID)
    behavState = trainingTable.behavState;
    % create a logical for each behavior bin
    for dd = 1:length(behavState)
        if strcmp(behavState{dd,1},'Not Sleep') == true
            awakeLogical(dd,1) = 1; %#ok<*AGROW>
            nremLogical(dd,1) = 0;
            remLogical(dd,1) = 0;
        elseif strcmp(behavState{dd,1},'NREM Sleep') == true
            awakeLogical(dd,1) = 0;
            nremLogical(dd,1) = 1;
            remLogical(dd,1) = 0;
        elseif strcmp(behavState{dd,1},'REM Sleep') == true
            awakeLogical(dd,1) = 0;
            nremLogical(dd,1) = 0;
            remLogical(dd,1) = 1;
        end
    end
    % save each logical under the sleep folder
    MergedData.sleep.logicals.(modelName).awakeLogical = logical(awakeLogical);
    MergedData.sleep.logicals.(modelName).nremLogical = logical(nremLogical);
    MergedData.sleep.logicals.(modelName).remLogical = logical(remLogical);
    save(mergedDataFileID,'MergedData')
end

end
