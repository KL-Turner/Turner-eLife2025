function [SpecData] = NormalizeSpectrograms_IOS_eLife2025(RestingBaselines)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
neuralDataTypes = {'cortical_LH','cortical_RH','hippocampus'};
% file list A
specDirectory = dir('*_SpecData.mat');
specDataFiles = {specDirectory.name}';
specDataFileIDs = char(specDataFiles);
for aa = 1:size(specDataFileIDs,1)
    disp(['Normalizing spectrogram file ' num2str(aa) ' of ' num2str(size(specDataFileIDs,1)) '...']); disp(' ')
    load(specDataFileIDs(aa,:),'-mat');
    [~,fileDate,~] = GetFileInfo_IOS_eLife2025(specDataFileIDs(aa,:));
    strDay = ConvertDate_IOS_eLife2025(fileDate);
    for bb = 1:length(neuralDataTypes)
        neuralDataType = neuralDataTypes{1,bb};
        baseLine = RestingBaselines.Spectrograms.(neuralDataType).(strDay);
        S = SpecData.(neuralDataType).S;
        holdMatrix = baseLine.*ones(size(S));
        SpecData.(neuralDataType).normS = (S - holdMatrix)./holdMatrix;
    end
    save(specDataFileIDs(aa,:),'SpecData')
end