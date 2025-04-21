function [] = CreateAllSpecDataStruct_IOS_eLife2025()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
neuralDataTypes = {'cortical_LH','cortical_RH','hippocampus'};
% character list of all SpecData files
specDataFileStruct = dir('*_SpecData.mat');
specDataFiles = {specDataFileStruct.name}';
specDataFileIDs = char(specDataFiles);
AllSpecData = [];
for aa = 1: size(specDataFileIDs,1)
    specDataFileID = specDataFileIDs(aa,:);
    [animalID,~,~] = GetFileInfo_IOS_eLife2025(specDataFileID);
    disp(['Adding spectrogram data to full struct for file number ' num2str(aa) ' of ' num2str(size(specDataFiles,1)) '...']); disp(' ')
    load(specDataFileID)
    for bb = 1:length(neuralDataTypes)
        neuralDataType = neuralDataTypes{1,bb};
        AllSpecData.(neuralDataType).fileIDs{aa,1} = specDataFileID;
        AllSpecData.(neuralDataType).S{aa,1} =  SpecData.(neuralDataType).S;
        AllSpecData.(neuralDataType).normS{aa,1} = SpecData.(neuralDataType).normS;
        AllSpecData.(neuralDataType).T{aa,1} =  SpecData.(neuralDataType).T;
        AllSpecData.(neuralDataType).F{aa,1} =  SpecData.(neuralDataType).F;
        AllSpecData.(neuralDataType).params =  SpecData.(neuralDataType).params;
        AllSpecData.(neuralDataType).movingwin =  SpecData.(neuralDataType).movingwin;
    end
end
AllSpecStructFileIDA = [animalID '_AllSpecStruct.mat'];
disp('Saving structure...'); disp(' ')
save(AllSpecStructFileIDA,'AllSpecData','-v7.3')