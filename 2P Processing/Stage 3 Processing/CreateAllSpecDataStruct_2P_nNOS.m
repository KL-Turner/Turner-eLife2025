function [] = CreateAllSpecDataStruct_2P(animalID,neuralDataTypes)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Normalizes each spectrogram by the resting baseline for that day.
%________________________________________________________________________________________________________________________

AllSpecStructFileID = [animalID '_AllSpecStruct.mat'];
% if ~exist(AllSpecStructFileID)
    % Character list of all SpecData files
    specDataFileStruct = dir('*_SpecData.mat');
    specDataFiles = {specDataFileStruct.name}';
    specDataFileIDs = char(specDataFiles);
    for c = 1: size(specDataFileIDs,1)
        specDataFileID = specDataFileIDs(c,:);
        disp(['Adding spectrogram data to full struct for file number ' num2str(c) ' of ' num2str(size(specDataFiles, 1)) '...']); disp(' ')
        load(specDataFileID)
        for d = 1:length(neuralDataTypes)
            neuralDataType = neuralDataTypes{1,d};
            AllSpecData.(neuralDataType).fileIDs{c,1} = specDataFileID;
            AllSpecData.(neuralDataType).fiveSec.S{c,1} =  SpecData.(neuralDataType).fiveSec.S;
            AllSpecData.(neuralDataType).fiveSec.normS{c,1} = SpecData.(neuralDataType).fiveSec.normS;
            AllSpecData.(neuralDataType).fiveSec.T{c,1} =  SpecData.(neuralDataType).fiveSec.T;
            AllSpecData.(neuralDataType).fiveSec.F{c,1} =  SpecData.(neuralDataType).fiveSec.F;
            AllSpecData.(neuralDataType).fiveSec.params =  SpecData.(neuralDataType).fiveSec.params;
            AllSpecData.(neuralDataType).fiveSec.movingwin =  SpecData.(neuralDataType).fiveSec.movingwin;
            AllSpecData.(neuralDataType).oneSec.S{c,1} = SpecData.(neuralDataType).oneSec.S;
            AllSpecData.(neuralDataType).oneSec.normS{c,1} = SpecData.(neuralDataType).oneSec.normS;
            AllSpecData.(neuralDataType).oneSec.T{c,1} = SpecData.(neuralDataType).oneSec.T;
            AllSpecData.(neuralDataType).oneSec.F{c,1} = SpecData.(neuralDataType).oneSec.F;
            AllSpecData.(neuralDataType).oneSec.params = SpecData.(neuralDataType).oneSec.params;
            AllSpecData.(neuralDataType).oneSec.movingwin = SpecData.(neuralDataType).oneSec.movingwin;
        end
    end
    disp('Saving structure...'); disp(' ')
    save(AllSpecStructFileID,'AllSpecData','-v7.3')
% end

end
