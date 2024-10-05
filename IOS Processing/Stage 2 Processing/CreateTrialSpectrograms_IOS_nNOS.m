function [] = CreateTrialSpectrograms_IOS_nNOS()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
neuralDataTypes = {'cortical_LH','cortical_RH','hippocampus'};
for aa = 1:size(rawDataFileIDs,1)
    rawDataFile = rawDataFileIDs(aa,:);
    clear RawData
    [animalID,~,fileID] = GetFileInfo_IOS_nNOS(rawDataFile);
    specDataFileID = [animalID '_' fileID '_SpecData.mat'];
    % 5 second spectrograms with 1/5 second step size
    if ~exist(specDataFileID,'file') == true
        SpecData = [];
        if exist('RawData','file') == false
            load(rawDataFile);
            duration = RawData.notes.trialDuration_sec;
            analogFs = RawData.notes.analogSamplingRate;
            expectedLength = duration*analogFs;
        end
        disp(['Analyzing spectrograms for file (' num2str(aa) '/' num2str(size(rawDataFileIDs,1)) ')']); disp(' ')
        for bb = 1:length(neuralDataTypes)
            neuralDataType = neuralDataTypes{1,bb};
            try
                rawNeuro = detrend(RawData.data.(neuralDataType)(1:expectedLength),'constant');
            catch
                sampleDiff = expectedLength - length(RawData.data.(neuralDataType));
                rawNeuro = detrend(horzcat(RawData.data.(neuralDataType),RawData.data.(neuralDataType)(end)*ones(1,sampleDiff)),'constant');
            end
            % parameters
            params.tapers = [5,9];
            params.Fs = analogFs;
            params.fpass = [1,100];
            movingwin = [5,1/5];
            % analyze each spectrogram based on parameters
            [S,T,F] = mtspecgramc(rawNeuro,movingwin,params);
            % save data ins tructure
            SpecData.(neuralDataType).S = S';
            SpecData.(neuralDataType).T = T;
            SpecData.(neuralDataType).F = F;
            SpecData.(neuralDataType).params = params;
            SpecData.(neuralDataType).movingwin = movingwin;
            save(specDataFileID,'SpecData');
        end
        save(specDataFileID,'SpecData');
    else
        disp(['Spectrograms for file (' num2str(aa) '/' num2str(size(rawDataFileIDs,1)) ') already analyzed.']); disp(' ')
    end
end