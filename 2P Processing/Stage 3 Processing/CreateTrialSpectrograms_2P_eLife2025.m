function CreateTrialSpectrograms_2P_nNOS(mergedDataFileIDs,specNeuralDataTypes)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Analyzes the raw neural data from each RawData.mat file and calculates two different spectrograms.
%________________________________________________________________________________________________________________________

for a = 1:size(mergedDataFileIDs, 1)
    specDataFileID = [mergedDataFileIDs(a,1:end - 15) '_SpecData.mat'];
    if ~exist(specDataFileID,'file') == true
        mergedDataFileID = mergedDataFileIDs(a,:);
        load(mergedDataFileID);
        duration = MergedData.notes.trialDuration_Sec;
        anFs = MergedData.notes.anFs;
        expectedLength = duration*anFs;
        for z = 1:length(specNeuralDataTypes)
            specNeuralDataType = specNeuralDataTypes{1,z};
            try
                rawNeuro = detrend(MergedData.data.(specNeuralDataType)(1:expectedLength),'constant');
            catch
                sampleDiff = expectedLength - length(MergedData.data.(specNeuralDataType));
                rawNeuro = detrend(horzcat(MergedData.data.(specNeuralDataType),MergedData.data.(specNeuralDataType)(end)*ones(1,sampleDiff)),'constant');
            end
            w0 = 60/(anFs/2);
            bw = w0/35;
            [num,den] = iirnotch(w0,bw);
            rawNeuro2 = filtfilt(num,den,rawNeuro);
            % Spectrogram parameters
            % 1 second spectrogram
            params1.tapers = [5,9];
            params1.Fs = anFs;
            params1.fpass = [1,100];
            movingwin1 = [1,1/30];
            % 5 second spectrogram
            params5.tapers = [5,9];
            params5.Fs = anFs;
            params5.fpass = [1,100];
            movingwin5 = [5,1/5];
            % analyze each spectrogram based on parameters
            disp(['Creating spectrogram for file number ' num2str(a) ' of ' num2str(size(mergedDataFileIDs,1)) '...']); disp(' ')
            [S1,T1,F1] = mtspecgramc(rawNeuro2,movingwin1,params1);
            [S5,T5,F5] = mtspecgramc(rawNeuro2,movingwin5,params5);
            % save data ins tructure
            if strcmp(specNeuralDataType,'rawCorticalNeural') == true
                neuralDataType = 'corticalNeural';
            elseif strcmp(specNeuralDataType,'rawHippocampalNeural') == true
                neuralDataType = 'hippocampalNeural';
            end
            SpecData.(neuralDataType).fiveSec.S = S5';
            SpecData.(neuralDataType).fiveSec.T = T5;
            SpecData.(neuralDataType).fiveSec.F = F5;
            SpecData.(neuralDataType).fiveSec.params = params5;
            SpecData.(neuralDataType).fiveSec.movingwin = movingwin5;
            SpecData.(neuralDataType).oneSec.S = S1';
            SpecData.(neuralDataType).oneSec.T = T1;
            SpecData.(neuralDataType).oneSec.F = F1;
            SpecData.(neuralDataType).oneSec.params = params1;
            SpecData.(neuralDataType).oneSec.movingwin = movingwin1;
            save(specDataFileID,'SpecData');
        end
    else
        disp('File exists. continuing...'); disp(' ')
    end
end

end
