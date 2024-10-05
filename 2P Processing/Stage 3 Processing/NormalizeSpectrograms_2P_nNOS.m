function [SpecData] = NormalizeSpectrograms_2P(specDataFiles,neuralDataTypes,RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Normalizes each spectrogram by the resting baseline for that day.
%________________________________________________________________________________________________________________________

for a = 1:size(specDataFiles,1)
    disp(['Normalizing spectrogram file ' num2str(a) ' of ' num2str(size(specDataFiles,1)) '...']); disp(' ')
    load(specDataFiles(a,:),'-mat');
    [~,~,fileDate,~,~,~] = GetFileInfo2_2P(specDataFiles(a,:));
    strDay = ConvertDate_2P(fileDate);
    for b = 1:length(neuralDataTypes)
        neuralDataType = neuralDataTypes{1,b};
        baseLine1 = RestingBaselines.Spectrograms.(neuralDataType).oneSec.(strDay);
        baseLine5 = RestingBaselines.Spectrograms.(neuralDataType).fiveSec.(strDay);
        S1 = SpecData.(neuralDataType).oneSec.S;
        S5 = SpecData.(neuralDataType).fiveSec.S;       
        holdMatrix1 = baseLine1.*ones(size(S1));
        holdMatrix5 = baseLine5.*ones(size(S5));        
        normS1 = (S1 - holdMatrix1)./holdMatrix1;
        normS5 = (S5 - holdMatrix5)./holdMatrix5;
        SpecData.(neuralDataType).oneSec.normS = normS1;
        SpecData.(neuralDataType).fiveSec.normS = normS5;
    end
    save(specDataFiles(a,:),'SpecData')
end

end
