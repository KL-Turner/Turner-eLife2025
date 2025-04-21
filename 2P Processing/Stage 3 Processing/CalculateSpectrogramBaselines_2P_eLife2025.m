function [RestingBaselines] = CalculateSpectrogramBaselines_2P_eLife2025(animal,neuralDataTypes,trialDuration_sec,specDataFiles,RestingBaselines,baselineType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Uses the resting time indeces to extract the average resting power in each frequency bin during periods of
%            rest to normalize the spectrogram data.
%________________________________________________________________________________________________________________________

for a = 1:length(neuralDataTypes)
    neuralDataType = neuralDataTypes{1,a};
    restFileList = unique(RestingBaselines.(baselineType).baselineFileInfo.fileIDs);      % Obtain the list of unique fileIDs
    restS1 = cell(size(restFileList,1),1);
    restS5 = cell(size(restFileList,1),1);
    % Obtain the spectrogram information from all the resting files
    for b = 1:length(restFileList)
        clear SpecData S1 T1 S5 T5
        fileID = restFileList{b,1};   % FileID of currently loaded file
        % Load in neural data from current file
        for c = 1:size(specDataFiles,1)
            [~,~,~,specDataFile,~,~] = GetFileInfo2_2P_eLife2025(specDataFiles(c,:));
            if strcmp(fileID,specDataFile)
                load(specDataFiles(c,:),'-mat')
                S1 = SpecData.(neuralDataType).oneSec.S;
                T1 = round(SpecData.(neuralDataType).oneSec.T,1);
                S5 = SpecData.(neuralDataType).fiveSec.S;
                T5 = round(SpecData.(neuralDataType).fiveSec.T,3);
            end
        end
        restS1{b,1} = S1;
        restS5{b,1} = S5;
    end
    for d = 1:length(restFileList)
        fileID = restFileList{d,1};
        strDay = ConvertDate_2P_eLife2025(fileID(1:6));
        S1_data = restS1{d,1};
        S5_data = restS5{d,1};
        binSize1 = 30;
        binSize5 = 5;
        S1_trialRest = [];
        S5_trialRest = [];
        for e = 1:length(RestingBaselines.(baselineType).baselineFileInfo.fileIDs)
            restFileID = RestingBaselines.(baselineType).baselineFileInfo.fileIDs{e,1};
            if strcmp(fileID,restFileID)
                restDuration1 = round(RestingBaselines.(baselineType).baselineFileInfo.durations(e,1),1);
                restDuration5 = round(RestingBaselines.(baselineType).baselineFileInfo.durations(e,1),1);
                startTime1 = double(round(RestingBaselines.(baselineType).baselineFileInfo.eventTimes(e,1),1));
                startTime5 = round(RestingBaselines.(baselineType).baselineFileInfo.eventTimes(e,1),1);
                % 1 second spectrogram conditions and indexing
                if startTime1 >= 0.5 && (startTime1 + restDuration1) <= (trialDuration_sec - 0.5)
                    startTime1_index = find(T1 == startTime1);
                    startTime1_index = startTime1_index(1);
                    restDuration1_Index = restDuration1*binSize1 - 1;
                    restDuration1_Index = restDuration1_Index(1);
                    S1_single_rest = S1_data(:,(startTime1_index:(startTime1_index + restDuration1_Index)));
                    S1_trialRest = [S1_single_rest,S1_trialRest]; %#ok<*AGROW>
                end
                % 5 second spectrogram conditions and indexing
                if startTime5 >= 2.5 && (startTime5 + restDuration5) <= (trialDuration_sec - 2.5)
                    [~,startTime5_index] = min(abs(T5 - startTime5));
                    restDuration5_Index = floor(restDuration5*binSize5) - 1;
                    S5_single_rest = S5_data(:,(startTime5_index:(startTime5_index + restDuration5_Index)));
                    S5_trialRest = [S5_single_rest,S5_trialRest];
                end
            end
        end
        S_trialAvg1 = mean(S1_trialRest,2);
        S_trialAvg5 = mean(S5_trialRest,2);
        trialRestData.([strDay '_' fileID]).oneSec.S_avg = S_trialAvg1;
        trialRestData.([strDay '_' fileID]).fiveSec.S_avg = S_trialAvg5;
    end
    fields = fieldnames(trialRestData);
    uniqueDays = GetUniqueDays_2P_eLife2025(RestingBaselines.(baselineType).baselineFileInfo.fileIDs);
    for f = 1:length(uniqueDays)
        g = 1;
        for field = 1:length(fields)
            if strcmp(fields{field}(7:12), uniqueDays{f})
                stringDay = ConvertDate_2P_eLife2025(uniqueDays{f});
                S_avgs.oneSec.(stringDay){g,1} = trialRestData.(fields{field}).oneSec.S_avg;
                S_avgs.fiveSec.(stringDay){g,1} = trialRestData.(fields{field}).fiveSec.S_avg;
                g = g + 1;
            end
        end
    end
    dayFields = fieldnames(S_avgs.oneSec);
    for h = 1:length(dayFields)
        dayVals1 = [];
        dayVals5 = [];
        for j = 1:length(S_avgs.oneSec.(dayFields{h}))
            dayVals1 = [dayVals1,S_avgs.oneSec.(dayFields{h}){j,1}];
            dayVals5 = [dayVals5,S_avgs.fiveSec.(dayFields{h}){j,1}];
        end
        disp(['Adding spectrogram baseline to baseline file for ' neuralDataType ' on ' dayFields{h} '...']); disp(' ')
        RestingBaselines.Spectrograms.(neuralDataType).oneSec.(dayFields{h}) = mean(dayVals1,2);
        RestingBaselines.Spectrograms.(neuralDataType).fiveSec.(dayFields{h}) = mean(dayVals5,2);
    end
end
save([animal '_RestingBaselines.mat'],'RestingBaselines');

end