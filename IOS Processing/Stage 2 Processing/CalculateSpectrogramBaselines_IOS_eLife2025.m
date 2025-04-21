function [RestingBaselines] = CalculateSpectrogramBaselines_IOS_eLife2025(RestingBaselines,baselineType)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
neuralDataTypes = {'cortical_LH','cortical_RH','hippocampus'};
animalID = RestingBaselines.setDuration.baselineFileInfo.animalID;
trialDuration_sec = RestingBaselines.setDuration.baselineFileInfo.trialDuration_sec;
for a = 1:length(neuralDataTypes)
    neuralDataType = neuralDataTypes{1,a};
    restFileList = unique(RestingBaselines.(baselineType).baselineFileInfo.fileIDs); % obtain the list of unique fileIDs
    restS5A = cell(size(restFileList,1),1);
    % obtain the spectrogram information from all the resting files
    for b = 1:length(restFileList)
        fileID = restFileList{b,1};   % FileID of currently loaded file
        % load in neural data from current file
        clear SpecData
        specDataFileIDA = [animalID '_' fileID '_SpecData.mat'];
        load(specDataFileIDA,'-mat')
        S5A = SpecData.(neuralDataType).S;
        T5A = round(SpecData.(neuralDataType).T,3);
        restS5A{b,1} = S5A;
    end
    for d = 1:length(restFileList)
        fileID = restFileList{d,1};
        strDay = ConvertDate_IOS_eLife2025(fileID(1:6));
        S5C_data = restS5A{d,1};
        binSize5A = 5;
        S5A_trialRest = [];
        for e = 1:length(RestingBaselines.(baselineType).baselineFileInfo.fileIDs)
            restFileID = RestingBaselines.(baselineType).baselineFileInfo.fileIDs{e,1};
            if strcmp(fileID,restFileID)
                restDuration5A = round(RestingBaselines.(baselineType).baselineFileInfo.durations(e,1),1);
                startTime5A = round(RestingBaselines.(baselineType).baselineFileInfo.eventTimes(e,1),1);
                % 5 second spectrogram conditions and indexing
                if startTime5A >= 2.5 && (startTime5A + restDuration5A) <= (trialDuration_sec - 2.5)
                    [~,startTime5A_index] = min(abs(T5A - startTime5A));
                    restDuration5A_Index = floor(restDuration5A*binSize5A) - 1;
                    S5A_single_rest = S5C_data(:,(startTime5A_index:(startTime5A_index + restDuration5A_Index)));
                    S5A_trialRest = [S5A_single_rest,S5A_trialRest];
                end
            end
        end
        S_trialAvg5A = mean(S5A_trialRest,2);
        trialRestData.([strDay '_' fileID]).fiveSecA.S_avg = S_trialAvg5A;
    end
    fields = fieldnames(trialRestData);
    uniqueDays = GetUniqueDays_IOS_eLife2025(RestingBaselines.(baselineType).baselineFileInfo.fileIDs);
    for f = 1:length(uniqueDays)
        g = 1;
        for field = 1:length(fields)
            if strcmp(fields{field}(7:12),uniqueDays{f})
                stringDay = ConvertDate_IOS_eLife2025(uniqueDays{f});
                S_avgs.fiveSecA.(stringDay){g,1} = trialRestData.(fields{field}).fiveSecA.S_avg;
                g = g + 1;
            end
        end
    end
    dayFields = fieldnames(S_avgs.fiveSecA);
    for h = 1:length(dayFields)
        dayVals5A = [];
        for j = 1:length(S_avgs.fiveSecA.(dayFields{h}))
            dayVals5A = [dayVals5A,S_avgs.fiveSecA.(dayFields{h}){j,1}];
        end
        disp(['Adding spectrogram baseline to baseline file for ' neuralDataType ' on ' dayFields{h} '...']); disp(' ')
        RestingBaselines.Spectrograms.(neuralDataType).(dayFields{h}) = mean(dayVals5A,2);
    end
end
save([animalID '_RestingBaselines.mat'],'RestingBaselines');

end
