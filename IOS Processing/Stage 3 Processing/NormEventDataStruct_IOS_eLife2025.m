function [EventData] = NormEventDataStruct_IOS_eLife2025(animalID,EventData,RestingBaselines,baselineType)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataTypes = fieldnames(EventData);
for dT = 1:length(dataTypes)
    dataType = char(dataTypes(dT));
    hemisphereDataTypes = fieldnames(EventData.(dataType));
    for hDT = 1:length(hemisphereDataTypes)
        hemDataType = char(hemisphereDataTypes(hDT));
        behaviorFields = fieldnames(EventData.(dataType).(hemDataType));
        for bF = 1:length(behaviorFields)
            behavField = char(behaviorFields(bF));
            normData = [];
            if isempty(EventData.(dataType).(hemDataType).(behavField).data) == false
                [uniqueDays,~,~] = GetUniqueDays_IOS_eLife2025(EventData.(dataType).(hemDataType).(behavField).fileDates);
                for uD = 1:length(uniqueDays)
                    date = uniqueDays{uD};
                    strDay = ConvertDate_IOS_eLife2025(date);
                    [~,dayInds] = GetDayInds_IOS_eLife2025(EventData.(dataType).(hemDataType).(behavField).fileDates,date);
                    disp(['Normalizing ' (hemDataType) ' ' (dataType) ' ' (behavField) ' for ' (strDay) '...']); disp(' ')
                    % calculate the baseline differently depending on data type
                    try
                        dayBaseline = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                    catch
                        if strcmp(hemDataType,'LH_gammaBandPower') == true
                            dayBaseline = RestingBaselines.(baselineType).cortical_LH.gammaBandPower.(strDay).mean;
                        elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                            dayBaseline = RestingBaselines.(baselineType).cortical_RH.gammaBandPower.(strDay).mean;
                        else
                            dayBaseline = NaN;
                        end
                    end
                    % pre-allocate array and use for permutation
                    normDayData = EventData.(dataType).(hemDataType).(behavField).data(dayInds,:,:);
                    % permute norm_session_data to handle both matrix and array (squeeze
                    % causes a matrix dimension error if not permuted)
                    dayData = permute(normDayData,unique([2,1,ndims(normDayData)],'stable'));
                    for dD = 1:size(dayData,2)
                        normDayData(dD,:,:) = squeeze(dayData(:,dD,:))./(ones(size(dayData,1),1)*dayBaseline) - 1;
                    end
                    normData(dayInds,:,:) = normDayData;
                end
                EventData.(dataType).(hemDataType).(behavField).NormData = normData;
            end
        end
    end
end
save([animalID '_EventData.mat'],'EventData','-v7.3')