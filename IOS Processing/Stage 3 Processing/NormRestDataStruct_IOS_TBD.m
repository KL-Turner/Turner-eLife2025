function [RestData] = NormRestDataStruct_IOS(animalID,RestData,RestingBaselines,baselineType)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataTypes = fieldnames(RestData);
for dT = 1:length(dataTypes)
    dataType = char(dataTypes(dT));
    hemisphereDataTypes = fieldnames(RestData.(dataType));
    for hDT = 1:length(hemisphereDataTypes)
        hemDataType = char(hemisphereDataTypes(hDT));
        normData = {};
        [uniqueDays,~,~] = GetUniqueDays_IOS(RestData.(dataType).(hemDataType).fileDates);
        for uD = 1:length(uniqueDays)
            date = uniqueDays{uD};
            strDay = ConvertDate_IOS(date);
            [~,dayInds] = GetDayInds_IOS(RestData.(dataType).(hemDataType).fileDates,date);
            disp(['Normalizing ' (hemDataType) ' ' (dataType) ' for ' (strDay) '...']); disp(' ')
            % calculate the baseline differently depending on data type
            dayData = RestData.(dataType).(hemDataType).data(dayInds);
            normDayData = cell(size(dayData));
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
            for dD = 1:size(dayData,1)
                cellBase = dayBaseline*ones(1,size(dayData{dD},2));
                normDayData{dD} = dayData{dD}./cellBase - 1;
            end
            normData(dayInds) = normDayData;
        end
        RestData.(dataType).(hemDataType).NormData = normData';
    end
end
save([animalID '_RestData.mat'],'RestData','-v7.3')