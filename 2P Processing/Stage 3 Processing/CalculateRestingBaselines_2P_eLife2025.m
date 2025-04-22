function [RestingBaselines] = CalculateRestingBaselines_2P_eLife2025(animalID,targetMinutes,trialDuration_sec,RestData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: This function finds the resting baseline for all fields of the RestData.mat structure, for each unique day
%________________________________________________________________________________________________________________________

disp(['Calculating the resting baselines for the first ' num2str(targetMinutes) ' minutes of each unique day...']); disp(' ')
% The RestData.mat struct has all resting events, regardless of duration. We want to set the threshold for rest as anything
% that is greater than 10 seconds.
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {2.5};
PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};
dataTypes = fieldnames(RestData);
for aa = 1:length(dataTypes)
    dataType = char(dataTypes(aa));   % Load each loop iteration's fieldname as a character string
    if strcmp(dataType,'vesselDiameter') == true
        subDataTypes = fieldnames(RestData.(dataType));   % Find the hemisphere dataTypes. These are typically LH, RH
        % Loop through each hemisphere dataType (LH, RH) because they are subfields and will have unique baselines
        for bb = 1:length(subDataTypes)
            subDataType = char(subDataTypes(bb));   % Load each loop iteration's hemisphere fieldname as a character string
            [restLogical] = FilterEvents_2P_eLife2025(RestData.(dataType).(subDataType),RestCriteria);
            [puffLogical] = FilterEvents_2P_eLife2025(RestData.(dataType).(subDataType),PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            allRestData = RestData.(dataType).(subDataType).data(combRestLogical,:);
            allRestFileIDs = RestData.(dataType).(subDataType).fileIDs(combRestLogical,:);
            allRestDurations = RestData.(dataType).(subDataType).durations(combRestLogical,:);
            allRestEventTimes = RestData.(dataType).(subDataType).eventTimes(combRestLogical,:);
            allRestVesselIDs = RestData.(dataType).(subDataType).vesselIDs(combRestLogical,:);
            uniqueVessels = unique(allRestVesselIDs);   % Total number of unique days in this folder under a single animal
            combDirectory = dir('*_MergedData.mat');
            mergedDataFiles = {combDirectory.name}';
            mergedDataFiles = char(mergedDataFiles);
            for cc = 1:size(mergedDataFiles)
                mergedDataFile = mergedDataFiles(cc,:);
                [animalID,~,~,fileID,~,vesselID] = GetFileInfo2_2P_eLife2025(mergedDataFile);
                allFileDates{cc,1} = fileID; %#ok<*AGROW>
                allVesselIDs{cc,1} = vesselID;
            end
            uniqueDates = GetUniqueDays_2P_eLife2025(allFileDates);
            for dd = 1:length(uniqueDates)
                uniqueDays{dd,1} = ConvertDate_2P_eLife2025(uniqueDates{dd,1});
            end
            tempBaseFileIDs = {};
            tempBaseEventTimes = [];
            tempBaseDurations = [];
            % Find the fieldnames of RestData and loop through each field. Each fieldname should be a different dataType of interest.
            % These will typically be CBV, Delta, Theta, Gamma, and MUA
            for ee = 1:length(uniqueVessels)
                uV = uniqueVessels{ee,1};
                for ff = 1:length(allVesselIDs)
                    checkID = allVesselIDs{ff,1};
                    if strcmp(uV,checkID)
                        validFiles(ff,1) = 1;
                    else
                        validFiles(ff,1) = 0;
                    end
                end
                validFiles = logical(validFiles);
                validDays = allFileDates(validFiles);
                for gg = 1:length(uniqueDays)
                    uD = uniqueDays{gg};
                    vesselFilterLogical = zeros(length(allRestVesselIDs),1);   % Find all the vessels that correspond to this loop's vessel (A1, A2, A3, etc)
                    for hh = 1:length(allRestVesselIDs)
                        vessel = char(allRestVesselIDs(hh,1));
                        if strcmp(uV,vessel)
                            vesselFilterLogical(hh,1) = 1;
                        end
                    end
                    vesselFilterLogical = logical(vesselFilterLogical);
                    singleVesselFileIDs = allRestFileIDs(vesselFilterLogical);
                    singleVesselDurations = allRestDurations(vesselFilterLogical);
                    singleVesselEventTimes = allRestEventTimes(vesselFilterLogical);
                    singleVesselData = allRestData(vesselFilterLogical);
                    dayFilterLogical = zeros(length(singleVesselFileIDs),1);
                    for ii = 1:length(singleVesselFileIDs)
                        day = ConvertDate_2P_eLife2025(singleVesselFileIDs{ii,1}(1:6));
                        if strcmp(uD,day)
                            dayFilterLogical(ii,1) = 1;
                        end
                    end
                    dayFilterLogical = logical(dayFilterLogical);
                    uniqueDayVesselFileIDs = singleVesselFileIDs(dayFilterLogical);
                    uniqueDayVesselDurations = singleVesselDurations(dayFilterLogical);
                    uniqueDayVesselEventTimes = singleVesselEventTimes(dayFilterLogical);
                    uniqueDayVesselData = singleVesselData(dayFilterLogical);
                    if ~isempty(uniqueDayVesselFileIDs)
                        uniqueDayFiles = unique(uniqueDayVesselFileIDs);
                        cutOffTime = targetMinutes/(trialDuration_sec/60);
                        dayLog = zeros(length(validDays),1);
                        for jj = 1:length(validDays)
                            for kk = 1:length(uniqueDayFiles)
                                if strcmp(validDays{jj,1},uniqueDayFiles{kk,1})
                                    dayLog(jj,1) = 1;
                                end
                            end
                        end
                        dayLog = logical(dayLog);
                        valDays = validDays(dayLog,:);
                        clear baselineFiles
                        try
                            baselineFiles = {valDays(1:cutOffTime,1)}';
                        catch
                            baselineFiles = {valDays(1:end,1)}';
                        end
                        timeFilterLogical = zeros(length(uniqueDayVesselFileIDs),1);
                        for ll = 1:length(uniqueDayVesselFileIDs)
                            uniqueDayVesselFileID = uniqueDayVesselFileIDs(ll,1);
                            for mm = 1:length(baselineFiles{1,1})
                                baselineFile = baselineFiles{1,1}{mm,1};
                                if strcmp(baselineFile, uniqueDayVesselFileID)
                                    timeFilterLogical(ll,1) = 1;
                                end
                            end
                        end
                        timeFilterLogical = logical(timeFilterLogical);
                        baselineVesselFileIDs = uniqueDayVesselFileIDs(timeFilterLogical);
                        baselineVesselDurations = uniqueDayVesselDurations(timeFilterLogical);
                        baselineVesselEventTimes = uniqueDayVesselEventTimes(timeFilterLogical);
                        baselineVesselData = uniqueDayVesselData(timeFilterLogical);
                        % find the means of each the unique vessel from the unique day
                        % for valid files
                        clear tempData_means
                        for nn = 1:length(baselineVesselData)
                            tempData_means(nn,1) = mean(baselineVesselData{nn,1});
                        end
                    else
                        tempData_means = [];
                        baselineVesselFileIDs = [];
                        baselineVesselDurations = [];
                        baselineVesselEventTimes = [];
                    end
                    RestingBaselines.setDuration.(dataType).(subDataType).(uV).(uD) = mean(tempData_means);
                end
                tempBaseFileIDs = vertcat(tempBaseFileIDs, baselineVesselFileIDs);
                tempBaseEventTimes = vertcat(tempBaseEventTimes, baselineVesselEventTimes);
                tempBaseDurations = vertcat(tempBaseDurations, baselineVesselDurations);
            end
        end
    else
        clear tempData_means
        subDataTypes = fieldnames(RestData.(dataType));   % Find the hemisphere dataTypes. These are typically LH, RH
        % Loop through each hemisphere dataType (LH, RH) because they are subfields and will have unique baselines
        for oo = 1:length(subDataTypes)
            subDataType = char(subDataTypes(oo));   % Load each loop iteration's hemisphere fieldname as a character string
            [restLogical] = FilterEvents_2P_eLife2025(RestData.(dataType).(subDataType),RestCriteria);
            [puffLogical] = FilterEvents_2P_eLife2025(RestData.(dataType).(subDataType),PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            allRestFiles = RestData.(dataType).(subDataType).fileIDs(combRestLogical,:);
            allRestDurations = RestData.(dataType).(subDataType).durations(combRestLogical,:);
            allRestEventTimes = RestData.(dataType).(subDataType).eventTimes(combRestLogical,:);
            restingData = RestData.(dataType).(subDataType).data(combRestLogical,:);
            uniqueDays = GetUniqueDays_2P_eLife2025(RestData.(dataType).(subDataType).fileIDs);   % Find the unique days of imaging
            uniqueFiles = unique(RestData.(dataType).(subDataType).fileIDs);   % Find the unique files from the filelist. This removes duplicates
            numberOfFiles = length(unique(RestData.(dataType).(subDataType).fileIDs));   % Find the number of unique files
            fileTarget = targetMinutes/(trialDuration_sec/60);   % Divide that number of unique files by 5 (minutes) to get the number of files
            % Loop through each unique day in order to create a logical to filter the file list so that it only includes the first
            % x number of files that fall within the targetMinutes requirement
            for pp = 1:length(uniqueDays)
                day = uniqueDays(pp);
                rr = 1;
                for qq = 1:numberOfFiles
                    file = uniqueFiles(qq);
                    fileID = file{1}(1:6);
                    if strcmp(day, fileID) && rr <= fileTarget
                        filtLogical{pp,1}(qq,1) = 1; %#ok<*AGROW>
                        rr = rr + 1;
                    else
                        filtLogical{pp,1}(qq,1) = 0;
                    end
                end
            end
            % Combine the 3 logicals so that it reflects the first "x" number of files from each day
            finalLogical = any(sum(cell2mat(filtLogical'),2),2);
            % Now that the appropriate files from each day are identified, loop through each file name with respect to the original
            % list of ALL resting files, only keeping the ones that fall within the first targetMinutes of each day.
            filtRestFiles = uniqueFiles(finalLogical,:);
            for ss = 1:length(allRestFiles)
                logic = strcmp(allRestFiles{ss},filtRestFiles);
                logicSum = sum(logic);
                if logicSum == 1
                    fileFilter(ss,1) = 1;
                else
                    fileFilter(ss,1) = 0;
                end
            end
            finalFileFilter = logical(fileFilter);
            finalFileIDs = allRestFiles(finalFileFilter,:);
            finalFileDurations = allRestDurations(finalFileFilter,:);
            finalFileEventTimes = allRestEventTimes(finalFileFilter,:);
            finalRestData = restingData(finalFileFilter,:);
            % Loop through each unique day and pull out the data that corresponds to the resting files
            for tt = 1:length(uniqueDays)
                vv = 1;
                for uu = 1:length(finalFileIDs)
                    fileID = finalFileIDs{uu,1}(1:6);
                    date{tt,1} = ConvertDate_2P_eLife2025(uniqueDays{tt,1});
                    if strcmp(fileID, uniqueDays{tt,1}) == 1
                        tempData.(date{tt,1}){vv,1} = finalRestData{uu,1};
                        vv = vv + 1;
                    end
                end
            end
            % find the means of each unique day
            for ww = 1:size(date,1)
                tempData_means{ww,1} = cellfun(@(x) mean(x),tempData.(date{ww,1}));    % LH date-specific means
            end
            % Save the means into the Baseline struct under the current loop iteration with the associated dates
            for xx = 1:length(uniqueDays)
                RestingBaselines.setDuration.(dataType).(subDataType).(date{xx,1}) = mean(tempData_means{xx,1});    % LH date-specific means
            end
        end
    end
end
RestingBaselines.setDuration.baselineFileInfo.fileIDs = finalFileIDs;
RestingBaselines.setDuration.baselineFileInfo.eventTimes = finalFileEventTimes;
RestingBaselines.setDuration.baselineFileInfo.durations = finalFileDurations;
RestingBaselines.setDuration.targetMinutes = targetMinutes;
save([animalID '_RestingBaselines.mat'],'RestingBaselines');

end
