function [RestingBaselines] = CalculateManualRestingBaselinesTimeIndeces_2P_eLife2025
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Manually designate files with event times that correspond to appropriate rest
%________________________________________________________________________________________________________________________

disp('Calculating the resting baselines using manually selected files from each unique day...'); disp(' ')
% character list of all MergedData files
mergedDataFileStruct = dir('*_MergedData.mat');
mergedDataFiles = {mergedDataFileStruct.name}';
mergedDataFileIDs = char(mergedDataFiles);
% find and load RestingBaselines.mat struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFiles = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFiles);
load(baselineDataFileID)
% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFiles = {restDataFileStruct.name}';
restDataFileID = char(restDataFiles);
load(restDataFileID)
% determine the animal's ID use the RestData.mat file's name for the current folder
fileBreaks = strfind(restDataFileID,'_');
animalID = restDataFileID(1:fileBreaks(1)-1);
% the RestData.mat struct has all resting events, regardless of duration. We want to set the threshold for rest as anything
% that is greater than a certain amount of time
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {5};
PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};
% loop through each file and manually designate which files have appropriate amounts of rest
% if this is already completed, load the struct and skip
for aa = 1:size(mergedDataFileIDs,1)
    disp(['Loading file ' num2str(aa) ' of ' num2str(size(mergedDataFileIDs,1)) '...']); disp(' ')
    mergedDataFileID = mergedDataFileIDs(aa,:);
    load(mergedDataFileID,'-mat')
    ManualDecisions.fileIDs{aa,1} = mergedDataFileIDs(aa,:);
    if isfield(MergedData,'manualBaselineInfo') == false
        saveFigs = 'n';
        baselineType = 'setDuration';
        bb = false;
        while bb == false
            % load a figure with the data to visualize which periods are rest. Note that this data is, by default, normalized
            % by the first 30 minutes of data which may or may not reflect accurate normalizations
            [singleTrialFig] = GenerateSingleFigures_2P_eLife2025(mergedDataFileID,baselineType,saveFigs,RestingBaselines);
            fileDecision = input(['Use data from ' mergedDataFileID ' for resting baseline calculation? (y/n): '], 's'); disp(' ')
            if strcmp(fileDecision,'y') || strcmp(fileDecision,'n')
                bb = true;
                ManualDecisions.validFiles{aa,1} = fileDecision;
                if strcmp(fileDecision,'y') == true
                    startTime = input('Input the start time for resting data: '); disp(' ')
                    endTime = input('Input the end time for resting data: '); disp(' ')
                    ManualDecisions.startTimes{aa,1} = startTime;
                    ManualDecisions.endTimes{aa,1} = endTime;
                    if isempty(startTime) == true || isempty(endTime) == true
                        b = false;
                    end
                else
                    startTime = 0;
                    endTime = 0;
                    ManualDecisions.startTimes{aa,1} = 0;
                    ManualDecisions.endTimes{aa,1} = 0;
                end
                close(singleTrialFig)
                load(mergedDataFileID)
                MergedData.manualBaselineInfo.fileDecision = fileDecision;
                MergedData.manualBaselineInfo.startTime = startTime;
                MergedData.manualBaselineInfo.endTime = endTime;
                save(mergedDataFileID,'MergedData')
            else
                bb = false;
                close(singleTrialFig)
            end
        end
    else
        ManualDecisions.validFiles{aa,1} = MergedData.manualBaselineInfo.fileDecision;
        ManualDecisions.startTimes{aa,1} = MergedData.manualBaselineInfo.startTime;
        ManualDecisions.endTimes{aa,1} = MergedData.manualBaselineInfo.endTime;
    end
end
% save structure with each file's decision and start/end time
save([animalID '_ManualBaselineFileList.mat'],'ManualDecisions')

% decimate the unique file list to identify which files we decided to use for resting time
% pull the start/end times along with each desired file
cc = 1;
for dd = 1:length(ManualDecisions.validFiles)
    if strcmp(ManualDecisions.validFiles{dd,1},'y') == true
        fileBreaks = strfind(ManualDecisions.fileIDs{dd,1},'_');
        subFilterFileList_A{cc,1} = ManualDecisions.fileIDs{dd,1}(fileBreaks(2) + 1:fileBreaks(6) - 1); %#ok<*AGROW>
        subFilterFileList_B{cc,1} = ManualDecisions.fileIDs{dd,1};
        subFilterStartTimes{cc,1} = ManualDecisions.startTimes{dd,1};
        subFilterEndTimes{cc,1} = ManualDecisions.endTimes{dd,1};
        cc = cc + 1;
    end
end
% find the fieldnames of RestData and loop through each field. Each fieldname should be a different dataType of interest.
% these will typically be CBV, Delta, Theta, Gamma, and MUA, etc
dataTypes = fieldnames(RestData);
for ee = 1:length(dataTypes)
    dataType = char(dataTypes(ee));
    clear tempData tempDataMeans
    if strcmp(dataType,'vesselDiameter') == false
        % find any sub-dataTypes. These are typically LH, RH
        subDataTypes = fieldnames(RestData.(dataType));
        for ff = 1:length(subDataTypes)
            subDataType = char(subDataTypes(ff));
            % use the criteria we specified earlier to find all resting events that are greater than the criteria
            [restLogical] = FilterEvents_2P_eLife2025(RestData.(dataType).(subDataType),RestCriteria);
            [puffLogical] = FilterEvents_2P_eLife2025(RestData.(dataType).(subDataType),PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            allRestFileIDs = RestData.(dataType).(subDataType).fileIDs(combRestLogical,:);
            allRestDurations = RestData.(dataType).(subDataType).durations(combRestLogical,:);
            allRestEventTimes = RestData.(dataType).(subDataType).eventTimes(combRestLogical,:);
            allRestingData = RestData.(dataType).(subDataType).data(combRestLogical,:);
            % find the unique days and unique file IDs
            uniqueDays = GetUniqueDays_2P_eLife2025(RestData.(dataType).(subDataType).fileIDs);
            uniqueFiles = unique(RestData.(dataType).(subDataType).fileIDs);
            numberOfFiles = length(unique(RestData.(dataType).(subDataType).fileIDs));
            % loop through each unique day in order to create a logical to filter the file list
            for gg = 1:length(uniqueDays)
                uniqueDay = uniqueDays(gg);
                hh = 1;
                for ii = 1:numberOfFiles
                    uniqueFileID = uniqueFiles(ii);
                    uniqueFileID_short = uniqueFileID{1}(1:6);
                    goodFile = 'n';
                    % determine if the file occurs during this specific day - this is for all files
                    for jj = 1:length(subFilterFileList_A)
                        subFilterFile = subFilterFileList_A{jj,1};
                        if strcmp(uniqueFileID,subFilterFile) == true
                            goodFile = 'y';
                        end
                    end
                    % determine whether the approved files are part of the 'approved' file list
                    if strcmp(uniqueDay,uniqueFileID_short) == true && strcmp(goodFile,'y') == true
                        uniqueDayFiltLogical{gg,1}(ii,1) = 1;
                        hh = hh + 1;
                    else
                        uniqueDayFiltLogical{gg,1}(ii,1) = 0;
                    end
                end
            end
            finalUniqueDayFiltLogical = any(sum(cell2mat(uniqueDayFiltLogical'),2),2);
            % now that the appropriate files from each day are identified, loop through each file name with respect to the original
            % list of ALL resting files, only keeping the ones that fall within the first targetMinutes of each day.
            filtRestFiles = uniqueFiles(finalUniqueDayFiltLogical,:);
            for kk = 1:length(allRestFileIDs)
                fileCompare = strcmp(allRestFileIDs{kk},filtRestFiles);
                includeFile = sum(fileCompare);
                if includeFile == 1
                    allFileFilter(kk,1) = 1;
                else
                    allFileFilter(kk,1) = 0;
                end
            end
            AllFileFilter = logical(allFileFilter);
            filtFileIDs = allRestFileIDs(AllFileFilter,:);
            filtDurations = allRestDurations(AllFileFilter,:);
            filtEventTimes = allRestEventTimes(AllFileFilter,:);
            filtRestData = allRestingData(AllFileFilter,:);
            % now that we have decimated the original list to only reflect the proper unique day, approved files
            % we want to only take events that occur during our approved time duration
            for ll = 1:length(filtFileIDs)
                finalFileID = filtFileIDs{ll,1};
                for mm = 1:length(subFilterFileList_A)
                    sFile = subFilterFileList_A{mm,1};
                    if strcmp(finalFileID,sFile) == true
                        sTime = subFilterStartTimes{mm,1};
                        eTime = subFilterEndTimes{mm,1};
                    end
                end
                eventTime = filtEventTimes(ll,1);
                % 2 seconds before event time, 10 seconds after
                if (eventTime - 2) >= sTime && (eventTime + 10) <= eTime
                    eventTimeFilter(ll,1) = 1;
                else
                    eventTimeFilter(ll,1) = 0;
                end
            end
            EventTimeFilter = logical(eventTimeFilter);
            finalEventFileIDs = filtFileIDs(EventTimeFilter,:);
            finalEventDurations = filtDurations(EventTimeFilter,:);
            finalEventTimes = filtEventTimes(EventTimeFilter,:);
            finalEventRestData = filtRestData(EventTimeFilter,:);
            % again loop through each unique day and pull out the data that corresponds to the final resting files
            for nn = 1:length(uniqueDays)
                oo = 1;
                for pp = 1:length(finalEventFileIDs)
                    uniqueFileID_short = finalEventFileIDs{pp,1}(1:6);
                    uniqueDate{nn,1} = ConvertDate_2P_eLife2025(uniqueDays{nn,1});
                    if strcmp(uniqueFileID_short,uniqueDays{nn,1}) == 1
                        tempData.(uniqueDate{nn,1}){oo,1} = finalEventRestData{pp,1};
                        oo = oo + 1;
                    end
                end
            end
            % find the means of each unique day
            for qq = 1:size(uniqueDate,1)
                tempDataMeans{qq,1} = cellfun(@(x)mean(x),tempData.(uniqueDate{qq,1}));
            end
            % save the means into the Baseline struct under the current loop iteration with the associated dates
            for rr = 1:length(uniqueDays)
                RestingBaselines.manualSelection.(dataType).(subDataType).(uniqueDate{rr,1}) = mean(tempDataMeans{rr,1});
            end
        end
    else
        % find any sub-dataTypes. These are typically LH, RH
        subDataTypes = fieldnames(RestData.(dataType));
        for ss = 1:length(subDataTypes)
            subDataType = char(subDataTypes(ss));
            % use the criteria we specified earlier to find all resting events that are greater than the criteria
            [restLogical] = FilterEvents_2P_eLife2025(RestData.(dataType).(subDataType),RestCriteria);
            [puffLogical] = FilterEvents_2P_eLife2025(RestData.(dataType).(subDataType),PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            allRestFileIDs = RestData.(dataType).(subDataType).fileIDs(combRestLogical,:);
            allRestDurations = RestData.(dataType).(subDataType).durations(combRestLogical,:);
            allRestEventTimes = RestData.(dataType).(subDataType).eventTimes(combRestLogical,:);
            allRestingData = RestData.(dataType).(subDataType).data(combRestLogical,:);
            % find the unique days and unique file IDs
            uniqueDays = GetUniqueDays_2P_eLife2025(RestData.(dataType).(subDataType).fileIDs);
            uniqueFiles = unique(RestData.(dataType).(subDataType).fileIDs);
            numberOfFiles = length(unique(RestData.(dataType).(subDataType).fileIDs));
            % loop through each unique day in order to create a logical to filter the file list
            for tt = 1:length(uniqueDays)
                uniqueDay = uniqueDays(tt);
                uu = 1;
                for vv = 1:numberOfFiles
                    uniqueFileID = uniqueFiles(vv);
                    uniqueFileID_short = uniqueFileID{1}(1:6);
                    goodFile = 'n';
                    % determine if the file occurs during this specific day - this is for all files
                    for ww = 1:length(subFilterFileList_A)
                        subFilterFile = subFilterFileList_A{ww,1};
                        if strcmp(uniqueFileID,subFilterFile) == true
                            goodFile = 'y';
                        end
                    end
                    % determine whether the approved files are part of the 'approved' file list
                    if strcmp(uniqueDay,uniqueFileID_short) == true && strcmp(goodFile,'y') == true
                        uniqueDayFiltLogical{tt,1}(vv,1) = 1;
                        uu = uu + 1;
                    else
                        uniqueDayFiltLogical{tt,1}(vv,1) = 0;
                    end
                end
            end
            finalUniqueDayFiltLogical = any(sum(cell2mat(uniqueDayFiltLogical'),2),2);
            % now that the appropriate files from each day are identified, loop through each file name with respect to the original
            % list of ALL resting files, only keeping the ones that fall within the first targetMinutes of each day.
            filtRestFiles = uniqueFiles(finalUniqueDayFiltLogical,:);
            for xx = 1:length(allRestFileIDs)
                fileCompare = strcmp(allRestFileIDs{xx},filtRestFiles);
                includeFile = sum(fileCompare);
                if includeFile == 1
                    allFileFilter(xx,1) = 1;
                else
                    allFileFilter(xx,1) = 0;
                end
            end
            AllFileFilter = logical(allFileFilter);
            filtFileIDs = allRestFileIDs(AllFileFilter,:);
            filtDurations = allRestDurations(AllFileFilter,:);
            filtEventTimes = allRestEventTimes(AllFileFilter,:);
            filtRestData = allRestingData(AllFileFilter,:);
            % now that we have decimated the original list to only reflect the proper unique day, approved files
            % we want to only take events that occur during our approved time duration
            for yy = 1:length(filtFileIDs)
                finalFileID = filtFileIDs{yy,1};
                for zz = 1:length(subFilterFileList_A)
                    sFile = subFilterFileList_A{zz,1};
                    if strcmp(finalFileID,sFile) == true
                        sTime = subFilterStartTimes{zz,1};
                        eTime = subFilterEndTimes{zz,1};
                    end
                end
                eventTime = filtEventTimes(yy,1);
                % 2 seconds before event time, 10 seconds after
                if (eventTime - 2) >= sTime && (eventTime + 10) <= eTime
                    eventTimeFilter(yy,1) = 1;
                else
                    eventTimeFilter(yy,1) = 0;
                end
            end
            EventTimeFilter = logical(eventTimeFilter);
            finalEventFileIDs = filtFileIDs(EventTimeFilter,:);
            finalEventDurations = filtDurations(EventTimeFilter,:);
            finalEventTimes = filtEventTimes(EventTimeFilter,:);
            finalEventRestData = filtRestData(EventTimeFilter,:);
            % again loop through each unique day and pull out the data that corresponds to the final resting files
            for aaa = 1:length(uniqueDays)
                ccc = 1;
                for bbb = 1:length(finalEventFileIDs)
                    uniqueFileID_short = finalEventFileIDs{bbb,1}(1:6);
                    uniqueDate{aaa,1} = ConvertDate_2P_eLife2025(uniqueDays{aaa,1});
                    if strcmp(uniqueFileID_short,uniqueDays{aaa,1}) == 1
                        tempData.(uniqueDate{aaa,1}).data{ccc,1} = finalEventRestData{bbb,1};
                        tempData.(uniqueDate{aaa,1}).fileIDs{ccc,1} = finalEventFileIDs{bbb,1};
                        ccc = ccc + 1;
                    end
                end
            end
            % determine the vessel IDs per date
            for ddd = 1:length(uniqueDays)
                uniqueDate{ddd,1} = ConvertDate_2P_eLife2025(uniqueDays{ddd,1});
                for eee = 1:length(subFilterFileList_B)
                    [~,~,fileDate,~,~,~] = GetFileInfo2_2P_eLife2025(subFilterFileList_B{eee,1});
                    if strcmp(fileDate,uniqueDays{ddd,1}) == true
                        fileListFilter(eee,1) = 1;
                    else
                        fileListFilter(eee,1) = 0;
                    end
                end
                fileListFilter = logical(fileListFilter);
                validFileList = subFilterFileList_B(fileListFilter,:);
                for fff = 1:length(validFileList)
                    [~,~,~,vesselFileIDs{fff,1},~,vesselIDs{fff,1}] = GetFileInfo2_2P_eLife2025(validFileList{fff,1});
                end
                tempData.(uniqueDate{ddd,1}).vesselIDs = vesselIDs;
                tempData.(uniqueDate{ddd,1}).vesselFileIDs = vesselFileIDs;
            end
            %
            for ggg = 1:length(uniqueDays)
                uniqueDate{ggg,1} = ConvertDate_2P_eLife2025(uniqueDays{ggg,1});
                for hhh = 1:length(tempData.(uniqueDate{ggg,1}).vesselIDs)
                    vID = tempData.(uniqueDate{ggg,1}).vesselIDs{hhh,1};
                    vFile = tempData.(uniqueDate{ggg,1}).vesselFileIDs{hhh,1};
                    if isfield(tempData.(uniqueDate{ggg,1}),vID) == false
                        jjj = 1;
                        for iii = 1:length(tempData.(uniqueDate{ggg,1}).fileIDs)
                            if strcmp(vFile,tempData.(uniqueDate{ggg,1}).fileIDs{iii,1}) == true
                                tempData.(uniqueDate{ggg,1}).(vID).data{jjj,1} = tempData.(uniqueDate{ggg,1}).data{iii,1};
                                jjj = jjj + 1;
                            end
                        end
                    else
                        lll = length(tempData.(uniqueDate{ggg,1}).(vID).data) + 1;
                        for mmm = 1:length(tempData.(uniqueDate{ggg,1}).fileIDs)
                            if strcmp(vFile,tempData.(uniqueDate{ggg,1}).fileIDs{mmm,1}) == true
                                tempData.(uniqueDate{ggg,1}).(vID).data{lll,1} = tempData.(uniqueDate{ggg,1}).data{mmm,1};
                                lll = lll + 1;
                            end
                        end
                    end
                end
            end
            %
            for nnn = 1:length(subFilterFileList_B)
                [~,~,~,~,~,vessels{nnn,1}] = GetFileInfo2_2P_eLife2025(subFilterFileList_B{nnn,1});
            end
            vIDs = unique(vessels);
            % find the means of each unique day
            for ooo = 1:size(uniqueDate,1)
                fields = fieldnames(tempData.(uniqueDate{ooo,1}));
                for ppp = 1:length(fields)
                    field = fields{ppp,1};
                    if ~strcmp(field,'data') && ~strcmp(field,'fileIDs') && ~strcmp(field,'vesselIDs') && ~strcmp(field,'vesselFileIDs')
                        tempDataMeans.(field).(uniqueDate{ooo,1}) = cellfun(@(x)mean(x),tempData.(uniqueDate{ooo,1}).(field).data);
                    end
                end
            end
            % save the means into the Baseline struct under the current loop iteration with the associated dates
            for qqq = 1:length(vIDs)
                dates = fieldnames(tempDataMeans.(vIDs{qqq,1}));
                for rrr = 1:length(dates)
                    date = dates{rrr,1};
                    RestingBaselines.manualSelection.(dataType).(subDataType).(vIDs{qqq,1}).(date) = mean(tempDataMeans.(vIDs{qqq,1}).(date));
                end
            end
        end
    end
end
% save results
RestingBaselines.manualSelection.baselineFileInfo.fileIDs = finalEventFileIDs;
RestingBaselines.manualSelection.baselineFileInfo.eventTimes = finalEventTimes;
RestingBaselines.manualSelection.baselineFileInfo.durations = finalEventDurations;
RestingBaselines.manualSelection.baselineFileInfo.selections = ManualDecisions.validFiles;
RestingBaselines.manualSelection.baselineFileInfo.selectionFiles = ManualDecisions.fileIDs;
save(baselineDataFileID,'RestingBaselines')

end
