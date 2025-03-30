function [] = AddSleepParameters_IOS_nNOS(procDataFileIDs)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% load the baseline structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)
binTime = 5; % seconds
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    [~,fileDate,~] = GetFileInfo_IOS_nNOS(procDataFileID);
    strDay = ConvertDate_IOS_nNOS(fileDate);
    load(procDataFileID)
    % if isfield(ProcData,'sleep') == false
    disp(['Adding sleep parameters to ProcData file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    specDataFileID = [procDataFileID(1:end - 12) 'SpecData.mat'];
    load(specDataFileID)
    [dataTypes] = DetermineWavelengthDatatypes_IOS_nNOS(ProcData.notes.imagingWavelengths,3);
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        if any(strcmp(dataType,{'HbT','HbO','HbR','GCaMP'})) == true
            samplingRate = ProcData.notes.wavelengthSamplingRate;
        elseif strcmp(dataType,'heartRate') == true
            samplingRate = 1;
        else
            samplingRate = ProcData.notes.dsFs;
        end
        subDataTypes = fieldnames(ProcData.data.(dataType));
        for cc = 1:length(subDataTypes)
            subDataType = subDataTypes{cc,1};
            if any(strcmp(dataType,{'cortical_LH','cortical_RH','hippocampus'})) == true
                data.(dataType).(subDataType) = (ProcData.data.(dataType).(subDataType) - RestingBaselines.manualSelection.(dataType).(subDataType).(strDay).mean)/RestingBaselines.manualSelection.(dataType).(subDataType).(strDay).mean;
            elseif strcmp(dataType,'EMG')
                data.(dataType).(subDataType) = (ProcData.data.(dataType).(subDataType) - RestingBaselines.manualSelection.(dataType).(subDataType).(strDay).mean);
            else
                data.(dataType).(subDataType) = ProcData.data.(dataType).(subDataType);
            end
            % divide the neural signals into five second bins and put them in a cell array
            tempData.(dataType).(subDataType).tempStruct = cell(ProcData.notes.trialDuration_sec/binTime,1);
            % loop through all samples across the 15 minutes in 5 second bins (180 total)
            for dd = 1:length(tempData.(dataType).(subDataType).tempStruct)
                if dd == 1
                    tempData.(dataType).(subDataType).tempStruct(dd,1) = {data.(dataType).(subDataType)(dd:samplingRate)};
                elseif dd == length(tempData.(dataType).(subDataType).tempStruct)
                    tempData.(dataType).(subDataType).tempStruct(dd,1) = {data.(dataType).(subDataType)((((samplingRate*(dd - 1)) + 1)):end)};
                else
                    tempData.(dataType).(subDataType).tempStruct(dd,1) = {data.(dataType).(subDataType)((((samplingRate*(dd - 1)) + 1)):(samplingRate*binTime*dd))};
                end
            end
            ProcData.sleep.parameters.(dataType).(subDataType) = tempData.(dataType).(subDataType).tempStruct;
        end
        specDataTypes = {'cortical_LH','cortical_RH','hippocampus'};
        LFP_bands = {'delta','theta','alpha','beta','gamma'};
        for ee = 1:length(specDataTypes)
            specDataType = specDataTypes{1,ee};
            offset = 2.5; % sec
            binWidth = 5; % sec
            T = round(SpecData.(specDataType).T,1);
            F = SpecData.(specDataType).F;
            specData = SpecData.(specDataType).normS;
            freqFloor = floor(F);
            for ff = 1:length(LFP_bands)
                LFP_band = LFP_bands{1,ff};
                switch LFP_band
                    case 'delta'
                        freqStartIdx = freqFloor == 1;
                        freqStopIdx = freqFloor == 4;
                    case 'theta'
                        freqStartIdx = freqFloor == 4;
                        freqStopIdx = freqFloor == 10;
                    case 'alpha'
                        freqStartIdx = freqFloor == 10;
                        freqStopIdx = freqFloor == 13;
                    case 'beta'
                        freqStartIdx = freqFloor == 13;
                        freqStopIdx = freqFloor == 30;
                    case 'gamma'
                        freqStartIdx = freqFloor == 30;
                        freqStopIdx = freqFloor == 100;
                end
                bandStart = find(freqStartIdx,1,'first');
                bandStop = find(freqStopIdx,1,'last');
                data.(specDataType).(LFP_band) = mean(specData(bandStart:bandStop,:),1);
                tempData.(specDataType).(LFP_band).tempStruct = cell(ProcData.notes.trialDuration_sec/binTime,1);
                for gg = 1:length(tempData.(specDataType).(LFP_band).tempStruct)
                    if gg == 1
                        startTime = offset;
                        startTimeIdx = find(T == startTime);
                        endTime = 5;
                        [~,endTimeIdx] = min(abs(T - endTime));
                        tempData.(specDataType).(LFP_band).tempStruct{gg,1} = {data.(specDataType).(LFP_band)(startTimeIdx:endTimeIdx)};
                    elseif dd == length(tempData.(specDataType).(LFP_band).tempStruct)
                        startTime = ProcData.notes.trialDuration_sec - 5;
                        [~,startTimeIdx] = min(abs(T - startTime));
                        endTime = ProcData.notes.trialDuration_sec - offset;
                        [~,endTimeIdx] = min(abs(T - endTime));
                        tempData.(specDataType).(LFP_band).tempStruct{gg,1} = {data.(specDataType).(LFP_band)(startTimeIdx:endTimeIdx)};
                    else
                        startTime = binWidth*(dd - 1);
                        [~,startTimeIdx] = min(abs(T - startTime));
                        endTime = binWidth*dd;
                        [~,endTimeIdx] = min(abs(T - endTime));
                        tempData.(specDataType).(LFP_band).tempStruct{gg,1} = {data.(specDataType).(LFP_band)(startTimeIdx + 1:endTimeIdx + 1)};
                    end
                end
                ProcData.sleep.parameters.(specDataType).(LFP_band) = tempData.(specDataType).(LFP_band).tempStruct;
            end
        end
    end
    % save data structure
    save(procDataFileID,'ProcData');
    % else
    %     disp(['Sleep parameters for file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ') already added.']); disp(' ')
    % end
end