function [Results_Transitions_GCaMP] = AnalyzeArousalTransitions_GCaMP_eLife2025(animalID,group,set,rootFolder,delim,Results_Transitions_GCaMP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
modelName = [animalID '_IOS_eLife2025_RF_SleepScoringModel.mat'];
load(modelName)
modelDataFileStruct = dir('*_ModelData.mat');
modelDataFile = {modelDataFileStruct.name}';
modelDataFileIDs = char(modelDataFile);
baselineFileStruct = dir('*_RestingBaselines.mat');
baselineFile = {baselineFileStruct.name}';
baselineFileID = char(baselineFile);
load(baselineFileID)
samplingRate = 30;
iosFs = 10;
specSamplingRate = 10;
hemispheres = {'LH','RH','fLH','fRH'};
hemoDataTypes = {'HbT','HbO','HbR','GCaMP'};
fileDates = fieldnames(RestingBaselines.manualSelection.CBV.LH);
% go through each file and sleep score the data
for aa = 1:size(modelDataFileIDs,1)
    modelDataFileID = modelDataFileIDs(aa,:);
    if aa == 1
        load(modelDataFileID)
        dataLength = size(paramsTable,1);
        joinedTable = paramsTable;
        joinedFileList = cell(size(paramsTable,1),1);
        joinedFileList(:) = {modelDataFileID};
    else
        load(modelDataFileID)
        fileIDCells = cell(size(paramsTable,1),1);
        fileIDCells(:) = {modelDataFileID};
        joinedTable = vertcat(joinedTable,paramsTable);
        joinedFileList = vertcat(joinedFileList,fileIDCells);
    end
end
scoringTable = joinedTable;
[labels,~] = predict(RF_MDL,scoringTable);
% apply a logical patch on the REM events
REMindex = strcmp(labels,'REM Sleep');
numFiles = length(labels)/dataLength;
reshapedREMindex = reshape(REMindex,dataLength,numFiles);
patchedREMindex = [];
% patch missing REM indeces due to theta band falling off
for bb = 1:size(reshapedREMindex,2)
    remArray = reshapedREMindex(:,bb);
    patchedREMarray = LinkBinaryEvents_IOS_eLife2025(remArray',[5,0]);
    patchedREMindex = vertcat(patchedREMindex,patchedREMarray');
end
% change labels for each event
for cc = 1:length(labels)
    if patchedREMindex(cc,1) == 1
        labels{cc,1} = 'REM Sleep';
    end
end
% convert strings to numbers for easier comparisons
labelNumbers = zeros(length(labels),1);
for dd = 1:length(labels)
    if strcmp(labels{dd,1},'Not Sleep') == true
        labelNumbers(dd,1) = 1;
    elseif strcmp(labels{dd,1},'NREM Sleep') == true
        labelNumbers(dd,1) = 2;
    elseif strcmp(labels{dd,1},'REM Sleep') == true
        labelNumbers(dd,1) = 3;
    end
end
% reshape
fileIDs = unique(joinedFileList);
fileLabels = reshape(labelNumbers,dataLength,size(modelDataFileIDs,1))';
stringArray = zeros(1,12);
for dd = 1:length(transitions)
    transition = transitions{1,dd};
    if strcmp(transition,'AWAKEtoNREM') == true
        for ee = 1:12
            if ee <= 6
                stringArray(1,ee) = 1;
            else
                stringArray(1,ee) = 2;
            end
        end
    elseif strcmp(transition,'NREMtoAWAKE') == true
        for ee = 1:12
            if ee <= 6
                stringArray(1,ee) = 2;
            else
                stringArray(1,ee) = 1;
            end
        end
    elseif strcmp(transition,'NREMtoREM') == true
        for ee = 1:12
            if ee <= 6
                stringArray(1,ee) = 2;
            else
                stringArray(1,ee) = 3;
            end
        end
    elseif strcmp(transition,'REMtoAWAKE') == true
        for ee = 1:12
            if ee <= 6
                stringArray(1,ee) = 3;
            else
                stringArray(1,ee) = 1;
            end
        end
    end
    % go through and pull out all indeces of matching strings
    idx = 1;
    for f = 1:length(fileIDs)
        fileID = fileIDs{f,1};
        labelArray = fileLabels(f,:);
        indeces = strfind(labelArray,stringArray);
        if isempty(indeces) == false
            for g1 = 1:length(indeces)
                data.(transition).files{idx,1} = fileID;
                data.(transition).startInd(idx,1) = indeces(1,g1);
                idx = idx + 1;
            end
        end
    end
end
% extract data
for hh = 1:length(transitions)
    transition = transitions{1,hh};
    iqx = 1;
    for ii = 1:length(data.(transition).files)
        file = data.(transition).files{ii,1};
        startBin = data.(transition).startInd(ii,1);
        if startBin > 1 && startBin < (180 - 12)
            [animalID,fileDate,fileID] = GetFileInfo_IOS_eLife2025(file);
            strDay = ConvertDate_IOS_eLife2025(fileDate);
            procDataFileID = [animalID '_' fileID '_ProcData.mat'];
            load(procDataFileID)
            specDataFileID = [animalID '_' fileID '_SpecDataB.mat'];
            load(specDataFileID)
            startTime = (startBin - 1)*5; % sec
            endTime = startTime + (12*5); % sec
            % whisking data
            [z1,p1,k1] = butter(4,10/(samplingRate/2),'low');
            [sos1,g1] = zp2sos(z1,p1,k1);
            filtWhiskAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle(startTime*samplingRate + 1:endTime*samplingRate));
            % heart rate data
            heartRate = ProcData.data.heartRate(startTime + 1:endTime);
            % EMG
            EMG = (ProcData.data.EMG.emg(startTime*samplingRate + 1:endTime*samplingRate) - RestingBaselines.manualSelection.EMG.emg.(strDay).mean);
            % spectrogram data
            cortical_LHnormS = SpecData.cortical_LH.normS;
            hippocampusNormS = SpecData.hippocampus.normS;
            T = round(SpecData.cortical_LH.T,1);
            F = SpecData.cortical_LH.F;
            specStartIndex = find(T == startTime);
            specStartIndex = specStartIndex(1);
            specEndIndex = find(T == endTime);
            specEndIndex = specEndIndex(end);
            LH_cortSpec = cortical_LHnormS(:,specStartIndex + 1:specEndIndex);
            Hip_spec = hippocampusNormS(:,specStartIndex + 1:specEndIndex);
            T_short = T(1:size(LH_cortSpec,2));
            [z2,p2,k2] = butter(4,1/(iosFs/2),'low');
            [sos2,g2] = zp2sos(z2,p2,k2);
            for qq = 1:length(hemispheres)
                hemisphere = hemispheres{1,qq};
                for zz = 1:length(hemoDataTypes)
                    hemoDataType = hemoDataTypes{1,zz};
                    tempData = ProcData.data.(hemoDataType).(hemisphere);
                    filtData.(hemisphere).(hemoDataType) = filtfilt(sos2,g2,tempData(startTime*iosFs + 1:endTime*iosFs));
                end
            end
            % set transition
            data.(transition).fileDate{iqx,1} = strDay;
            data.(transition).whisk(iqx,:) = filtWhiskAngle;
            data.(transition).HR(iqx,:) = heartRate;
            data.(transition).EMG(iqx,:) = EMG;
            data.(transition).LH_cort(:,:,iqx) = LH_cortSpec(:,1:specSamplingRate*60);
            data.(transition).Hip(:,:,iqx) = Hip_spec(:,1:specSamplingRate*60);
            data.(transition).T_short = T_short(1:specSamplingRate*60);
            data.(transition).F = F;
            for qq = 1:length(hemispheres)
                hemisphere = hemispheres{1,qq};
                for zz = 1:length(hemoDataTypes)
                    hemoDataType = hemoDataTypes{1,zz};
                    data.(hemisphere).(transition).(hemoDataType)(iqx,:) = filtData.(hemisphere).(hemoDataType);
                end
            end
            iqx = iqx + 1;
        end
    end
end
for qq = 1:length(hemispheres)
    hemisphere = hemispheres{1,qq};
    % take averages of each behavior
    for dd = 1:length(transitions)
        transition = transitions{1,dd};
        % save results
        Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).whisk = mean(data.(transition).whisk,1);
        Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).HR = mean(data.(transition).HR,1);
        Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).EMG = mean(data.(transition).EMG,1);
        Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).Hip = mean(data.(transition).Hip,3);
        Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).T = data.(transition).T_short;
        Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).F = data.(transition).F;
        Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).indFileDate = data.(transition).fileDate;
        Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).fileDates = fileDates;
        Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).Cort = mean(data.(transition).LH_cort,3);
        for zz = 1:length(hemoDataTypes)
            hemoDataType = hemoDataTypes{1,zz};
            Results_Transitions_GCaMP.(group).(animalID).(hemisphere).(transition).(hemoDataType) = mean(data.(hemisphere).(transition).(hemoDataType),1);
        end
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_Transitions_GCaMP.mat','Results_Transitions_GCaMP')
cd([rootFolder delim 'Data'])