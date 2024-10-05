function [SleepData] = CreateSleepData_2P(NREMsleepTime,REMsleepTime,modelName,SleepData)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% character list of all MergedData files in current directory
mergedDataFileStruct = dir('*_MergedData.mat');
mergedDataFiles = {mergedDataFileStruct.name}';
mergedDataFileIDs = char(mergedDataFiles);
%% create sleep scored data structure for NREM data
% Identify sleep epochs and place in SleepData.mat structure
sleepBins = NREMsleepTime/5;
for aa = 1:size(mergedDataFileIDs,1)
    clearvars -except aa mergedDataFileIDs sleepBins NREMsleepTime REMsleepTime modelName SleepData
    mergedDataFileID = mergedDataFileIDs(aa,:);
    load(mergedDataFileID);
    [~,~,~,fileID,~,vesselID] = GetFileInfo2_2P(mergedDataFileID);
    nremLogical = MergedData.sleep.logicals.(modelName).nremLogical;
    targetTime = ones(1,sleepBins);
    % find the periods of time where there are the desired number of consecutive bins
    sleepIndex = find(conv(nremLogical,targetTime) >= sleepBins) - (sleepBins - 1);
    if ~isempty(sleepIndex)
        sleepCriteria = (0:(sleepBins - 1));
        % sleep Index now has the proper time stamps from sleep logical
        fixedSleepIndex = unique(sleepIndex + sleepCriteria);
        % loop through the length of sleep Index, and pull out associated data
        for bb = 1:length(fixedSleepIndex)
            % cortical neural bands
            Cort_deltaPower{bb,1} = MergedData.sleep.parameters.corticalNeural.deltaBandPower{fixedSleepIndex(bb),1}; %#ok<*AGROW>
            Cort_thetaPower{bb,1} = MergedData.sleep.parameters.corticalNeural.thetaBandPower{fixedSleepIndex(bb),1};
            Cort_alphaPower{bb,1} = MergedData.sleep.parameters.corticalNeural.alphaBandPower{fixedSleepIndex(bb),1};
            Cort_betaPower{bb,1} = MergedData.sleep.parameters.corticalNeural.betaBandPower{fixedSleepIndex(bb),1};
            Cort_gammaPower{bb,1} = MergedData.sleep.parameters.corticalNeural.gammaBandPower{fixedSleepIndex(bb),1};
            Cort_muaPower{bb,1} = MergedData.sleep.parameters.corticalNeural.muaPower{fixedSleepIndex(bb),1};
            % hippocampal neural bands
            Hip_deltaPower{bb,1} = MergedData.sleep.parameters.hippocampalNeural.deltaBandPower{fixedSleepIndex(bb),1};
            Hip_thetaPower{bb,1} = MergedData.sleep.parameters.hippocampalNeural.thetaBandPower{fixedSleepIndex(bb),1};
            Hip_alphaPower{bb,1} = MergedData.sleep.parameters.hippocampalNeural.alphaBandPower{fixedSleepIndex(bb),1};
            Hip_betaPower{bb,1} = MergedData.sleep.parameters.hippocampalNeural.betaBandPower{fixedSleepIndex(bb),1};
            Hip_gammaPower{bb,1} = MergedData.sleep.parameters.hippocampalNeural.gammaBandPower{fixedSleepIndex(bb),1};
            Hip_muaPower{bb,1} = MergedData.sleep.parameters.hippocampalNeural.muaPower{fixedSleepIndex(bb),1};
            % vessel diameter
            vesselDiameter{bb,1} = MergedData.sleep.parameters.vesselDiameter.data{fixedSleepIndex(bb),1};
            % whisker angle
            WhiskerAcceleration{bb,1} = MergedData.sleep.parameters.whiskerAcceleration{fixedSleepIndex(bb),1};
            % bin times
            BinTimes{bb,1} = 5*fixedSleepIndex(bb);
        end
        indexBreaks = find(fixedSleepIndex(2:end) - fixedSleepIndex(1:end - 1) > 1);    % Find if there are numerous sleep periods
        % if there is only one period of sleep in this file and not multiple events
        if isempty(indexBreaks)
            % cortical neural bands
            matCort_DeltaPower = cell2mat(Cort_deltaPower);
            arrayCort_DeltaPower = reshape(matCort_DeltaPower',[1,size(matCort_DeltaPower,2)*size(matCort_DeltaPower,1)]);
            cellCort_DeltaPower = {arrayCort_DeltaPower};
            matCort_ThetaPower = cell2mat(Cort_thetaPower);
            arrayCort_ThetaPower = reshape(matCort_ThetaPower',[1,size(matCort_ThetaPower,2)*size(matCort_ThetaPower,1)]);
            cellCort_ThetaPower = {arrayCort_ThetaPower};
            matCort_AlphaPower = cell2mat(Cort_alphaPower);
            arrayCort_AlphaPower = reshape(matCort_AlphaPower',[1,size(matCort_AlphaPower,2)*size(matCort_AlphaPower,1)]);
            cellCort_AlphaPower = {arrayCort_AlphaPower};
            matCort_BetaPower = cell2mat(Cort_betaPower);
            arrayCort_BetaPower = reshape(matCort_BetaPower',[1,size(matCort_BetaPower,2)*size(matCort_BetaPower,1)]);
            cellCort_BetaPower = {arrayCort_BetaPower};
            matCort_GammaPower = cell2mat(Cort_gammaPower);
            arrayCort_GammaPower = reshape(matCort_GammaPower',[1,size(matCort_GammaPower,2)*size(matCort_GammaPower,1)]);
            cellCort_GammaPower = {arrayCort_GammaPower};
            matCort_MUAPower = cell2mat(Cort_muaPower);
            arrayCort_MUAPower = reshape(matCort_MUAPower',[1,size(matCort_MUAPower,2)*size(matCort_MUAPower,1)]);
            cellCort_MUAPower = {arrayCort_MUAPower};
            % hippocampal neural bands
            matHip_DeltaPower = cell2mat(Hip_deltaPower);
            arrayHip_DeltaPower = reshape(matHip_DeltaPower',[1,size(matHip_DeltaPower,2)*size(matHip_DeltaPower,1)]);
            cellHip_DeltaPower = {arrayHip_DeltaPower};
            matHip_ThetaPower = cell2mat(Hip_thetaPower);
            arrayHip_ThetaPower = reshape(matHip_ThetaPower',[1,size(matHip_ThetaPower,2)*size(matHip_ThetaPower,1)]);
            cellHip_ThetaPower = {arrayHip_ThetaPower};
            matHip_AlphaPower = cell2mat(Hip_alphaPower);
            arrayHip_AlphaPower = reshape(matHip_AlphaPower',[1,size(matHip_AlphaPower,2)*size(matHip_AlphaPower,1)]);
            cellHip_AlphaPower = {arrayHip_AlphaPower};
            matHip_BetaPower = cell2mat(Hip_betaPower);
            arrayHip_BetaPower = reshape(matHip_BetaPower',[1,size(matHip_BetaPower,2)*size(matHip_BetaPower,1)]);
            cellHip_BetaPower = {arrayHip_BetaPower};
            matHip_GammaPower = cell2mat(Hip_gammaPower);
            arrayHip_GammaPower = reshape(matHip_GammaPower',[1,size(matHip_GammaPower,2)*size(matHip_GammaPower,1)]);
            cellHip_GammaPower = {arrayHip_GammaPower};
            matHip_MUAPower = cell2mat(Hip_muaPower);
            arrayHip_MUAPower = reshape(matHip_MUAPower',[1,size(matHip_MUAPower,2)*size(matHip_MUAPower,1)]);
            cellHip_MUAPower = {arrayHip_MUAPower};
            % vessel diameter
            matVesselDiameter = cell2mat(vesselDiameter);
            arrayVesselDiameter = reshape(matVesselDiameter',[1,size(matVesselDiameter,2)*size(matVesselDiameter,1)]);
            cellVesselDiameter = {arrayVesselDiameter};
            % whisker acceleration
            for cc = 1:length(WhiskerAcceleration)
                targetPoints = size(WhiskerAcceleration{1,1},2);
                if size(WhiskerAcceleration{cc,1},2) ~= targetPoints
                    maxLength = size(WhiskerAcceleration{cc,1},2);
                    difference = targetPoints - size(WhiskerAcceleration{cc,1},2);
                    for dd = 1:difference
                        WhiskerAcceleration{cc,1}(maxLength + dd) = 0;
                    end
                end
            end
            matWhiskerAcceleration = cell2mat(WhiskerAcceleration);
            arrayWhiskerAcceleration = reshape(matWhiskerAcceleration',[1,size(matWhiskerAcceleration,2)*size(matWhiskerAcceleration,1)]);
            cellWhiskerAcceleration = {arrayWhiskerAcceleration};
            % bin times
            matBinTimes = cell2mat(BinTimes);
            arrayBinTimes = reshape(matBinTimes',[1,size(matBinTimes,2)*size(matBinTimes,1)]);
            cellBinTimes = {arrayBinTimes};
        else
            count = length(fixedSleepIndex);
            holdIndex = zeros(1,(length(indexBreaks) + 1));
            for ee = 1:length(indexBreaks) + 1
                if ee == 1
                    holdIndex(ee) = indexBreaks(ee);
                elseif ee == length(indexBreaks) + 1
                    holdIndex(ee) = count - indexBreaks(ee - 1);
                else
                    holdIndex(ee)= indexBreaks(ee) - indexBreaks(ee - 1);
                end
            end
            splitCounter = 1:length(Cort_deltaPower);
            convertedMat2Cell = mat2cell(splitCounter',holdIndex);
            for ff = 1:length(convertedMat2Cell)
                % cortical neural bands
                mat2CellCort_DeltaPower{ff,1} = Cort_deltaPower(convertedMat2Cell{ff,1});
                mat2CellCort_ThetaPower{ff,1} = Cort_thetaPower(convertedMat2Cell{ff,1});
                mat2CellCort_AlphaPower{ff,1} = Cort_alphaPower(convertedMat2Cell{ff,1});
                mat2CellCort_BetaPower{ff,1} = Cort_betaPower(convertedMat2Cell{ff,1});
                mat2CellCort_GammaPower{ff,1} = Cort_gammaPower(convertedMat2Cell{ff,1});
                mat2CellCort_MUAPower{ff,1} = Cort_muaPower(convertedMat2Cell{ff,1});
                % hippocampal neural bands
                mat2CellHip_DeltaPower{ff,1} = Hip_deltaPower(convertedMat2Cell{ff,1});
                mat2CellHip_ThetaPower{ff,1} = Hip_thetaPower(convertedMat2Cell{ff,1});
                mat2CellHip_AlphaPower{ff,1} = Hip_alphaPower(convertedMat2Cell{ff,1});
                mat2CellHip_BetaPower{ff,1} = Hip_betaPower(convertedMat2Cell{ff,1});
                mat2CellHip_GammaPower{ff,1} = Hip_gammaPower(convertedMat2Cell{ff,1});
                mat2CellHip_MUAPower{ff,1} = Hip_muaPower(convertedMat2Cell{ff,1});
                % vessel diameter
                mat2CellVesselDiameter{ff,1} = vesselDiameter(convertedMat2Cell{ff,1});
                % whisker angle
                mat2CellWhiskerAcceleration{ff,1} = WhiskerAcceleration(convertedMat2Cell{ff,1});
                % bin times
                mat2CellBinTimes{ff,1} = BinTimes(convertedMat2Cell{ff,1});
            end
            for gg = 1:length(mat2CellCort_DeltaPower)
                % cortical neural bands
                matCort_DeltaPower = cell2mat(mat2CellCort_DeltaPower{gg,1});
                arrayCort_DeltaPower = reshape(matCort_DeltaPower',[1,size(matCort_DeltaPower,2)*size(matCort_DeltaPower,1)]);
                cellCort_DeltaPower{gg,1} = arrayCort_DeltaPower;
                matCort_ThetaPower = cell2mat(mat2CellCort_ThetaPower{gg,1});
                arrayCort_ThetaPower = reshape(matCort_ThetaPower',[1,size(matCort_ThetaPower,2)*size(matCort_ThetaPower,1)]);
                cellCort_ThetaPower{gg,1} = arrayCort_ThetaPower;
                matCort_AlphaPower = cell2mat(mat2CellCort_AlphaPower{gg,1});
                arrayCort_AlphaPower = reshape(matCort_AlphaPower',[1,size(matCort_AlphaPower,2)*size(matCort_AlphaPower,1)]);
                cellCort_AlphaPower{gg,1} = arrayCort_AlphaPower;
                matCort_BetaPower = cell2mat(mat2CellCort_BetaPower{gg,1});
                arrayCort_BetaPower = reshape(matCort_BetaPower',[1,size(matCort_BetaPower,2)*size(matCort_BetaPower,1)]);
                cellCort_BetaPower{gg,1} = arrayCort_BetaPower;
                matCort_GammaPower = cell2mat(mat2CellCort_GammaPower{gg,1});
                arrayCort_GammaPower = reshape(matCort_GammaPower',[1,size(matCort_GammaPower,2)*size(matCort_GammaPower,1)]);
                cellCort_GammaPower{gg,1} = arrayCort_GammaPower;
                matCort_MUAPower = cell2mat(mat2CellCort_MUAPower{gg,1});
                arrayCort_MUAPower = reshape(matCort_MUAPower',[1,size(matCort_MUAPower,2)*size(matCort_MUAPower,1)]);
                cellCort_MUAPower{gg,1} = arrayCort_MUAPower;
                % hippocampal neural bands
                matHip_DeltaPower = cell2mat(mat2CellHip_DeltaPower{gg,1});
                arrayHip_DeltaPower = reshape(matHip_DeltaPower',[1,size(matHip_DeltaPower,2)*size(matHip_DeltaPower,1)]);
                cellHip_DeltaPower{gg,1} = arrayHip_DeltaPower;
                matHip_ThetaPower = cell2mat(mat2CellHip_ThetaPower{gg,1});
                arrayHip_ThetaPower = reshape(matHip_ThetaPower',[1,size(matHip_ThetaPower,2)*size(matHip_ThetaPower,1)]);
                cellHip_ThetaPower{gg,1} = arrayHip_ThetaPower;
                matHip_AlphaPower = cell2mat(mat2CellHip_AlphaPower{gg,1});
                arrayHip_AlphaPower = reshape(matHip_AlphaPower',[1,size(matHip_AlphaPower,2)*size(matHip_AlphaPower,1)]);
                cellHip_AlphaPower{gg,1} = arrayHip_AlphaPower;
                matHip_BetaPower = cell2mat(mat2CellHip_BetaPower{gg,1});
                arrayHip_BetaPower = reshape(matHip_BetaPower',[1,size(matHip_BetaPower,2)*size(matHip_BetaPower,1)]);
                cellHip_BetaPower{gg,1} = arrayHip_BetaPower;
                matHip_GammaPower = cell2mat(mat2CellHip_GammaPower{gg,1});
                arrayHip_GammaPower = reshape(matHip_GammaPower',[1,size(matHip_GammaPower,2)*size(matHip_GammaPower,1)]);
                cellHip_GammaPower{gg,1} = arrayHip_GammaPower;
                matHip_MUAPower = cell2mat(mat2CellHip_MUAPower{gg,1});
                arrayHip_MUAPower = reshape(matHip_MUAPower',[1,size(matHip_MUAPower,2)*size(matHip_MUAPower,1)]);
                cellHip_MUAPower{gg,1} = arrayHip_MUAPower;
                % vessel diameter
                matVesselDiameter = cell2mat(mat2CellVesselDiameter{gg,1});
                arrayVesselDiameter = reshape(matVesselDiameter',[1,size(matVesselDiameter,2)*size(matVesselDiameter,1)]);
                cellVesselDiameter{gg,1} = arrayVesselDiameter;
                % whisker acceleration
                for hh = 1:size(mat2CellWhiskerAcceleration{gg,1},1)
                    targetPoints = size(mat2CellWhiskerAcceleration{gg,1}{1,1},2);
                    if size(mat2CellWhiskerAcceleration{gg,1}{hh,1},2) ~= targetPoints
                        maxLength = size(mat2CellWhiskerAcceleration{gg,1}{hh,1},2);
                        difference = targetPoints - size(mat2CellWhiskerAcceleration{gg,1}{hh,1},2);
                        for dd = 1:difference
                            mat2CellWhiskerAcceleration{gg,1}{hh,1}(maxLength + dd) = 0;
                        end
                    end
                end
                matWhiskerAcceleration = cell2mat(mat2CellWhiskerAcceleration{gg,1});
                arrayWhiskerAcceleration = reshape(matWhiskerAcceleration',[1,size(matWhiskerAcceleration,2)*size(matWhiskerAcceleration,1)]);
                cellWhiskerAcceleration{gg,1} = arrayWhiskerAcceleration;
                % bin times
                matBinTimes = cell2mat(mat2CellBinTimes{gg,1});
                arrayBinTimes = reshape(matBinTimes',[1,size(matBinTimes,2)*size(matBinTimes,1)]);
                cellBinTimes{gg,1} = arrayBinTimes;
            end
        end        
        % if the structure is empty, we need a special case to format the struct properly
        if isfield(SleepData,(modelName)) == false
            % loop through however many sleep epochs this file has
            for ii = 1:size(cellCort_DeltaPower,2)
                % cortical neural bands
                SleepData.(modelName).NREM.data.corticalNeural.deltaBandPower{ii,1} = cellCort_DeltaPower{1,1};
                SleepData.(modelName).NREM.data.corticalNeural.thetaBandPower{ii,1} = cellCort_ThetaPower{1,1};
                SleepData.(modelName).NREM.data.corticalNeural.alphaBandPower{ii,1} = cellCort_AlphaPower{1,1};
                SleepData.(modelName).NREM.data.corticalNeural.betaBandPower{ii,1} = cellCort_BetaPower{1,1};
                SleepData.(modelName).NREM.data.corticalNeural.gammaBandPower{ii,1} = cellCort_GammaPower{1,1};
                SleepData.(modelName).NREM.data.corticalNeural.muaPower{ii,1} = cellCort_MUAPower{1,1};
                % hippocampal neural bands
                SleepData.(modelName).NREM.data.hippocampalNeural.deltaBandPower{ii,1} = cellHip_DeltaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.thetaBandPower{ii,1} = cellHip_ThetaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.alphaBandPower{ii,1} = cellHip_AlphaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.betaBandPower{ii,1} = cellHip_BetaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.gammaBandPower{ii,1} = cellHip_GammaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.muaPower{ii,1} = cellHip_MUAPower{1,1};
                % vessel diameter
                SleepData.(modelName).NREM.data.vesselDiameter.data{ii,1} = cellVesselDiameter{1,1};
                % whisker acceleration
                SleepData.(modelName).NREM.data.WhiskerAcceleration{ii,1} = cellWhiskerAcceleration{1,1};
                % fileIDs and bin times
                SleepData.(modelName).NREM.FileIDs{ii,1} = fileID;
                SleepData.(modelName).NREM.VesselIDs{ii,1} = vesselID;
                SleepData.(modelName).NREM.BinTimes{ii,1} = cellBinTimes{1,1};
            end
            % if the struct is not empty, add each new iteration after previous data
        else
            % loop through however many sleep epochs this file has
            for jj = 1:size(cellCort_DeltaPower,1)
                % cortical neural bands
                SleepData.(modelName).NREM.data.corticalNeural.deltaBandPower{size(SleepData.(modelName).NREM.data.corticalNeural.deltaBandPower,1) + 1,1} = cellCort_DeltaPower{jj,1};
                SleepData.(modelName).NREM.data.corticalNeural.thetaBandPower{size(SleepData.(modelName).NREM.data.corticalNeural.thetaBandPower,1) + 1,1} = cellCort_ThetaPower{jj,1};
                SleepData.(modelName).NREM.data.corticalNeural.alphaBandPower{size(SleepData.(modelName).NREM.data.corticalNeural.alphaBandPower,1) + 1,1} = cellCort_AlphaPower{jj,1};
                SleepData.(modelName).NREM.data.corticalNeural.betaBandPower{size(SleepData.(modelName).NREM.data.corticalNeural.betaBandPower,1) + 1,1} = cellCort_BetaPower{jj,1};
                SleepData.(modelName).NREM.data.corticalNeural.gammaBandPower{size(SleepData.(modelName).NREM.data.corticalNeural.gammaBandPower,1) + 1,1} = cellCort_GammaPower{jj,1};
                SleepData.(modelName).NREM.data.corticalNeural.muaPower{size(SleepData.(modelName).NREM.data.corticalNeural.muaPower,1) + 1,1} = cellCort_MUAPower{jj,1};
                % hippocampal neural bands
                SleepData.(modelName).NREM.data.hippocampalNeural.deltaBandPower{size(SleepData.(modelName).NREM.data.hippocampalNeural.deltaBandPower,1) + 1,1} = cellHip_DeltaPower{jj,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.thetaBandPower{size(SleepData.(modelName).NREM.data.hippocampalNeural.thetaBandPower,1) + 1,1} = cellHip_ThetaPower{jj,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.alphaBandPower{size(SleepData.(modelName).NREM.data.hippocampalNeural.alphaBandPower,1) + 1,1} = cellHip_AlphaPower{jj,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.betaBandPower{size(SleepData.(modelName).NREM.data.hippocampalNeural.betaBandPower,1) + 1,1} = cellHip_BetaPower{jj,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.gammaBandPower{size(SleepData.(modelName).NREM.data.hippocampalNeural.gammaBandPower,1) + 1,1} = cellHip_GammaPower{jj,1};
                SleepData.(modelName).NREM.data.hippocampalNeural.muaPower{size(SleepData.(modelName).NREM.data.hippocampalNeural.muaPower,1) + 1,1} = cellHip_MUAPower{jj, 1};
                % vessel diameter
                SleepData.(modelName).NREM.data.vesselDiameter.data{size(SleepData.(modelName).NREM.data.vesselDiameter.data,1) + 1,1} = cellVesselDiameter{jj,1};
                % whisker acceleration
                SleepData.(modelName).NREM.data.WhiskerAcceleration{size(SleepData.(modelName).NREM.data.WhiskerAcceleration,1) + 1,1} = cellWhiskerAcceleration{jj,1};
                % file IDs and bin times
                SleepData.(modelName).NREM.FileIDs{size(SleepData.(modelName).NREM.FileIDs,1) + 1,1} = fileID;
                SleepData.(modelName).NREM.VesselIDs{size(SleepData.(modelName).NREM.VesselIDs,1) + 1,1} = vesselID;
                SleepData.(modelName).NREM.BinTimes{size(SleepData.(modelName).NREM.BinTimes,1) + 1,1} = cellBinTimes{jj,1};
            end
        end
    end
    disp(['Adding NREM sleeping epochs from MergedData file ' num2str(aa) ' of ' num2str(size(mergedDataFileIDs,1)) '...']); disp(' ')
end

%% create sleep scored data structure for REM data
% Identify sleep epochs and place in SleepData.mat structure
sleepBins = REMsleepTime/5;
for kk = 1:size(mergedDataFileIDs,1)
    clearvars -except kk mergedDataFileIDs sleepBins NREMsleepTime REMsleepTime modelName SleepData
    mergedDataFileID = mergedDataFileIDs(kk,:);
    load(mergedDataFileID);
    [~,~,~,fileID,~,vesselID] = GetFileInfo2_2P(mergedDataFileID);
    remLogical = MergedData.sleep.logicals.(modelName).remLogical;
    targetTime = ones(1,sleepBins);
    % find the periods of time where there are the desired number of consecutive bins
    sleepIndex = find(conv(remLogical,targetTime) >= sleepBins) - (sleepBins - 1);
    if ~isempty(sleepIndex)
        sleepCriteria = (0:(sleepBins - 1));
        % sleep Index now has the proper time stamps from sleep logical
        fixedSleepIndex = unique(sleepIndex + sleepCriteria);
        % loop through the length of sleep Index, and pull out associated data
        for ll = 1:length(fixedSleepIndex)
            % cortical neural bands
            Cort_deltaPower{ll,1} = MergedData.sleep.parameters.corticalNeural.deltaBandPower{fixedSleepIndex(ll),1}; %#ok<*AGROW>
            Cort_thetaPower{ll,1} = MergedData.sleep.parameters.corticalNeural.thetaBandPower{fixedSleepIndex(ll),1};
            Cort_alphaPower{ll,1} = MergedData.sleep.parameters.corticalNeural.alphaBandPower{fixedSleepIndex(ll),1};
            Cort_betaPower{ll,1} = MergedData.sleep.parameters.corticalNeural.betaBandPower{fixedSleepIndex(ll),1};
            Cort_gammaPower{ll,1} = MergedData.sleep.parameters.corticalNeural.gammaBandPower{fixedSleepIndex(ll),1};
            Cort_muaPower{ll,1} = MergedData.sleep.parameters.corticalNeural.muaPower{fixedSleepIndex(ll),1};
            % hippocampal neural bands
            Hip_deltaPower{ll,1} = MergedData.sleep.parameters.hippocampalNeural.deltaBandPower{fixedSleepIndex(ll),1};
            Hip_thetaPower{ll,1} = MergedData.sleep.parameters.hippocampalNeural.thetaBandPower{fixedSleepIndex(ll),1};
            Hip_alphaPower{ll,1} = MergedData.sleep.parameters.hippocampalNeural.alphaBandPower{fixedSleepIndex(ll),1};
            Hip_betaPower{ll,1} = MergedData.sleep.parameters.hippocampalNeural.betaBandPower{fixedSleepIndex(ll),1};
            Hip_gammaPower{ll,1} = MergedData.sleep.parameters.hippocampalNeural.gammaBandPower{fixedSleepIndex(ll),1};
            Hip_muaPower{ll,1} = MergedData.sleep.parameters.hippocampalNeural.muaPower{fixedSleepIndex(ll),1};
            % vessel diameter
            vesselDiameter{ll,1} = MergedData.sleep.parameters.vesselDiameter.data{fixedSleepIndex(ll),1};
            % whisker angle
            WhiskerAcceleration{ll,1} = MergedData.sleep.parameters.whiskerAcceleration{fixedSleepIndex(ll),1};
            % bin times
            BinTimes{ll,1} = 5*fixedSleepIndex(ll);
        end
        indexBreaks = find(fixedSleepIndex(2:end) - fixedSleepIndex(1:end - 1) > 1);    % Find if there are numerous sleep periods
        % if there is only one period of sleep in this file and not multiple events
        if isempty(indexBreaks)
            % cortical neural bands
            matCort_DeltaPower = cell2mat(Cort_deltaPower);
            arrayCort_DeltaPower = reshape(matCort_DeltaPower',[1,size(matCort_DeltaPower,2)*size(matCort_DeltaPower,1)]);
            cellCort_DeltaPower = {arrayCort_DeltaPower};
            matCort_ThetaPower = cell2mat(Cort_thetaPower);
            arrayCort_ThetaPower = reshape(matCort_ThetaPower',[1,size(matCort_ThetaPower,2)*size(matCort_ThetaPower,1)]);
            cellCort_ThetaPower = {arrayCort_ThetaPower};
            matCort_AlphaPower = cell2mat(Cort_alphaPower);
            arrayCort_AlphaPower = reshape(matCort_AlphaPower',[1,size(matCort_AlphaPower,2)*size(matCort_AlphaPower,1)]);
            cellCort_AlphaPower = {arrayCort_AlphaPower};
            matCort_BetaPower = cell2mat(Cort_betaPower);
            arrayCort_BetaPower = reshape(matCort_BetaPower',[1,size(matCort_BetaPower,2)*size(matCort_BetaPower,1)]);
            cellCort_BetaPower = {arrayCort_BetaPower};
            matCort_GammaPower = cell2mat(Cort_gammaPower);
            arrayCort_GammaPower = reshape(matCort_GammaPower',[1,size(matCort_GammaPower,2)*size(matCort_GammaPower,1)]);
            cellCort_GammaPower = {arrayCort_GammaPower};
            matCort_MUAPower = cell2mat(Cort_muaPower);
            arrayCort_MUAPower = reshape(matCort_MUAPower',[1,size(matCort_MUAPower,2)*size(matCort_MUAPower,1)]);
            cellCort_MUAPower = {arrayCort_MUAPower};
            % hippocampal neural bands
            matHip_DeltaPower = cell2mat(Hip_deltaPower);
            arrayHip_DeltaPower = reshape(matHip_DeltaPower',[1,size(matHip_DeltaPower,2)*size(matHip_DeltaPower,1)]);
            cellHip_DeltaPower = {arrayHip_DeltaPower};
            matHip_ThetaPower = cell2mat(Hip_thetaPower);
            arrayHip_ThetaPower = reshape(matHip_ThetaPower',[1,size(matHip_ThetaPower,2)*size(matHip_ThetaPower,1)]);
            cellHip_ThetaPower = {arrayHip_ThetaPower};
            matHip_AlphaPower = cell2mat(Hip_alphaPower);
            arrayHip_AlphaPower = reshape(matHip_AlphaPower',[1,size(matHip_AlphaPower,2)*size(matHip_AlphaPower,1)]);
            cellHip_AlphaPower = {arrayHip_AlphaPower};
            matHip_BetaPower = cell2mat(Hip_betaPower);
            arrayHip_BetaPower = reshape(matHip_BetaPower',[1,size(matHip_BetaPower,2)*size(matHip_BetaPower,1)]);
            cellHip_BetaPower = {arrayHip_BetaPower};
            matHip_GammaPower = cell2mat(Hip_gammaPower);
            arrayHip_GammaPower = reshape(matHip_GammaPower',[1,size(matHip_GammaPower,2)*size(matHip_GammaPower,1)]);
            cellHip_GammaPower = {arrayHip_GammaPower};
            matHip_MUAPower = cell2mat(Hip_muaPower);
            arrayHip_MUAPower = reshape(matHip_MUAPower',[1,size(matHip_MUAPower,2)*size(matHip_MUAPower,1)]);
            cellHip_MUAPower = {arrayHip_MUAPower};
            % vessel diameter
            matVesselDiameter = cell2mat(vesselDiameter);
            arrayVesselDiameter = reshape(matVesselDiameter',[1,size(matVesselDiameter,2)*size(matVesselDiameter,1)]);
            cellVesselDiameter = {arrayVesselDiameter};
            % whisker acceleration
            for mm = 1:length(WhiskerAcceleration)
                targetPoints = size(WhiskerAcceleration{1,1},2);
                if size(WhiskerAcceleration{mm,1},2) ~= targetPoints
                    maxLength = size(WhiskerAcceleration{mm,1},2);
                    difference = targetPoints - size(WhiskerAcceleration{mm,1},2);
                    for ss = 1:difference
                        WhiskerAcceleration{mm,1}(maxLength + ss) = 0;
                    end
                end
            end
            matWhiskerAcceleration = cell2mat(WhiskerAcceleration);
            arrayWhiskerAcceleration = reshape(matWhiskerAcceleration',[1,size(matWhiskerAcceleration,2)*size(matWhiskerAcceleration,1)]);
            cellWhiskerAcceleration = {arrayWhiskerAcceleration};
            % bin times
            matBinTimes = cell2mat(BinTimes);
            arrayBinTimes = reshape(matBinTimes',[1,size(matBinTimes,2)*size(matBinTimes,1)]);
            cellBinTimes = {arrayBinTimes};
        else
            count = length(fixedSleepIndex);
            holdIndex = zeros(1,(length(indexBreaks) + 1));
            for oo = 1:length(indexBreaks) + 1
                if oo == 1
                    holdIndex(oo) = indexBreaks(oo);
                elseif oo == length(indexBreaks) + 1
                    holdIndex(oo) = count - indexBreaks(oo - 1);
                else
                    holdIndex(oo)= indexBreaks(oo) - indexBreaks(oo - 1);
                end
            end
            splitCounter = 1:length(Cort_deltaPower);
            convertedMat2Cell = mat2cell(splitCounter',holdIndex);
            for pp = 1:length(convertedMat2Cell)
                % cortical neural bands
                mat2CellCort_DeltaPower{pp,1} = Cort_deltaPower(convertedMat2Cell{pp,1});
                mat2CellCort_ThetaPower{pp,1} = Cort_thetaPower(convertedMat2Cell{pp,1});
                mat2CellCort_AlphaPower{pp,1} = Cort_alphaPower(convertedMat2Cell{pp,1});
                mat2CellCort_BetaPower{pp,1} = Cort_betaPower(convertedMat2Cell{pp,1});
                mat2CellCort_GammaPower{pp,1} = Cort_gammaPower(convertedMat2Cell{pp,1});
                mat2CellCort_MUAPower{pp,1} = Cort_muaPower(convertedMat2Cell{pp,1});
                % hippocampal neural bands
                mat2CellHip_DeltaPower{pp,1} = Hip_deltaPower(convertedMat2Cell{pp,1});
                mat2CellHip_ThetaPower{pp,1} = Hip_thetaPower(convertedMat2Cell{pp,1});
                mat2CellHip_AlphaPower{pp,1} = Hip_alphaPower(convertedMat2Cell{pp,1});
                mat2CellHip_BetaPower{pp,1} = Hip_betaPower(convertedMat2Cell{pp,1});
                mat2CellHip_GammaPower{pp,1} = Hip_gammaPower(convertedMat2Cell{pp,1});
                mat2CellHip_MUAPower{pp,1} = Hip_muaPower(convertedMat2Cell{pp,1});
                % vessel diameter
                mat2CellVesselDiameter{pp,1} = vesselDiameter(convertedMat2Cell{pp,1});
                % whisker angle
                mat2CellWhiskerAcceleration{pp,1} = WhiskerAcceleration(convertedMat2Cell{pp,1});
                % bin times
                mat2CellBinTimes{pp,1} = BinTimes(convertedMat2Cell{pp,1});
            end
            for qq = 1:length(mat2CellCort_DeltaPower)
                % cortical neural bands
                matCort_DeltaPower = cell2mat(mat2CellCort_DeltaPower{qq,1});
                arrayCort_DeltaPower = reshape(matCort_DeltaPower',[1,size(matCort_DeltaPower,2)*size(matCort_DeltaPower,1)]);
                cellCort_DeltaPower{qq,1} = arrayCort_DeltaPower;
                matCort_ThetaPower = cell2mat(mat2CellCort_ThetaPower{qq,1});
                arrayCort_ThetaPower = reshape(matCort_ThetaPower',[1,size(matCort_ThetaPower,2)*size(matCort_ThetaPower,1)]);
                cellCort_ThetaPower{qq,1} = arrayCort_ThetaPower;
                matCort_AlphaPower = cell2mat(mat2CellCort_AlphaPower{qq,1});
                arrayCort_AlphaPower = reshape(matCort_AlphaPower',[1,size(matCort_AlphaPower,2)*size(matCort_AlphaPower,1)]);
                cellCort_AlphaPower{qq,1} = arrayCort_AlphaPower;
                matCort_BetaPower = cell2mat(mat2CellCort_BetaPower{qq,1});
                arrayCort_BetaPower = reshape(matCort_BetaPower',[1,size(matCort_BetaPower,2)*size(matCort_BetaPower,1)]);
                cellCort_BetaPower{qq,1} = arrayCort_BetaPower;
                matCort_GammaPower = cell2mat(mat2CellCort_GammaPower{qq,1});
                arrayCort_GammaPower = reshape(matCort_GammaPower',[1,size(matCort_GammaPower,2)*size(matCort_GammaPower,1)]);
                cellCort_GammaPower{qq,1} = arrayCort_GammaPower;
                matCort_MUAPower = cell2mat(mat2CellCort_MUAPower{qq,1});
                arrayCort_MUAPower = reshape(matCort_MUAPower',[1,size(matCort_MUAPower,2)*size(matCort_MUAPower,1)]);
                cellCort_MUAPower{qq,1} = arrayCort_MUAPower;
                % hippocampal neural bands
                matHip_DeltaPower = cell2mat(mat2CellHip_DeltaPower{qq,1});
                arrayHip_DeltaPower = reshape(matHip_DeltaPower',[1,size(matHip_DeltaPower,2)*size(matHip_DeltaPower,1)]);
                cellHip_DeltaPower{qq,1} = arrayHip_DeltaPower;
                matHip_ThetaPower = cell2mat(mat2CellHip_ThetaPower{qq,1});
                arrayHip_ThetaPower = reshape(matHip_ThetaPower',[1,size(matHip_ThetaPower,2)*size(matHip_ThetaPower,1)]);
                cellHip_ThetaPower{qq,1} = arrayHip_ThetaPower;
                matHip_AlphaPower = cell2mat(mat2CellHip_AlphaPower{qq,1});
                arrayHip_AlphaPower = reshape(matHip_AlphaPower',[1,size(matHip_AlphaPower,2)*size(matHip_AlphaPower,1)]);
                cellHip_AlphaPower{qq,1} = arrayHip_AlphaPower;
                matHip_BetaPower = cell2mat(mat2CellHip_BetaPower{qq,1});
                arrayHip_BetaPower = reshape(matHip_BetaPower',[1,size(matHip_BetaPower,2)*size(matHip_BetaPower,1)]);
                cellHip_BetaPower{qq,1} = arrayHip_BetaPower;
                matHip_GammaPower = cell2mat(mat2CellHip_GammaPower{qq,1});
                arrayHip_GammaPower = reshape(matHip_GammaPower',[1,size(matHip_GammaPower,2)*size(matHip_GammaPower,1)]);
                cellHip_GammaPower{qq,1} = arrayHip_GammaPower;
                matHip_MUAPower = cell2mat(mat2CellHip_MUAPower{qq,1});
                arrayHip_MUAPower = reshape(matHip_MUAPower',[1,size(matHip_MUAPower,2)*size(matHip_MUAPower,1)]);
                cellHip_MUAPower{qq,1} = arrayHip_MUAPower;
                % vessel diameter
                matVesselDiameter = cell2mat(mat2CellVesselDiameter{qq,1});
                arrayVesselDiameter = reshape(matVesselDiameter',[1,size(matVesselDiameter,2)*size(matVesselDiameter,1)]);
                cellVesselDiameter{qq,1} = arrayVesselDiameter;
                % whisker acceleration
                for rr = 1:size(mat2CellWhiskerAcceleration{qq,1},1)
                    targetPoints = size(mat2CellWhiskerAcceleration{qq,1}{1,1},2);
                    if size(mat2CellWhiskerAcceleration{qq,1}{rr,1},2) ~= targetPoints
                        maxLength = size(mat2CellWhiskerAcceleration{qq,1}{rr,1},2);
                        difference = targetPoints - size(mat2CellWhiskerAcceleration{qq,1}{rr,1},2);
                        for ss = 1:difference
                            mat2CellWhiskerAcceleration{qq,1}{rr,1}(maxLength + ss) = 0;
                        end
                    end
                end
                matWhiskerAcceleration = cell2mat(mat2CellWhiskerAcceleration{qq,1});
                arrayWhiskerAcceleration = reshape(matWhiskerAcceleration',[1,size(matWhiskerAcceleration,2)*size(matWhiskerAcceleration,1)]);
                cellWhiskerAcceleration{qq,1} = arrayWhiskerAcceleration;
                % bin times
                matBinTimes = cell2mat(mat2CellBinTimes{qq,1});
                arrayBinTimes = reshape(matBinTimes',[1,size(matBinTimes,2)*size(matBinTimes,1)]);
                cellBinTimes{qq,1} = arrayBinTimes;
            end
        end       
        % if the structure is empty, we need a special case to format the struct properly
        if isfield(SleepData.(modelName),'REM') == false
            % loop through however many sleep epochs this file has
            for tt = 1:size(cellCort_DeltaPower,2)
                % cortical neural bands
                SleepData.(modelName).REM.data.corticalNeural.deltaBandPower{tt,1} = cellCort_DeltaPower{1,1};
                SleepData.(modelName).REM.data.corticalNeural.thetaBandPower{tt,1} = cellCort_ThetaPower{1,1};
                SleepData.(modelName).REM.data.corticalNeural.alphaBandPower{tt,1} = cellCort_AlphaPower{1,1};
                SleepData.(modelName).REM.data.corticalNeural.betaBandPower{tt,1} = cellCort_BetaPower{1,1};
                SleepData.(modelName).REM.data.corticalNeural.gammaBandPower{tt,1} = cellCort_GammaPower{1,1};
                SleepData.(modelName).REM.data.corticalNeural.muaPower{tt,1} = cellCort_MUAPower{1,1};
                % hippocampal neural bands
                SleepData.(modelName).REM.data.hippocampalNeural.deltaBandPower{tt,1} = cellHip_DeltaPower{1,1};
                SleepData.(modelName).REM.data.hippocampalNeural.thetaBandPower{tt,1} = cellHip_ThetaPower{1,1};
                SleepData.(modelName).REM.data.hippocampalNeural.alphaBandPower{tt,1} = cellHip_AlphaPower{1,1};
                SleepData.(modelName).REM.data.hippocampalNeural.betaBandPower{tt,1} = cellHip_BetaPower{1,1};
                SleepData.(modelName).REM.data.hippocampalNeural.gammaBandPower{tt,1} = cellHip_GammaPower{1,1};
                SleepData.(modelName).REM.data.hippocampalNeural.muaPower{tt,1} = cellHip_MUAPower{1,1};
                % vessel diameter
                SleepData.(modelName).REM.data.vesselDiameter.data{tt,1} = cellVesselDiameter{1,1};
                % whisker acceleration
                SleepData.(modelName).REM.data.WhiskerAcceleration{tt,1} = cellWhiskerAcceleration{1,1};
                % fileIDs and bin times
                SleepData.(modelName).REM.FileIDs{tt,1} = fileID;
                SleepData.(modelName).REM.VesselIDs{tt,1} = vesselID;
                SleepData.(modelName).REM.BinTimes{tt,1} = cellBinTimes{1,1};
            end
            % if the struct is not empty, add each new iteration after previous data
        else
            % loop through however many sleep epochs this file has
            for uu = 1:size(cellCort_DeltaPower,1)
                % cortical neural bands
                SleepData.(modelName).REM.data.corticalNeural.deltaBandPower{size(SleepData.(modelName).REM.data.corticalNeural.deltaBandPower,1) + 1,1} = cellCort_DeltaPower{uu,1};
                SleepData.(modelName).REM.data.corticalNeural.thetaBandPower{size(SleepData.(modelName).REM.data.corticalNeural.thetaBandPower,1) + 1,1} = cellCort_ThetaPower{uu,1};
                SleepData.(modelName).REM.data.corticalNeural.alphaBandPower{size(SleepData.(modelName).REM.data.corticalNeural.alphaBandPower,1) + 1,1} = cellCort_AlphaPower{uu,1};
                SleepData.(modelName).REM.data.corticalNeural.betaBandPower{size(SleepData.(modelName).REM.data.corticalNeural.betaBandPower,1) + 1,1} = cellCort_BetaPower{uu,1};
                SleepData.(modelName).REM.data.corticalNeural.gammaBandPower{size(SleepData.(modelName).REM.data.corticalNeural.gammaBandPower,1) + 1,1} = cellCort_GammaPower{uu,1};
                SleepData.(modelName).REM.data.corticalNeural.muaPower{size(SleepData.(modelName).REM.data.corticalNeural.muaPower,1) + 1,1} = cellCort_MUAPower{uu,1};
                % hippocampal neural bands
                SleepData.(modelName).REM.data.hippocampalNeural.deltaBandPower{size(SleepData.(modelName).REM.data.hippocampalNeural.deltaBandPower,1) + 1,1} = cellHip_DeltaPower{uu,1};
                SleepData.(modelName).REM.data.hippocampalNeural.thetaBandPower{size(SleepData.(modelName).REM.data.hippocampalNeural.thetaBandPower,1) + 1,1} = cellHip_ThetaPower{uu,1};
                SleepData.(modelName).REM.data.hippocampalNeural.alphaBandPower{size(SleepData.(modelName).REM.data.hippocampalNeural.alphaBandPower,1) + 1,1} = cellHip_AlphaPower{uu,1};
                SleepData.(modelName).REM.data.hippocampalNeural.betaBandPower{size(SleepData.(modelName).REM.data.hippocampalNeural.betaBandPower,1) + 1,1} = cellHip_BetaPower{uu,1};
                SleepData.(modelName).REM.data.hippocampalNeural.gammaBandPower{size(SleepData.(modelName).REM.data.hippocampalNeural.gammaBandPower,1) + 1,1} = cellHip_GammaPower{uu,1};
                SleepData.(modelName).REM.data.hippocampalNeural.muaPower{size(SleepData.(modelName).REM.data.hippocampalNeural.muaPower,1) + 1,1} = cellHip_MUAPower{uu, 1};
                % vessel diameter
                SleepData.(modelName).REM.data.vesselDiameter.data{size(SleepData.(modelName).REM.data.vesselDiameter.data,1) + 1,1} = cellVesselDiameter{uu,1};
                % whisker acceleration
                SleepData.(modelName).REM.data.WhiskerAcceleration{size(SleepData.(modelName).REM.data.WhiskerAcceleration,1) + 1,1} = cellWhiskerAcceleration{uu,1};
                % file IDs and bin times
                SleepData.(modelName).REM.FileIDs{size(SleepData.(modelName).REM.FileIDs,1) + 1,1} = fileID;
                SleepData.(modelName).REM.VesselIDs{size(SleepData.(modelName).REM.VesselIDs,1) + 1,1} = vesselID;
                SleepData.(modelName).REM.BinTimes{size(SleepData.(modelName).REM.BinTimes,1) + 1,1} = cellBinTimes{uu,1};
            end
        end
    end
    disp(['Adding REM sleeping epochs from MergedData file ' num2str(kk) ' of ' num2str(size(mergedDataFileIDs,1)) '...']); disp(' ')
end
disp([modelName ' model data added to SleepData structure.']); disp(' ')

end
