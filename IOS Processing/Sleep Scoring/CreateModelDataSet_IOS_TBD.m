function [] = CreateModelDataSet_IOS(procDataFileIDs)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
binTime = 5; % sec
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    disp(['Creating model data set for ProcData file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    modelDataSetID = [procDataFileID(1:end - 12) 'ModelData.mat'];
    load(procDataFileID)
    % extract relevant parameters from each epoch
    for bb = 1:length(ProcData.notes.trialDuration_sec/binTime)
        % average cortical delta
        LH_cortDelta = mean(cell2mat(ProcData.sleep.parameters.cortical_LH.delta{bb,1}));
        RH_cortDelta = mean(cell2mat(ProcData.sleep.parameters.cortical_RH.delta{bb,1}));
        if LH_cortDelta >= RH_cortDelta
            cortDeltaColumn(bb,1) = LH_cortDelta;
        else
            cortDeltaColumn(bb,1) = RH_cortDelta;
        end
    end
    % average cortical beta
    LH_cortBeta = mean(cell2mat(ProcData.sleep.parameters.cortical_LH.beta{bb,1}));
    RH_cortBeta = mean(cell2mat(ProcData.sleep.parameters.cortical_RH.beta{bb,1}));
    if LH_cortBeta >= RH_cortBeta
        cortBetaColumn(bb,1) = LH_cortBeta;
    else
        cortBetaColumn(bb,1) = RH_cortBeta;
    end
    % average cortical gamma
    LH_cortGamma = mean(cell2mat(ProcData.sleep.parameters.cortical_LH.gamma{bb,1}));
    RH_cortGamma = mean(cell2mat(ProcData.sleep.parameters.cortical_RH.gamma{bb,1}));
    if LH_cortGamma >= RH_cortGamma
        cortGammaColumn(bb,1) = LH_cortGamma;
    else
        cortGammaColumn(bb,1) = RH_cortGamma;
    end
    % average hippocampal theta
    hippThetaColumn(bb,1) = mean(cell2mat(ProcData.sleep.parameters.hippocampus.theta{bb,1}));
    % number of binarized whisking events
    whiskEventsColumn(bb,1) = sum(ProcData.sleep.parameters.whiskerAngle.binarization{bb,1});
    % median of the EMG power
    EMGColumn(bb,1) = median(ProcData.sleep.parameters.EMG.power{bb,1});
    % average heart rate
    heartRateColumn(bb,1) = round(mean(ProcData.sleep.parameters.heartRate.frequency{bb,1}),1);
    variableNames = {'maxCortDelta','maxCortBeta','maxCortGamma','maxHippTheta','numWhiskEvents','avgEMG','avgHeartRate'};
    paramsTable = table(cortDeltaColumn,cortBetaColumn,cortGammaColumn,hippThetaColumn,whiskEventsColumn,EMGColumn,heartRateColumn,'VariableNames',variableNames);
    save(modelDataSetID,'paramsTable')
end