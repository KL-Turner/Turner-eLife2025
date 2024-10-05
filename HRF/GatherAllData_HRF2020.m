function [NeuralData,HemoData] = GatherAllData_HRF2020(neuralBand,hemisphere,behavior,RestingBaselines,ScoringResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner - adapted from code written by Aaron T. Winder
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: calculate the hemodynamic response function from neural data
%________________________________________________________________________________________________________________________

% Character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
NeuralData = {};
HemoData = {};
if strcmp(behavior,'All') == true
    % load each file and put processed data into each structure
    for aa = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(aa,:);
        load(procDataFileID,'-mat')
        [~,fileDate,~] = GetFileInfo_HRF2020(procDataFileID);
        strDay = ConvertDate_HRF2020(fileDate);
        NeuralData{aa,1} = (ProcData.data.(['cortical_' (hemisphere)]).(neuralBand) - RestingBaselines.manualSelection.(['cortical_' (hemisphere)]).(neuralBand).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' (hemisphere)]).(neuralBand).(strDay).mean;
        HemoData{aa,1} = ProcData.data.HbT.(hemisphere);
    end
elseif strcmp(behavior,'Alert') == true
    dd = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,alertDataFileDate,alertDataFileID] = GetFileInfo_HRF2020(procDataFileID);
        strDay = ConvertDate_HRF2020(alertDataFileDate);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(alertDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of sleep
            load(procDataFileID)
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                NeuralData{dd,1} = (ProcData.data.(['cortical_' (hemisphere)]).(neuralBand) - RestingBaselines.manualSelection.(['cortical_' (hemisphere)]).(neuralBand).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' (hemisphere)]).(neuralBand).(strDay).mean;
                HemoData{dd,1} = ProcData.data.HbT.(hemisphere);
                dd = dd + 1;
            end
        end
    end
elseif strcmp(behavior,'Asleep') == true
    gg = 1;
    for ee = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(ee,:);
        [~,alertDataFileDate,alertDataFileID] = GetFileInfo_HRF2020(procDataFileID);
        strDay = ConvertDate_HRF2020(alertDataFileDate);
        scoringLabels = [];
        for ff = 1:length(ScoringResults.fileIDs)
            if strcmp(alertDataFileID,ScoringResults.fileIDs{ff,1}) == true
                scoringLabels = ScoringResults.labels{ff,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
            load(procDataFileID)
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                NeuralData{gg,1} = (ProcData.data.(['cortical_' (hemisphere)]).(neuralBand) - RestingBaselines.manualSelection.(['cortical_' (hemisphere)]).(neuralBand).(strDay).mean)./RestingBaselines.manualSelection.(['cortical_' (hemisphere)]).(neuralBand).(strDay).mean;
                HemoData{gg,1} = ProcData.data.HbT.(hemisphere);
                gg = gg + 1;
            end
        end
    end
end

end
