function [Results_PowerSpec_LFP] = AnalyzePowerSpectrum_LFP(animalID,group,set,rootFolder,delim,Results_PowerSpec_LFP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% character list of RawData file IDs
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
% character list of ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% find and load Forest_ScoringResults struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% analyze power spectra during periods of alert/asleep/all
behaviors = {'Alert','Asleep','All'};
hemispheres = {'LH','RH'};
xx = 1; yy = 1; zz = 1;
analogFs = 20000;
dsFs = 1000;
params.tapers = [5,9]; % Tapers [n, 2n - 1]
params.pad = 1;
params.Fs = dsFs;
params.fpass = [1,100]; % Pass band [0, nyquist]
params.trialave = 1;
params.err = [2,0.05];
data.LH.Alert = []; data.LH.Asleep = []; data.LH.All = [];
data.RH.Alert = []; data.RH.Asleep = []; data.RH.All = [];
for bb = 1:size(rawDataFileIDs,1)
    rawDataFileID = rawDataFileIDs(bb,:);
    procDataFileID = procDataFileIDs(bb,:);
    [~,~,allDataFileID] = GetFileInfo_IOS(procDataFileID);
    scoringLabels = [];
    for cc = 1:length(ScoringResults.fileIDs)
        if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
            scoringLabels = ScoringResults.labels{cc,1};
        end
    end
    load(procDataFileID,'-mat')
    puffs = ProcData.data.stimulations.LPadSol;
    % don't include trials with stimulation
    if isempty(puffs) == true
        load(rawDataFileID,'-mat')
        motionArtifact = ProcData.notes.motionArtifact;
        if motionArtifact == false
            data.LH.All{xx,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
            data.RH.All{xx,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
            xx = xx + 1;
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % 36 bins (180 total) or 3 minutes of sleep
            if motionArtifact == false
                data.LH.Alert{yy,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
                data.RH.Alert{yy,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
                yy = yy + 1;
            end
        elseif sum(strcmp(scoringLabels,'Not Sleep')) < 36 % 36 bins (180 total) or 3 minutes of awake
            if motionArtifact == false
                data.LH.Asleep{zz,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
                data.RH.Asleep{zz,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
                zz = zz + 1;
            end
        end
    end
end
% calculate LFP power spectrum
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        procData = [];
        if isempty(data.(hemisphere).(behavior)) == false
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            for cc = 1:length(data.(hemisphere).(behavior))
                procData(:,cc) = detrend(data.(hemisphere).(behavior){cc,1},'constant');
            end
            % calculate the power spectra of the desired signals
            [S,f,sErr] = mtspectrumc(procData,params);
            % save results
            Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S = S;
            Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).f = f;
            Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).sErr = sErr;
        else
            % save results
            Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S = [];
            Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).f = [];
            Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).sErr = [];
        end
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PowerSpec_LFP.mat','Results_PowerSpec_LFP')
cd([rootFolder delim 'Data'])