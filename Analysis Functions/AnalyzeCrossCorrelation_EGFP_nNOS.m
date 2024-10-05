function [Results_CrossCorr_EGFP] = AnalyzeCrossCorrelation_EGFP_nNOS(animalID,group,set,rootFolder,delim,Results_CrossCorr_EGFP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% character list of ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% find and load RestData struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% find and load RestingBaselines struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
hemispheres = {'LH','RH'};
% go through each valid data type for arousal-based cross-correlation analysis
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    % lowpass filter
    samplingRate = 10;
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    lagTime = 5;
    maxLag = lagTime*samplingRate;
    %% All
    zz = 1;
    allHbT = []; allNeural = [];
    for cc = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(cc,:);
        load(procDataFileID,'-mat')
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.CBV_HbT.(hemisphere)(1:300*samplingRate),'constant'));
            allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP7s.(['cor' hemisphere])(1:300*samplingRate),'constant'));
            zz = zz + 1;
            allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.CBV_HbT.(hemisphere)(300*samplingRate + 1:600*samplingRate),'constant'));
            allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP7s.(['cor' hemisphere])(300*samplingRate + 1:600*samplingRate),'constant'));
            zz = zz + 1;
            allHbT{zz,1} = filtfilt(sos,g,detrend(ProcData.data.CBV_HbT.(hemisphere)(600*samplingRate + 1:end),'constant'));
            allNeural{zz,1} = filtfilt(sos,g,detrend(ProcData.data.GCaMP7s.(['cor' hemisphere])(600*samplingRate + 1:end),'constant'));
            zz = zz + 1;
        end
    end
    if isempty(allHbT) == false
        % run cross-correlation analysis - average through time
        for dd = 1:length(allHbT)
            allHbTarray = allHbT{dd,1};
            allNeuralarray = allNeural{dd,1};
            [allHbTvNeuralxcVals(dd,:),allNeural_lags] = xcorr(allHbTarray,allNeuralarray,maxLag,'coeff');
        end
        allMeanHbTvNeuralxcVals = mean(allHbTvNeuralxcVals,1);
        % save results
        Results_CrossCorr_EGFP.(group).(animalID).(hemisphere).All.lags = allNeural_lags;
        Results_CrossCorr_EGFP.(group).(animalID).(hemisphere).All.xcVals = allMeanHbTvNeuralxcVals;
    else
        % save results
        Results_CrossCorr_EGFP.(group).(animalID).(hemisphere).All.lags = [];
        Results_CrossCorr_EGFP.(group).(animalID).(hemisphere).All.xcVals = [];
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_CrossCorr_EGFP.mat','Results_CrossCorr_EGFP')
cd([rootFolder delim 'Data'])