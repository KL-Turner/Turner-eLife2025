function [AnalysisResults] = EvaluateCBVPredictionPower_HRF2020(animalID,neuralBand,hemisphere,behavior,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner - adapted from code written by Aaron T. Winder
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: convolve HRF with neural activity and estimate accuracy of predictions via R2
%________________________________________________________________________________________________________________________

%% setup parameters and load initial data structures
EventIndex.fitStart = 1;
EventIndex.testStart = 2;
EventIndex.increment = 2;
modelType = 'Forest';
samplingRate = 30;
% load the Resting baselines structure
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% find and load Manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID)
% find and load Forest_ScoringResults.mat struct
forestScoringResultsFileID = 'Forest_ScoringResults.mat';
load(forestScoringResultsFileID,'-mat')
%% load the data structure relevent to the analysis
if strcmp(behavior,'Rest')
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    BehData = RestData;
elseif strcmp(behavior,'Whisk') || strcmp(behavior,'Contra')
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
    BehData = EventData;
elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM')
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
end
%% get the valid neural and hemodynamic data from the data structures
if strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true || strcmp(behavior,'Rest') == true
    % extract neural and hemodynamic data from event structure
    if strcmp(hemisphere(1:4),'cort') == true
        [NeuralDataStruct,NeuralFiltArray] = SelectConvolutionBehavioralEvents_HRF2020(BehData.(['cortical_' hemisphere(end - 1:end)]).(neuralBand),behavior,hemisphere);
    elseif strcmp(hemisphere(1:4),'hipp') == true
        [NeuralDataStruct,NeuralFiltArray] = SelectConvolutionBehavioralEvents_HRF2020(BehData.hippocampus.(neuralBand),behavior,hemisphere);
    end
    [HemoDataStruct,HemoFiltArray] = SelectConvolutionBehavioralEvents_HRF2020(BehData.CBV_HbT.(['adj' hemisphere(end - 1:end)]),behavior,hemisphere);
    % remove events that don't meet criteria
    [NeuralData,~,~,~] = RemoveInvalidData_HRF2020(NeuralDataStruct.NormData(NeuralFiltArray,:),NeuralDataStruct.fileIDs(NeuralFiltArray,:),NeuralDataStruct.duration(NeuralFiltArray,:),NeuralDataStruct.eventTime(NeuralFiltArray,:),ManualDecisions);
    [HemoData,~,~,~] = RemoveInvalidData_HRF2020(HemoDataStruct.data(HemoFiltArray,:),HemoDataStruct.fileIDs(HemoFiltArray,:),HemoDataStruct.duration(HemoFiltArray,:),HemoDataStruct.eventTime(HemoFiltArray,:),ManualDecisions);
elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
    % extract neural and hemodynamic data from sleep structure
    if strcmp(hemisphere(1:4),'cort') == true
        [NeuralData,~,~] = RemoveStimSleepData_HRF2020(animalID,SleepData.(modelType).(behavior).data.(['cortical_' hemisphere(end - 1:end)]).(neuralBand),SleepData.(modelType).(behavior).FileIDs,SleepData.(modelType).(behavior).BinTimes);
    elseif strcmp(hemisphere(1:4),'hipp') == true
        [NeuralData,~,~] = RemoveStimSleepData_HRF2020(animalID,SleepData.(modelType).(behavior).data.hippocampus.(neuralBand),SleepData.(modelType).(behavior).FileIDs,SleepData.(modelType).(behavior).BinTimes);
    end
    [HemoData,~,~] = RemoveStimSleepData_HRF2020(animalID,SleepData.(modelType).(behavior).data.CBV_HbT.(hemisphere(end - 1:end)),SleepData.(modelType).(behavior).FileIDs,SleepData.(modelType).(behavior).BinTimes);
elseif strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
    [NeuralData,HemoData] = GatherAllData_HRF2020(neuralBand,hemisphere,behavior,RestingBaselines,ScoringResults);
end
%% setup the neural data based on behavior
fitIndex = EventIndex.fitStart:EventIndex.increment:size(NeuralData,1);
testIndex = EventIndex.testStart:EventIndex.increment:size(NeuralData,1);
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true || strcmp(behavior,'All') == true || strcmp(behavior,'Asleep') == true || strcmp(behavior,'Alert') == true
    % adjust index for alert/asleep behaviors
    if strcmp(behavior,'Asleep') == true || strcmp(behavior,'Alert') == true
        fitIndex = EventIndex.fitStart:1:size(NeuralData,1);
        testIndex = EventIndex.testStart:1:size(NeuralData,1);
    end
    % fit data
    neuralFitData = NeuralData(fitIndex);
    processedNeuralFitData = cell(size(neuralFitData));
    for aa = 1:length(neuralFitData)
        template = zeros(size(neuralFitData{aa}));
        strt = 2*samplingRate;
        stp = size(template,2);
        template(:,strt:stp) = neuralFitData{aa}(:,strt:stp) - mean(neuralFitData{aa}(:,strt:stp));
        processedNeuralFitData{aa} = template;
    end
    ProcNeuralFitData = processedNeuralFitData;
    % test data
    neuralTestData = NeuralData(testIndex);
    processedNeuralTestData = cell(size(neuralTestData));
    for bb = 1:length(neuralTestData)
        template = zeros(size(neuralTestData{bb}));
        strt = 2*samplingRate;
        stp = size(template,2);
        template(:,strt:stp) = neuralTestData{bb}(:,strt:stp) - mean(neuralTestData{bb}(:,strt:stp));
        processedNeuralTestData{bb} = template;
    end
    ProcNeuralTestData = processedNeuralTestData;
elseif strcmp(behavior,'Whisk') == true
    neuralDataEnd = 5;   % seconds
    strt = round((NeuralDataStruct.epoch.offset - 1)*samplingRate);
    stp = strt + round(neuralDataEnd*samplingRate);
    % fit data
    ProcNeuralFitData = zeros(size(NeuralData(fitIndex,:)));
    fitNeuralOffset = mean(NeuralData(fitIndex,1:strt),2)*ones(1,stp - strt+1);
    ProcNeuralFitData(:,strt:stp) = NeuralData(fitIndex,strt:stp) - fitNeuralOffset;
    % test data
    ProcNeuralTestData = zeros(size(NeuralData(testIndex,:)));
    testNeuralOffset = mean(NeuralData(testIndex,1:strt),2)*ones(1,stp - strt+1);
    ProcNeuralTestData(:,strt:stp) = NeuralData(testIndex,strt:stp) - testNeuralOffset;
elseif strcmp(behavior,'Contra') == true
    neuralDataEnd = 1.5;   % seconds
    strt = round((NeuralDataStruct.epoch.offset)*samplingRate);
    stp = strt + (neuralDataEnd*samplingRate);
    % fit data
    ProcNeuralFitData = zeros(size(NeuralData(fitIndex,:)));
    fitNeuralOffset = mean(NeuralData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    ProcNeuralFitData(:,strt:stp) = NeuralData(fitIndex,strt:stp) - fitNeuralOffset;
    % test data
    ProcNeuralTestData = zeros(size(NeuralData(testIndex,:)));
    testNeuralOffset = mean(NeuralData(testIndex,1:strt),2)*ones(1,stp - strt + 1);
    ProcNeuralTestData(:,strt:stp) = NeuralData(testIndex,strt:stp) - testNeuralOffset;
end
%% setup the hemodynamic data based on behavior
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true || strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
    % fit data
    fitHemoData = HemoData(fitIndex);
    processedHemoFitData = cell(size(fitHemoData));
    for cc = 1:length(fitHemoData)
        template = zeros(size(fitHemoData{cc}));
        strt = 2*samplingRate;
        stp = size(template,2);
        offset = mean(fitHemoData{cc})*ones(1,stp - strt+1);
        template(:,strt:stp) = detrend(fitHemoData{cc}(:,strt:stp) - offset);
        processedHemoFitData{cc} = template;
    end
    ProcHemoFitData = processedHemoFitData;
    % test data
    testHemoData = HemoData(testIndex);
    processedHemoTestData = cell(size(testHemoData));
    for dd = 1:length(testHemoData)
        template = zeros(size(testHemoData{dd}));
        strt = 2*samplingRate;
        stp = size(template,2);
        offset = mean(testHemoData{dd})*ones(1,stp - strt+1);
        template(:,strt:stp) = detrend(testHemoData{dd}(:,strt:stp) - offset);
        processedHemoTestData{dd} = template;
    end
    ProcHemoTestData = processedHemoTestData;
elseif strcmp(behavior,'Whisk')
    hemoDataEnd = 7;   % seconds
    strt = round((HemoDataStruct.epoch.offset - 1)*samplingRate);
    stp = strt + round(hemoDataEnd*samplingRate);
    % fit data
    ProcHemoFitData = zeros(size(HemoData(fitIndex,:)));
    fitHemoOffset = mean(HemoData(fitIndex,1:strt),2)*ones(1,stp - strt+1);
    ProcHemoFitData(:,strt:stp) = HemoData(fitIndex,strt:stp) - fitHemoOffset;
    % test data
    ProcHemoTestData = zeros(size(HemoData(testIndex,:)));
    testHemoOffset = mean(HemoData(testIndex,1:strt),2)*ones(1,stp - strt+1);
    ProcHemoTestData(:,strt:stp) = HemoData(testIndex,strt:stp) - testHemoOffset;
elseif strcmp(behavior,'Contra')
    hemoDataEnd = 3;   % seconds
    strt = round(HemoDataStruct.epoch.offset*samplingRate);
    stp = strt + (hemoDataEnd*samplingRate);
    % fit data
    ProcHemoFitData = zeros(size(HemoData(fitIndex,:)));
    fitHemoOffset = mean(HemoData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    ProcHemoFitData(:,strt:stp) = HemoData(fitIndex,strt:stp) - fitHemoOffset;
    % test data
    ProcHemoTestData = zeros(size(HemoData(testIndex,:)));
    testHemoOffset = mean(HemoData(testIndex,1:strt),2)*ones(1,stp - strt + 1);
    ProcHemoTestData(:,strt:stp) = HemoData(testIndex,strt:stp) - testHemoOffset;
end
%% calculate spectral power
% parameters for coherencyc - information available in function
if strcmp(behavior,'Rest') == true || strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true
    params.tapers = [1,1];   % Tapers [n, 2n - 1]
elseif strcmp(behavior,'NREM') == true
    params.tapers = [3,5];   % Tapers [n, 2n - 1]
elseif strcmp(behavior,'REM') == true
    params.tapers = [5,9];   % Tapers [n, 2n - 1]
elseif strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
    params.tapers = [10,19];   % Tapers [n, 2n - 1]
end
params.pad = 1;
params.Fs = samplingRate;
params.fpass = [0,2];   % Pass band [0, nyquist]
params.trialave = 1;
params.err = [2,0.05];
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true || strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
    % test data
    for ff = 1:length(ProcHemoTestData)
        strt = 2*samplingRate;
        stp = length(ProcHemoTestData{ff});
        if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
            [IR_Act,IR_Pred] = ConvolveHRF_HRF2020(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_gammaFunction,detrend(ProcNeuralTestData{ff}),detrend(ProcHemoTestData{ff}),0);
            [FM_Act,FM_Pred] = ConvolveHRF_HRF2020(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).FM_gammaFunction,detrend(ProcNeuralTestData{ff}),detrend(ProcHemoTestData{ff}),0);
        elseif strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
            [IR_Act,IR_Pred] = ConvolveHRF_HRF2020(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).All.IR_gammaFunction,detrend(ProcNeuralTestData{ff}),detrend(ProcHemoTestData{ff}),0);
            [FM_Act,FM_Pred] = ConvolveHRF_HRF2020(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).All.FM_gammaFunction,detrend(ProcNeuralTestData{ff}),detrend(ProcHemoTestData{ff}),0);
        end
        IR_mPred = IR_Pred(strt:stp) - mean(IR_Pred(strt:stp));
        FM_mPred = FM_Pred(strt:stp) - mean(FM_Pred(strt:stp));
        IR_mAct = IR_Act(strt:stp) - mean(IR_Act(strt:stp));
        FM_mAct = FM_Act(strt:stp) - mean(FM_Act(strt:stp));
        TestData.IR_mPred{ff,1} = IR_mPred;
        TestData.FM_mPred{ff,1} = FM_mPred;
        TestData.IR_mAct{ff,1} = IR_mAct;
        TestData.FM_mAct{ff,1} = FM_mAct;
    end
    if isempty(ProcHemoTestData) == false
        % reshape data
        for jj = 1:length(TestData.IR_mPred)
            if strcmp(behavior,'Rest')
                minTime = 10;
            elseif strcmp(behavior,'NREM') == true
                minTime = 30;
            elseif strcmp(behavior,'REM') == true
                minTime = 60;
            elseif strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
                minTime = 900;
            end
            IR_truncPred(:,jj) = TestData.IR_mPred{jj,1}(1:(minTime - 2)*samplingRate + 1); %#ok<*AGROW,*NASGU>
            FM_truncPred(:,jj) = TestData.FM_mPred{jj,1}(1:(minTime - 2)*samplingRate + 1);
            IR_truncAct(:,jj) = TestData.IR_mAct{jj,1}(1:(minTime - 2)*samplingRate + 1);
            FM_truncAct(:,jj) = TestData.FM_mAct{jj,1}(1:(minTime - 2)*samplingRate + 1);
        end
        % calculate the power of the desired signals
        [IR_pred_S,IR_pred_f,IR_pred_sErr] = mtspectrumc_HRF2020(IR_truncPred,params);
        [FM_pred_S,FM_pred_f,FM_pred_sErr] = mtspectrumc_HRF2020(FM_truncPred,params);
        [IR_act_S,IR_act_f,IR_act_sErr] = mtspectrumc_HRF2020(IR_truncAct,params);
        [FM_act_S,FM_act_f,FM_act_sErr] = mtspectrumc_HRF2020(FM_truncAct,params);
        % save results
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_pred_S = IR_pred_S;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_pred_S = FM_pred_S;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_pred_f = IR_pred_f;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_pred_f = FM_pred_f;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_pred_sErr = IR_pred_sErr;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_pred_sErr = FM_pred_sErr;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_act_S = IR_act_S;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_act_S = FM_act_S;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_act_f = IR_act_f;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_act_f = FM_act_f;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_act_sErr = IR_act_sErr;
        AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_act_sErr = FM_act_sErr;
    end
elseif strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true
    % test data
    for ii = 1:size(ProcHemoTestData,1)
        [IR_Act,IR_Pred] = ConvolveHRF_HRF2020(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_gammaFunction,ProcNeuralTestData(ii,:),ProcHemoTestData(ii,:),0);
        [FM_Act,FM_Pred] = ConvolveHRF_HRF2020(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).FM_gammaFunction,ProcNeuralTestData(ii,:),ProcHemoTestData(ii,:),0);
        IR_mPred = IR_Pred(strt:stp) - mean(IR_Pred(strt:stp));
        FM_mPred = FM_Pred(strt:stp) - mean(FM_Pred(strt:stp));
        IR_mAct = IR_Act(strt:stp) - mean(IR_Act(strt:stp));
        FM_mAct = FM_Act(strt:stp) - mean(FM_Act(strt:stp));
        TestData.IR_mPred(:,ii) = IR_mPred;
        TestData.FM_mPred(:,ii) = FM_mPred;
        TestData.IR_mAct(:,ii) = IR_mAct;
        TestData.FM_mAct(:,ii) = FM_mAct;
    end
    % calculate the coherence between desired signals
    [IR_pred_S,IR_pred_f,IR_pred_sErr] = mtspectrumc_HRF2020(TestData.IR_mPred,params);
    [FM_pred_S,FM_pred_f,FM_pred_sErr] = mtspectrumc_HRF2020(TestData.FM_mPred,params);
    [IR_act_S,IR_act_f,IR_act_sErr] = mtspectrumc_HRF2020(TestData.IR_mAct,params);
    [FM_act_S,FM_act_f,FM_act_sErr] = mtspectrumc_HRF2020(TestData.FM_mAct,params);
    % save results
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_pred_S = IR_pred_S;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_pred_S = FM_pred_S;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_pred_f = IR_pred_f;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_pred_f = FM_pred_f;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_pred_sErr = IR_pred_sErr;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_pred_sErr = FM_pred_sErr;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_act_S = IR_act_S;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_act_S = FM_act_S;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_act_f = IR_act_f;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_act_f = FM_act_f;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).IR_act_sErr = IR_act_sErr;
    AnalysisResults.Power.(animalID).(neuralBand).(hemisphere).(behavior).FM_act_sErr = FM_act_sErr;
end

end
