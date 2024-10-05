function [Results_HRF_Ephys] = EvaluateCBVPredictionAccuracy_HRF2020(animalID,group,neuralBand,hemisphere,behavior,Results_HRF_Ephys)
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
% filter characteristics
[z,p,k] = butter(4,0.5/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
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
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
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
    [NeuralDataStruct,NeuralFiltArray] = SelectConvolutionBehavioralEvents_HRF2020(BehData.(['cortical_' (hemisphere)]).(neuralBand),behavior,hemisphere);
    [HemoDataStruct,HemoFiltArray] = SelectConvolutionBehavioralEvents_HRF2020(BehData.HbT.(hemisphere),behavior,hemisphere);
    % remove events that don't meet criteria
    [NeuralData,~,~,~] = RemoveInvalidData_HRF2020(NeuralDataStruct.NormData(NeuralFiltArray,:),NeuralDataStruct.fileIDs(NeuralFiltArray,:),NeuralDataStruct.duration(NeuralFiltArray,:),NeuralDataStruct.eventTime(NeuralFiltArray,:),ManualDecisions);
    [HemoData,~,~,~] = RemoveInvalidData_HRF2020(HemoDataStruct.data(HemoFiltArray,:),HemoDataStruct.fileIDs(HemoFiltArray,:),HemoDataStruct.duration(HemoFiltArray,:),HemoDataStruct.eventTime(HemoFiltArray,:),ManualDecisions);
elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
    % extract neural and hemodynamic data from sleep structure
    [NeuralData,~,~] = RemoveStimSleepData_HRF2020(animalID,SleepData.(modelType).(behavior).data.(['cortical_' (hemisphere)]).(neuralBand),SleepData.(modelType).(behavior).FileIDs,SleepData.(modelType).(behavior).BinTimes);
    [HemoData,~,~] = RemoveStimSleepData_HRF2020(animalID,SleepData.(modelType).(behavior).data.HbT.((hemisphere)),SleepData.(modelType).(behavior).FileIDs,SleepData.(modelType).(behavior).BinTimes);
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
%% calculate R-squared on individual data
% pre-allocate
IR_IndR = NaN*ones(1,size(ProcHemoTestData,1));
IR_IndR2 = NaN*ones(1,size(ProcHemoTestData,1));
FM_IndR = NaN*ones(1,size(ProcHemoTestData,1));
FM_IndR2 = NaN*ones(1,size(ProcHemoTestData,1));
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true || strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
    % fit data
    for ee = 1:length(ProcHemoFitData)
        strt = 2*samplingRate;
        stp = length(ProcHemoFitData{ee});
        if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
            [IR_Act,IR_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_gammaFunction,detrend(ProcNeuralFitData{ee}),detrend(ProcHemoFitData{ee}),0);
            [FM_Act,FM_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_gammaFunction,detrend(ProcNeuralFitData{ee}),detrend(ProcHemoFitData{ee}),0);
        elseif strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
            [IR_Act,IR_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).All.IR_gammaFunction,detrend(ProcNeuralFitData{ee}),detrend(ProcHemoFitData{ee}),0);
            [FM_Act,FM_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).All.FM_gammaFunction,detrend(ProcNeuralFitData{ee}),detrend(ProcHemoFitData{ee}),0);
        end
        % impulse and fminsearch gamma functions
        IR_mPred = IR_Pred(strt:stp) - mean(IR_Pred(strt:stp));
        FM_mPred = FM_Pred(strt:stp) - mean(FM_Pred(strt:stp));
        IR_mAct = IR_Act(strt:stp) - mean(IR_Act(strt:stp));
        FM_mAct = FM_Act(strt:stp) - mean(FM_Act(strt:stp));
        FitData.IR_mPred{ee} = IR_mPred;
        FitData.FM_mPred{ee} = FM_mPred;
        FitData.IR_mAct{ee} = IR_mAct;
        FitData.FM_mAct{ee} = FM_mAct;
    end
    % test data
    for ff = 1:length(ProcHemoTestData)
        strt = 2*samplingRate;
        stp = length(ProcHemoTestData{ff});
        if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
            [IR_Act,IR_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_gammaFunction,detrend(ProcNeuralTestData{ff}),detrend(ProcHemoTestData{ff}),0);
            [FM_Act,FM_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_gammaFunction,detrend(ProcNeuralTestData{ff}),detrend(ProcHemoTestData{ff}),0);
        elseif strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
            [IR_Act,IR_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).All.IR_gammaFunction,detrend(ProcNeuralTestData{ff}),detrend(ProcHemoTestData{ff}),0);
            [FM_Act,FM_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).All.FM_gammaFunction,detrend(ProcNeuralTestData{ff}),detrend(ProcHemoTestData{ff}),0);
        end
        % impulse and fminsearch gamma functions
        IR_mPred = IR_Pred(strt:stp) - mean(IR_Pred(strt:stp));
        FM_mPred = IR_Pred(strt:stp) - mean(FM_Pred(strt:stp));
        IR_mAct = IR_Act(strt:stp) - mean(IR_Act(strt:stp));
        FM_mAct = IR_Act(strt:stp) - mean(FM_Act(strt:stp));
        TestData.IR_mPred{ff} = IR_mPred;
        TestData.FM_mPred{ff} = FM_mPred;
        TestData.IR_mAct{ff} = IR_mAct;
        TestData.FM_mAct{ff} = FM_mAct;
    end
    % caculate R and R2 for impulse gamma function
    try
        for gg = 1:length(TestData.IR_mPred)
            try
                IR_IndR(gg) = CalculateR_HRF2020(filtfilt(sos,g,TestData.IR_mPred{1,gg}),filtfilt(sos,g,TestData.IR_mAct{1,gg}));
                IR_IndR2(gg) = CalculateRsquared_HRF2020(filtfilt(sos,g,TestData.IR_mPred{1,gg}),filtfilt(sos,g,TestData.IR_mAct{1,gg}));
            catch
                IR_IndR(gg) = CalculateR_HRF2020(TestData.IR_mPred{1,gg},TestData.IR_mAct{1,gg});
                IR_IndR2(gg) = CalculateRsquared_HRF2020(TestData.IR_mPred{1,gg},TestData.IR_mAct{1,gg});
            end
        end
    catch
        IR_IndR = NaN;
        IR_IndR2 = NaN;
    end
    % caculate R and R2 for fminsearch gamma function
    try
        for gg = 1:length(TestData.FM_mPred)
            try
                FM_IndR(gg) = CalculateR_HRF2020(filtfilt(sos,g,TestData.FM_mPred{1,gg}),filtfilt(sos,g,TestData.FM_mAct{1,gg}));
                FM_IndR2(gg) = CalculateRsquared_HRF2020(filtfilt(sos,g,TestData.FM_mPred{1,gg}),filtfilt(sos,g,TestData.FM_mAct{1,gg}));
            catch
                FM_IndR(gg) = CalculateR_HRF2020(TestData.FM_mPred{1,gg},TestData.FM_mAct{1,gg});
                FM_IndR2(gg) = CalculateRsquared_HRF2020(TestData.FM_mPred{1,gg},TestData.FM_mAct{1,gg});
            end
        end
    catch
        FM_IndR = NaN;
        FM_IndR2 = NaN;
    end
    % save results
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_R = IR_IndR;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_R = FM_IndR;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_R2 = IR_IndR2;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_R2 = FM_IndR2;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_meanR = NaN;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_meanR = NaN;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_meanR2 = NaN;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_meanR2 = NaN;
elseif strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true
    % fit data
    for hh = 1:size(ProcHemoFitData,1)
        [IR_Act,IR_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_gammaFunction,ProcNeuralFitData(hh,:),ProcHemoFitData(hh,:),0);
        [FM_Act,FM_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_gammaFunction,ProcNeuralFitData(hh,:),ProcHemoFitData(hh,:),0);
        IR_mPred = IR_Pred(strt:stp) - mean(IR_Pred(strt:stp));
        FM_mPred = FM_Pred(strt:stp) - mean(FM_Pred(strt:stp));
        IR_mAct = IR_Act(strt:stp) - mean(IR_Act(strt:stp));
        FM_mAct = FM_Act(strt:stp) - mean(FM_Act(strt:stp));
        FitData.IR_mPred{hh,1} = IR_mPred;
        FitData.FM_mPred{hh,1} = FM_mPred;
        FitData.IR_mAct{hh,1} = IR_mAct;
        FitData.FM_mAct{hh,1} = FM_mAct;
    end
    % test data
    for ii = 1:size(ProcHemoTestData,1)
        [IR_Act,IR_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_gammaFunction,ProcNeuralTestData(ii,:),ProcHemoTestData(ii,:),0);
        [FM_Act,FM_Pred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_gammaFunction,ProcNeuralTestData(ii,:),ProcHemoTestData(ii,:),0);
        IR_mPred = IR_Pred(strt:stp) - mean(IR_Pred(strt:stp));
        FM_mPred = FM_Pred(strt:stp) - mean(FM_Pred(strt:stp));
        IR_mAct = IR_Act(strt:stp) - mean(IR_Act(strt:stp));
        FM_mAct = FM_Act(strt:stp) - mean(FM_Act(strt:stp));
        TestData.IR_mPred{ii,1} = IR_mPred;
        TestData.FM_mPred{ii,1} = FM_mPred;
        TestData.IR_mAct{ii,1} = IR_mAct;
        TestData.FM_mAct{ii,1} = FM_mAct;
    end
    % calculate R and R2
    for jj = 1:length(TestData.IR_mPred)
        IR_IndR(jj) = CalculateR_HRF2020(filtfilt(sos,g,TestData.IR_mPred{jj,1}),filtfilt(sos,g,TestData.IR_mAct{jj,1}));
        FM_IndR(jj) = CalculateR_HRF2020(filtfilt(sos,g,TestData.FM_mPred{jj,1}),filtfilt(sos,g,TestData.FM_mAct{jj,1}));
        IR_IndR2(jj) = CalculateRsquared_HRF2020(filtfilt(sos,g,TestData.IR_mPred{jj,1}),filtfilt(sos,g,TestData.IR_mAct{jj,1}));
        FM_IndR2(jj) = CalculateRsquared_HRF2020(filtfilt(sos,g,TestData.FM_mPred{jj,1}),filtfilt(sos,g,TestData.FM_mAct{jj,1}));
    end
    % calculate R and R2 on average data
    [IR_meanAct,IR_meanPred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_gammaFunction,mean(ProcNeuralFitData(hh,:),1),mean(ProcHemoFitData(hh,:),1),0);
    [FM_meanAct,FM_meanPred] = ConvolveHRF_HRF2020(Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_gammaFunction,mean(ProcNeuralFitData(hh,:),1),mean(ProcHemoFitData(hh,:),1),0);
    IR_avgPred = IR_meanPred(strt:stp) - mean(IR_meanPred(strt:stp));
    FM_avgPred = FM_meanPred(strt:stp) - mean(FM_meanPred(strt:stp));
    IR_avgAct = IR_meanAct(strt:stp) - mean(IR_meanAct(strt:stp));
    FM_avgAct = FM_meanAct(strt:stp) - mean(FM_meanAct(strt:stp));
    IR_meanR = CalculateR_HRF2020(filtfilt(sos,g,IR_avgPred),filtfilt(sos,g,IR_avgAct));
    FM_meanR = CalculateR_HRF2020(filtfilt(sos,g,FM_avgPred),filtfilt(sos,g,FM_avgAct));
    IR_meanR2 = CalculateRsquared_HRF2020(filtfilt(sos,g,IR_avgPred),filtfilt(sos,g,IR_avgAct));
    FM_meanR2 = CalculateRsquared_HRF2020(filtfilt(sos,g,FM_avgPred),filtfilt(sos,g,FM_avgAct));
    % save results
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_R = IR_IndR;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_R = FM_IndR;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_R2 = IR_IndR2;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_R2 = FM_IndR2;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_meanR = IR_meanR;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_meanR = FM_meanR;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).IR_meanR2 = IR_meanR2;
    Results_HRF_Ephys.(group).(animalID).(hemisphere).(neuralBand).(behavior).FM_meanR2 = FM_meanR2;
end