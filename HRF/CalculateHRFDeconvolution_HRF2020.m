function [AnalysisResults] = CalculateHRFDeconvolution_HRF2020(animalID,neuralBand,hemisphere,behavior,saveFigs,rootFolder,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner - adapted from code written by Aaron T. Winder
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Calculate the hemodynamic response function from neural data using Toeplitz and fminsearch
%________________________________________________________________________________________________________________________

%% setup parameters and load initial data structures
if strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true || strcmp(behavior,'Rest') == true ||  strcmp(behavior,'All') == true
    HRFLims_long = [-3,5];
    HRFLims_short = [0,5];
    HRFParams.offset = 3;
    HRFParams.durIR = 10;
    HRFParams.durGam = 5;
elseif strcmp(behavior,'NREM') == true
    HRFLims_long = [-3,5];
    HRFLims_short = [0,5];
    HRFParams.offset = 3;
    HRFParams.durIR = 10;
    HRFParams.durGam = 5;
elseif strcmp(behavior,'REM') == true
    HRFLims_long = [-3,5];
    HRFLims_short = [0,5];
    HRFParams.offset = 3;
    HRFParams.durIR = 10;
    HRFParams.durGam = 5;
end
% event indecies for calc/test
EventInds.fitStart = 1;
EventInds.testStart = 2;
EventInds.increment = 2;
modelType = 'Forest';
samplingRate = 30;
% load the Resting baselines structure
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% load Manual baseline event information
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
% Insert padding of zeros with size equal to the HRF between individual
fitIndex = EventInds.fitStart:EventInds.increment:size(NeuralData,1);
zpadIndex = EventInds.testStart:EventInds.increment:size(NeuralData,1);
zpad = zeros(1,HRFParams.durIR*samplingRate);
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true || strcmp(behavior,'All') == true
    % set every other event to zeros for padding
    NeuralData(zpadIndex) = {zpad};
    % process each event individually since they are all different lengths
    processedNeuralFitData = cell(size(NeuralData));
    for aa = 1:length(NeuralData)
        template = zeros(size(NeuralData{aa}));
        strt = 2*samplingRate;
        stp = size(template,2);
        template(:,strt:stp) = detrend(NeuralData{aa}(:,strt:stp) - mean(NeuralData{aa}(:,strt:stp)));
        processedNeuralFitData{aa} = template;
    end
    ProcNeuralFitData = [processedNeuralFitData{:}];
elseif strcmp(behavior,'Whisk') == true
    neuralDataEnd = 5;   % seconds
    strt = round((NeuralDataStruct.epoch.offset - 1)*samplingRate);
    stp = strt + round(neuralDataEnd*samplingRate);
    template = zeros(size(NeuralData(fitIndex,:)));
    neuralOffset = mean(NeuralData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    template(:,strt:stp) = NeuralData(fitIndex,strt:stp) - neuralOffset;
    % padding between events
    neuralDataPad = [template,ones(length(fitIndex),1)*zpad];
    ProcNeuralFitData = reshape(neuralDataPad',1,numel(neuralDataPad));
elseif strcmp(behavior,'Contra') == true
    neuralDataEnd = 1.5;   % seconds
    strt = round((NeuralDataStruct.epoch.offset)*NeuralDataStruct.samplingRate);
    stp = strt + (neuralDataEnd*samplingRate);
    template = zeros(size(NeuralData(fitIndex,:)));
    neuralOffset = mean(NeuralData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    template(:,strt:stp) = NeuralData(fitIndex,strt:stp) - neuralOffset;
    % padding between events
    neuralDataPad = [template,ones(length(fitIndex),1)*zpad];
    ProcNeuralFitData = reshape(neuralDataPad',1,numel(neuralDataPad));
end
%% setup the hemodynamic data based on behavior
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true || strcmp(behavior,'All') == true
    % set every other event to zeros for padding
    HemoData(zpadIndex) = {zpad};
    % process each event individually since they are all different lengths
    processedNeuralFitData = cell(size(HemoData));
    for bb = 1:length(HemoData)
        template = zeros(size(HemoData{bb}));
        strt = 2*samplingRate;
        stp = size(template,2);
        hemoOffset = mean(HemoData{bb})*ones(1,stp - strt + 1);
        template(:,strt:stp) = detrend(HemoData{bb}(:,strt:stp) - hemoOffset);
        processedNeuralFitData{bb} = template;
    end
    ProcHemoFitData = [processedNeuralFitData{:}];
elseif strcmp(behavior,'Whisk')
    hemoDataEnd = 7;   % seconds
    strt = round((HemoDataStruct.epoch.offset - 1)*samplingRate);
    stp = strt + round(hemoDataEnd*samplingRate);
    template = zeros(size(HemoData(fitIndex,:)));
    hemoOffset = mean(HemoData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    template(:,strt:stp) = HemoData(fitIndex,strt:stp) - hemoOffset;
    % padding between events
    hemoDataPad = [template,ones(length(fitIndex),1)*zpad];
    ProcHemoFitData = reshape(hemoDataPad',1,numel(hemoDataPad));
elseif strcmp(behavior,'Contra')
    hemoDataEnd = 3;   % seconds
    strt = round(HemoDataStruct.epoch.offset*samplingRate);
    stp = strt + (hemoDataEnd*samplingRate);
    template = zeros(size(HemoData(fitIndex,:)));
    hemoOffset = mean(HemoData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    template(:,strt:stp) = HemoData(fitIndex,strt:stp) - hemoOffset;
    % padding between events
    hemoDataPad = [template,ones(length(fitIndex),1)*zpad];
    ProcHemoFitData = reshape(hemoDataPad',1,numel(hemoDataPad));
end
%% calculate HRF based on deconvolution
IR_est = ImpulseAnalytic_HRF2020(ProcNeuralFitData',ProcHemoFitData',HRFParams.offset*samplingRate,HRFParams.durIR*samplingRate);
HRF = sgolayfilt(IR_est.IR',3,samplingRate + 1);
timeVec = (1:length(HRF))/samplingRate - HRFParams.offset;
numEvents = length(fitIndex);
timeLims_long = timeVec >= HRFLims_long(1) & timeVec <= HRFLims_long(2);
timeLims_short = timeVec >= HRFLims_short(1) & timeVec <= HRFLims_short(2);
timeLimHRF_long = HRF(timeLims_long);
timeLimHRF_short = HRF(timeLims_short);
timeLimVec_long = timeVec(timeLims_long);
timeLimVec_short = timeVec(timeLims_short);
% save results
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_function_long = timeLimHRF_long;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_timeVec_long = timeLimVec_long;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_function_short = timeLimHRF_short;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_timeVec_short = timeLimVec_short;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).HRFParams = HRFParams;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_numEvents = numEvents;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_eventIndex = EventInds;
%% gamma HRF based on impulse deconvolution
[peak,peakIndex] = max(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_function_short);
peakTime = peakIndex/samplingRate;
threeQuarterMax = max(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_function_short)/(4/3);
index1 = find(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_function_short >= threeQuarterMax,1,'first');
% find where the data last rises above half the max.
index2 = find(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_function_short >= threeQuarterMax,1,'last');
threeQuarterWidth = (index2 - index1 + 1)/samplingRate; % FWHM in indexes.
initVals = [peak,peakTime,threeQuarterWidth];
if isempty(threeQuarterWidth) == true
    initVals = [1e-1,1,1];
end
% create gamma function based on impulse values
IR_t = 0:1/samplingRate:HRFParams.durGam;
IR_a = ((initVals(2)/initVals(3))^2*8*log10(2));
IR_beta = ((initVals(3)^2)/initVals(2)/8/log10(2));
IR_gamma = initVals(1)*(IR_t/initVals(2)).^IR_a.*exp((IR_t - initVals(2))/(-1*IR_beta));
% save results
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_gammaFunction = IR_gamma;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_gammaTimeVec = IR_t;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_gammaHRFParams = initVals;
%% calculate the gamma HRF using fminsearch
options = optimset('MaxFunEvals',2e4,'MaxIter',2e4,'TolFun',1e-7,'TolX',1e-7,'Display','off');
[GamParams,~,~] = fminsearch(@(x)GammaConvolve_HRF2020(x,ProcNeuralFitData,ProcHemoFitData,samplingRate,HRFParams.durGam),initVals,options);
FM_t = 0:1/samplingRate:HRFParams.durGam;
FM_a = ((GamParams(2)/GamParams(3))^2*8*log10(2));
FM_beta = ((GamParams(3)^2)/GamParams(2)/8/log10(2));
FM_gamma = GamParams(1)*(FM_t/GamParams(2)).^FM_a.*exp((FM_t - GamParams(2))/(-1*FM_beta));
% save results
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).FM_gammaFunction = FM_gamma;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).FM_gammaTimeVec = FM_t;
AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).FM_gammaHRFParams = GamParams;
%% Save figures if desired
if strcmp(saveFigs,'y') == true
    kernelFig = figure;
    sgtitle([animalID ' ' hemisphere ' ' neuralBand ' during ' behavior])
    p1 = plot(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_timeVec_long,AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_function_long);
    hold on
    p2 = plot(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_gammaTimeVec,AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).IR_gammaFunction);
    p3 = plot(AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).FM_gammaTimeVec,AnalysisResults.HRFs.(animalID).(neuralBand).(hemisphere).(behavior).FM_gammaFunction);
    yline(threeQuarterMax);
    try
        xline(index1/samplingRate);
        xline(index2/samplingRate);
    catch
        % none
    end
    ylabel('Amplitude (a.u.)')
    xlabel('Time (s)')
    legend([p1,p2,p3],'Impulse','Impulse gamma','fminsearch gamma','Location','NorthWest')
    axis square
    axis tight
    set(gca,'box','off')
    % save figures
     dirpath = [rootFolder delim animalID delim 'Figures' delim 'HRF Kernels' delim neuralBand delim behavior delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(kernelFig,[dirpath animalID '_' hemisphere '_' neuralBand '_' behavior '_HRFs']);
    close(kernelFig)
end

end
