function [Results_HRF_Ephys] = AnalyzeHRF_Ephys(animalID,group,set,rootFolder,delim,Results_HRF_Ephys)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
HRFLims_long = [-3,5];
HRFLims_short = [0,5];
HRFParams.offset = 3;
HRFParams.durIR = 10;
HRFParams.durGam = 5;
% event indecies for calc/test
EventInds.fitStart = 1;
EventInds.testStart = 2;
EventInds.increment = 2;
modelType = 'Forest';
samplingRate = 30;
% find and load RestData struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% find and load ManualDecisions struct
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% find and load EventData struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
% find and load RestingBaselines strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData strut
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% find and load Forest_ScoringResults.mat struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% loop variables
hemispheres = {'LH','RH'};
dataTypes = {'gammaBandPower','muaPower'};
behaviors = {'Contra','Whisk','Rest','NREM','REM','Alert','Asleep','All'};
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            %% get the valid neural and hemodynamic data from the data structures
            if strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true
                % extract neural and hemodynamic data from event structure
                [NeuralDataStruct,NeuralFiltArray] = SelectConvolutionBehavioralEvents_IOS(EventData.(['cortical_' hemisphere]).(dataType),behavior,hemisphere);
                [HemoDataStruct,HemoFiltArray] = SelectConvolutionBehavioralEvents_IOS(EventData.HbT.(hemisphere),behavior,hemisphere);
                % remove events that don't meet criteria
                [NeuralData,~,~,~] = RemoveInvalidData_IOS(NeuralDataStruct.NormData(NeuralFiltArray,:),NeuralDataStruct.fileIDs(NeuralFiltArray,:),NeuralDataStruct.duration(NeuralFiltArray,:),NeuralDataStruct.eventTime(NeuralFiltArray,:),ManualDecisions);
                [HemoData,~,~,~] = RemoveInvalidData_IOS(HemoDataStruct.data(HemoFiltArray,:),HemoDataStruct.fileIDs(HemoFiltArray,:),HemoDataStruct.duration(HemoFiltArray,:),HemoDataStruct.eventTime(HemoFiltArray,:),ManualDecisions);
            elseif strcmp(behavior,'Rest') == true
                % extract neural and hemodynamic data from event structure
                [NeuralDataStruct,NeuralFiltArray] = SelectConvolutionBehavioralEvents_IOS(RestData.(['cortical_' hemisphere]).(dataType),behavior,hemisphere);
                [HemoDataStruct,HemoFiltArray] = SelectConvolutionBehavioralEvents_IOS(RestData.HbT.(hemisphere),behavior,hemisphere);
                % remove events that don't meet criteria
                [NeuralData,~,~,~] = RemoveInvalidData_IOS(NeuralDataStruct.NormData(NeuralFiltArray,:),NeuralDataStruct.fileIDs(NeuralFiltArray,:),NeuralDataStruct.duration(NeuralFiltArray,:),NeuralDataStruct.eventTime(NeuralFiltArray,:),ManualDecisions);
                [HemoData,~,~,~] = RemoveInvalidData_IOS(HemoDataStruct.data(HemoFiltArray,:),HemoDataStruct.fileIDs(HemoFiltArray,:),HemoDataStruct.duration(HemoFiltArray,:),HemoDataStruct.eventTime(HemoFiltArray,:),ManualDecisions);
            elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
                % extract neural and hemodynamic data from sleep structure
                [NeuralData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).(behavior).data.(['cortical_' hemisphere]).(dataType),SleepData.(modelType).(behavior).FileIDs,SleepData.(modelType).(behavior).BinTimes);
                [HemoData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).(behavior).data.HbT.(hemisphere(end - 1:end)),SleepData.(modelType).(behavior).FileIDs,SleepData.(modelType).(behavior).BinTimes);
            elseif strcmp(behavior,'All') == true || strcmp(behavior,'Alert') == true || strcmp(behavior,'Asleep') == true
                [NeuralData,HemoData] = GatherLongEventData_IOS(dataType,hemisphere,behavior,RestingBaselines,ScoringResults);
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
                for dd = 1:length(NeuralData)
                    template = zeros(size(NeuralData{dd}));
                    strt = 2*samplingRate;
                    stp = size(template,2);
                    template(:,strt:stp) = detrend(NeuralData{dd}(:,strt:stp) - mean(NeuralData{dd}(:,strt:stp)));
                    processedNeuralFitData{dd} = template;
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
                for ee = 1:length(HemoData)
                    template = zeros(size(HemoData{ee}));
                    strt = 2*samplingRate;
                    stp = size(template,2);
                    hemoOffset = mean(HemoData{ee})*ones(1,stp - strt + 1);
                    template(:,strt:stp) = detrend(HemoData{ee}(:,strt:stp) - hemoOffset);
                    processedNeuralFitData{ee} = template;
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
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_function_long = timeLimHRF_long;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_timeVec_long = timeLimVec_long;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_function_short = timeLimHRF_short;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_timeVec_short = timeLimVec_short;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).HRFParams = HRFParams;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_numEvents = numEvents;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_eventIndex = EventInds;
            %% gamma HRF based on impulse deconvolution
            [peak,peakIndex] = max(Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_function_short);
            peakTime = peakIndex/samplingRate;
            threeQuarterMax = max(Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_function_short)/(4/3);
            index1 = find(Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_function_short >= threeQuarterMax,1,'first');
            % find where the data last rises above half the max.
            index2 = find(Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_function_short >= threeQuarterMax,1,'last');
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
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_gammaFunction = IR_gamma;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_gammaTimeVec = IR_t;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_gammaHRFParams = initVals;
            %% calculate the gamma HRF using fminsearch
            options = optimset('MaxFunEvals',2e4,'MaxIter',2e4,'TolFun',1e-7,'TolX',1e-7,'Display','off');
            [GamParams,~,~] = fminsearch(@(x)GammaConvolve_HRF2020(x,ProcNeuralFitData,ProcHemoFitData,samplingRate,HRFParams.durGam),initVals,options);
            FM_t = 0:1/samplingRate:HRFParams.durGam;
            FM_a = ((GamParams(2)/GamParams(3))^2*8*log10(2));
            FM_beta = ((GamParams(3)^2)/GamParams(2)/8/log10(2));
            FM_gamma = GamParams(1)*(FM_t/GamParams(2)).^FM_a.*exp((FM_t - GamParams(2))/(-1*FM_beta));
            % save results
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).FM_gammaFunction = FM_gamma;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).FM_gammaTimeVec = FM_t;
            Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).FM_gammaHRFParams = GamParams;
        end
    end
end
HRF_hemDataTypes = {'LH','RH'};
HRF_neuralBands = {'gammaBandPower','muaPower'};
behaviors = {'Contra','Whisk','Rest','NREM','REM','All','Alert','Asleep'};
for a = 1:length(HRF_hemDataTypes)
    hemDataType = HRF_hemDataTypes{1,a};
    for b = 1:length(HRF_neuralBands)
        neuralBand = HRF_neuralBands{1,b};
        for c = 1:length(behaviors)
            behavior = behaviors{1,c};
            [Results_HRF_Ephys] = EvaluateCBVPredictionAccuracy_HRF2020(animalID,group,neuralBand,hemDataType,behavior,Results_HRF_Ephys);
        end
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_HRF_Ephys.mat','Results_HRF_Ephys')
cd([rootFolder delim 'Data'])