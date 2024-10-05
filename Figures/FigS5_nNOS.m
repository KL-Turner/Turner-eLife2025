function [] = FigS5_nNOS(rootFolder,saveState,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% hemodynamic response function
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_HRF_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
dataTypes = {'gammaBandPower'};
behaviors = {'Contra','Whisk','Rest','NREM','REM','Alert','Asleep','All'};
variables = {'IR_function_long','IR_timeVec_long','IR_function_short','IR_timeVec_short','IR_gammaFunction','IR_gammaTimeVec','IR_R','IR_R2','group','animalID'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_HRF_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        % if any(strcmp(animalID,{'T180','T182'})) == false
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    ephysHRFData.(group).(hemisphere).(dataType).dummyCheck = 1;
                    if isfield(ephysHRFData.(group).(hemisphere).(dataType),(behavior)) == false
                        for ff = 1:length(variables)
                            variable = variables{1,ff};
                            if any(strcmp(variable,{'group','animalID'})) == true
                                ephysHRFData.(group).(hemisphere).(dataType).(behavior).(variable) = {};
                            else
                                ephysHRFData.(group).(hemisphere).(dataType).(behavior).(variable) = [];
                            end
                        end
                    end
                    if isempty(Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_R) == false
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_function_long = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_function_long,rescale(Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_function_long));
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_timeVec_long = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_timeVec_long,Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_timeVec_long);
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_function_short = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_function_short,rescale(Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_function_short));
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_timeVec_short = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_timeVec_short,Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_timeVec_short);
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_gammaFunction = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_gammaFunction,Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_gammaFunction);
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_gammaTimeVec = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_gammaTimeVec,Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_gammaTimeVec);
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_R = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_R,median(Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_R,'omitnan'));
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_R2 = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).IR_R2,median(Results_HRF_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).IR_R2,'omitnan'));
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).group = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).group,group);
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,ephysHRFData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
                    end
                end
                % end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    if any(strcmp(variable,{'group','animalID'})) == false
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(ephysHRFData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(ephysHRFData.(group).(hemisphere).(dataType).(behavior).(variable),0,1);
                    end
                    samplingRate = 30;
                    if strcmp(variable,{'IR_function_short'}) == true
                        [peak,peakIndex] = max(ephysHRFData.(group).(hemisphere).(dataType).(behavior).mean_IR_function_short);
                        peakTime = peakIndex/samplingRate;
                        threeQuarterMax = max(ephysHRFData.(group).(hemisphere).(dataType).(behavior).mean_IR_function_short)/(4/3);
                        index1 = find(ephysHRFData.(group).(hemisphere).(dataType).(behavior).mean_IR_function_short >= threeQuarterMax,1,'first');
                        % find where the data.Kernel last rises above half the max.
                        index2 = find(ephysHRFData.(group).(hemisphere).(dataType).(behavior).mean_IR_function_short >= threeQuarterMax,1,'last');
                        threeQuarterWidth = (index2 - index1 + 1)/samplingRate; % FWHM in indexes.
                        initVals = [peak,peakTime,threeQuarterWidth];
                        % create gamma function based on impulse values
                        t = 0:1/samplingRate:5;
                        IR_a = ((initVals(2)/initVals(3))^2*8*log10(2));
                        IR_beta = ((initVals(3)^2)/initVals(2)/8/log10(2));
                        ephysHRFData.(group).(hemisphere).(dataType).(behavior).repFunc = initVals(1)*(t/initVals(2)).^IR_a.*exp((t - initVals(2))/(-1*IR_beta));
                    end
                end
            end
        end
    end
end
% statistics - unpaired ttest
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            [HRFStats.(hemisphere).(dataType).(behavior).h,HRFStats.(hemisphere).(dataType).(behavior).p] = ttest2(ephysHRFData.Blank_SAP.(hemisphere).(dataType).(behavior).IR_R2,ephysHRFData.SSP_SAP.(hemisphere).(dataType).(behavior).IR_R2);
        end
    end
end

%% ephys power spectra
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PowerSpec_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'S','normS','f','binS','binf','group','animalID'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PowerSpec_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                ephysPowerData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(ephysPowerData.(group).(hemisphere).(dataType),behavior) == false
                        for ff = 1:length(variables)
                            variable = variables{1,ff};
                            if any(strcmp(variable,{'animalID','group','binf'})) == true
                                ephysPowerData.(group).(hemisphere).(dataType).(behavior).(variable) = {};
                            else
                                ephysPowerData.(group).(hemisphere).(dataType).(behavior).(variable) = [];
                            end
                        end
                    end
                    if isempty(Results_PowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).S) == false
                        ephysPowerData.(group).(hemisphere).(dataType).(behavior).S = cat(1,ephysPowerData.(group).(hemisphere).(dataType).(behavior).S,Results_PowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).S');
                        ephysPowerData.(group).(hemisphere).(dataType).(behavior).f = cat(1,ephysPowerData.(group).(hemisphere).(dataType).(behavior).f,Results_PowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).f);
                        freqBand = round(Results_PowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).f,2);
                        frequencyList = unique(freqBand);
                        for qq = 1:length(frequencyList)/2
                            freqIdx = find(freqBand == frequencyList(1,qq));
                            ephysPowerData.(group).(hemisphere).(dataType).(behavior).binS = cat(1,ephysPowerData.(group).(hemisphere).(dataType).(behavior).binS,mean(Results_PowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).S(freqIdx)));
                            ephysPowerData.(group).(hemisphere).(dataType).(behavior).binf = cat(1,ephysPowerData.(group).(hemisphere).(dataType).(behavior).binf,num2str(mean(freqBand(freqIdx))));
                            ephysPowerData.(group).(hemisphere).(dataType).(behavior).group = cat(1,ephysPowerData.(group).(hemisphere).(dataType).(behavior).group,group);
                            ephysPowerData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,ephysPowerData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
                        end
                    end
                end
            end
        end
    end
end
% find the peak of the resting PSD
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:size(ephysPowerData.(group).(hemisphere).(dataType).Rest.S,2)
                ephysPowerData.(group).(hemisphere).(dataType).baseline(dd,1) = max(ephysPowerData.(group).(hemisphere).(dataType).Rest.S(:,dd));
            end
        end
    end
end
% DC-shift each animal/hemisphere/behavior PSD with respect to the resting peak
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:size(ephysPowerData.(group).(hemisphere).(dataType).(behavior).S,1)
                    ephysPowerData.(group).(hemisphere).(dataType).(behavior).normS(ee,:) = (ephysPowerData.(group).(hemisphere).(dataType).(behavior).S(ee,:))*(1/(ephysPowerData.(group).(hemisphere).(dataType).baseline(ee,1)));
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    if any(strcmp(variable,{'normS','f'})) == true
                        ephysPowerData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(ephysPowerData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                        ephysPowerData.(group).(hemisphere).(dataType).(behavior).(['stdErr_' variable]) = std(ephysPowerData.(group).(hemisphere).(dataType).(behavior).(variable),0,1)./sqrt(size(ephysPowerData.(group).(hemisphere).(dataType).(behavior).(variable),1));
                    end
                end
            end
        end
    end
end
% GLME comparing peak correlation
for cc = 1:length(hemispheres)
    hemisphere = hemispheres{1,cc};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        for bb = 1:length(behaviors)
            behavior = behaviors{1,bb};
            EphysPowerStats.(hemisphere).(dataType).(behavior).tableSize = cat(1,ephysPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).binS,ephysPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).binS);
            EphysPowerStats.(hemisphere).(dataType).(behavior).Table = table('Size',[size(EphysPowerStats.(hemisphere).(dataType).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'AnimalID','Treatment','Frequency','Power'});
            EphysPowerStats.(hemisphere).(dataType).(behavior).Table.AnimalID = cat(1,ephysPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).animalID,ephysPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).animalID);
            EphysPowerStats.(hemisphere).(dataType).(behavior).Table.Treatment = cat(1,ephysPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).group,ephysPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).group);
            EphysPowerStats.(hemisphere).(dataType).(behavior).Table.Frequency = cat(1,ephysPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).binf,ephysPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).binf);
            EphysPowerStats.(hemisphere).(dataType).(behavior).Table.Power = cat(1,ephysPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).binS,ephysPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).binS);
            EphysPowerStats.(hemisphere).(dataType).(behavior).FitFormula = 'Power ~ 1 + Treatment + (1|Frequency) + (1|AnimalID)';
            EphysPowerStats.(hemisphere).(dataType).(behavior).Stats = fitglme(EphysPowerStats.(hemisphere).(dataType).(behavior).Table,EphysPowerStats.(hemisphere).(dataType).(behavior).FitFormula);
        end
    end
end

%% GCaMP power spectra
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PowerSpec_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'S','normS','f','binS','binf','group','animalID'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PowerSpec_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                gcampPowerData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(gcampPowerData.(group).(hemisphere).(dataType),behavior) == false
                        for ff = 1:length(variables)
                            variable = variables{1,ff};
                            if any(strcmp(variable,{'animalID','group','binf'})) == true
                                gcampPowerData.(group).(hemisphere).(dataType).(behavior).(variable) = {};
                            else
                                gcampPowerData.(group).(hemisphere).(dataType).(behavior).(variable) = [];
                            end
                        end
                    end
                    if isempty(Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).S) == false
                        gcampPowerData.(group).(hemisphere).(dataType).(behavior).S = cat(1,gcampPowerData.(group).(hemisphere).(dataType).(behavior).S,Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).S');
                        gcampPowerData.(group).(hemisphere).(dataType).(behavior).f = cat(1,gcampPowerData.(group).(hemisphere).(dataType).(behavior).f,Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).f);
                        freqBand = round(Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).f,2);
                        frequencyList = unique(freqBand);
                        for qq = 1:length(frequencyList)/2
                            freqIdx = find(freqBand == frequencyList(1,qq));
                            gcampPowerData.(group).(hemisphere).(dataType).(behavior).binS = cat(1,gcampPowerData.(group).(hemisphere).(dataType).(behavior).binS,mean(Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).S(freqIdx)));
                            gcampPowerData.(group).(hemisphere).(dataType).(behavior).binf = cat(1,gcampPowerData.(group).(hemisphere).(dataType).(behavior).binf,num2str(mean(freqBand(freqIdx))));
                            gcampPowerData.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampPowerData.(group).(hemisphere).(dataType).(behavior).group,group);
                            gcampPowerData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampPowerData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
                        end
                    end
                end
            end
        end
    end
end
% find the peak of the resting PSD
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:size(gcampPowerData.(group).(hemisphere).(dataType).Rest.S,2)
                gcampPowerData.(group).(hemisphere).(dataType).baseline(dd,1) = max(gcampPowerData.(group).(hemisphere).(dataType).Rest.S(:,dd));
            end
        end
    end
end
% DC-shift each animal/hemisphere/behavior PSD with respect to the resting peak
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:size(gcampPowerData.(group).(hemisphere).(dataType).(behavior).S,1)
                    gcampPowerData.(group).(hemisphere).(dataType).(behavior).normS(ee,:) = (gcampPowerData.(group).(hemisphere).(dataType).(behavior).S(ee,:))*(1/(gcampPowerData.(group).(hemisphere).(dataType).baseline(ee,1)));
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    if any(strcmp(variable,{'normS','f'})) == true
                        gcampPowerData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(gcampPowerData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                        gcampPowerData.(group).(hemisphere).(dataType).(behavior).(['stdErr_' variable]) = std(gcampPowerData.(group).(hemisphere).(dataType).(behavior).(variable),0,1)./sqrt(size(gcampPowerData.(group).(hemisphere).(dataType).(behavior).(variable),1));
                    end
                end
            end
        end
    end
end
% GLME comparing peak correlation
for cc = 1:length(hemispheres)
    hemisphere = hemispheres{1,cc};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        for bb = 1:length(behaviors)
            behavior = behaviors{1,bb};
            GCaMPPowerStats.(hemisphere).(dataType).(behavior).tableSize = cat(1,gcampPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).binS,gcampPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).binS);
            GCaMPPowerStats.(hemisphere).(dataType).(behavior).Table = table('Size',[size(GCaMPPowerStats.(hemisphere).(dataType).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'AnimalID','Treatment','Frequency','Power'});
            GCaMPPowerStats.(hemisphere).(dataType).(behavior).Table.AnimalID = cat(1,gcampPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).animalID,gcampPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).animalID);
            GCaMPPowerStats.(hemisphere).(dataType).(behavior).Table.Treatment = cat(1,gcampPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).group,gcampPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).group);
            GCaMPPowerStats.(hemisphere).(dataType).(behavior).Table.Frequency = cat(1,gcampPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).binf,gcampPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).binf);
            GCaMPPowerStats.(hemisphere).(dataType).(behavior).Table.Power = cat(1,gcampPowerData.Blank_SAP.(hemisphere).(dataType).(behavior).binS,gcampPowerData.SSP_SAP.(hemisphere).(dataType).(behavior).binS);
            GCaMPPowerStats.(hemisphere).(dataType).(behavior).FitFormula = 'Power ~ 1 + Treatment + (1|Frequency) + (1|AnimalID)';
            GCaMPPowerStats.(hemisphere).(dataType).(behavior).Stats = fitglme(GCaMPPowerStats.(hemisphere).(dataType).(behavior).Table,GCaMPPowerStats.(hemisphere).(dataType).(behavior).FitFormula);
        end
    end
end

%% Figure S5
FigS5 = figure('Name','Fig. S5');

% stim impulse HRF
subplot(2,3,1)
plot(ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.mean_IR_timeVec_long,ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.mean_IR_function_long,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.mean_IR_timeVec_long,ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.mean_IR_function_long,'color',colors('electric purple'),'LineWidth',2);
title('Stim impulse')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
set(gca,'box','off')
axis square

% stim gamma HRF
subplot(2,3,2)
plot(ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.mean_IR_timeVec_short,ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.repFunc,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.mean_IR_timeVec_short,ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.repFunc,'color',colors('electric purple'),'LineWidth',2);
title('Stim gamma-fit')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
set(gca,'box','off')
axis square

% stim R^2
subplot(2,3,3)
xInds = ones(1,length(ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.IR_R2));
scatter(xInds*1,ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.IR_R2,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.mean_IR_R2,ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.std_IR_R2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.IR_R2));
scatter(xInds*2,ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.IR_R2,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(2,ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.mean_IR_R2,ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.std_IR_R2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xlim([0,3])
title('Stim R^2')
ylabel('R^2')
set(gca,'box','off')
set(gca,'xtick',[])
axis square

% rest impulse HRF
subplot(2,3,4)
plot(ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.mean_IR_timeVec_long,ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.mean_IR_function_long,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.mean_IR_timeVec_long,ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.mean_IR_function_long,'color',colors('electric purple'),'LineWidth',2);
title('Rest impulse')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
set(gca,'box','off')
axis square

% rest gamma HRF
subplot(2,3,5)
plot(ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.mean_IR_timeVec_short,ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.repFunc,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.mean_IR_timeVec_short,ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.repFunc,'color',colors('electric purple'),'LineWidth',2);
title('Rest gamma-fit')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
set(gca,'box','off')
axis square

% rest R^2
subplot(2,3,6)
xInds = ones(1,length(ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.IR_R2));
scatter(xInds*1,ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.IR_R2,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.mean_IR_R2,ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.std_IR_R2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.IR_R2));
scatter(xInds*2,ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.IR_R2,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(2,ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.mean_IR_R2,ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.std_IR_R2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xlim([0,3])
title('Rest R^2')
ylabel('R^2')
set(gca,'box','off')
set(gca,'xtick',[])
axis square

%% save figure and stats
if saveState == true
    dirpath = [rootFolder delim 'MATLAB Figs/Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS5,[dirpath 'FigS5']);
    % statistical diary
    diaryFile = [dirpath 'FigS5_Stats.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on

    % stimulus-evoked HRF R2
    disp('======================================================================================================================')
    disp('stimulus-evoked HRF R2, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP: ' num2str(ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.mean_IR_R2) ' +/- ' num2str(ephysHRFData.Blank_SAP.RH.gammaBandPower.Contra.std_IR_R2)]); disp(' ')
    disp(['SSP-SAP: ' num2str(ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.mean_IR_R2) ' +/- ' num2str(ephysHRFData.SSP_SAP.RH.gammaBandPower.Contra.std_IR_R2)]); disp(' ')
    disp(['Blank vs. SAP ttest p = ' num2str(HRFStats.RH.gammaBandPower.Contra.p)]); disp(' ')

    % resting-state HRF R2
    disp('======================================================================================================================')
    disp('resting-state HRF R2, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP: ' num2str(ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.mean_IR_R2) ' +/- ' num2str(ephysHRFData.Blank_SAP.RH.gammaBandPower.Rest.std_IR_R2)]); disp(' ')
    disp(['SSP-SAP: ' num2str(ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.mean_IR_R2) ' +/- ' num2str(ephysHRFData.SSP_SAP.RH.gammaBandPower.Rest.std_IR_R2)]); disp(' ')
    disp(['Blank vs. SAP GLME p = ' num2str(HRFStats.RH.gammaBandPower.Rest.p)]); disp(' ')

    diary off
end
