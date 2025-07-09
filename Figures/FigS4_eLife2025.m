function [] = FigS4_eLife2025(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

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

%% Figure S4
FigS4 = figure('Name','Fig. S4');

% HbT power ('All')
subplot(2,6,1)
loglog(ephysPowerData.Blank_SAP.RH.HbT.All.mean_f,ephysPowerData.Blank_SAP.RH.HbT.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(ephysPowerData.Blank_SAP.RH.HbT.All.mean_f,ephysPowerData.Blank_SAP.RH.HbT.All.mean_normS + ephysPowerData.Blank_SAP.RH.HbT.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(ephysPowerData.Blank_SAP.RH.HbT.All.mean_f,ephysPowerData.Blank_SAP.RH.HbT.All.mean_normS - ephysPowerData.Blank_SAP.RH.HbT.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(ephysPowerData.SSP_SAP.RH.HbT.All.mean_f,ephysPowerData.SSP_SAP.RH.HbT.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(ephysPowerData.SSP_SAP.RH.HbT.All.mean_f,ephysPowerData.SSP_SAP.RH.HbT.All.mean_normS + ephysPowerData.SSP_SAP.RH.HbT.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(ephysPowerData.SSP_SAP.RH.HbT.All.mean_f,ephysPowerData.SSP_SAP.RH.HbT.All.mean_normS - ephysPowerData.SSP_SAP.RH.HbT.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Ipsi HbT Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% Gamma power ('All')
subplot(2,6,2)
loglog(ephysPowerData.Blank_SAP.RH.gammaBandPower.All.mean_f,ephysPowerData.Blank_SAP.RH.gammaBandPower.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(ephysPowerData.Blank_SAP.RH.gammaBandPower.All.mean_f,ephysPowerData.Blank_SAP.RH.gammaBandPower.All.mean_normS + ephysPowerData.Blank_SAP.RH.gammaBandPower.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(ephysPowerData.Blank_SAP.RH.gammaBandPower.All.mean_f,ephysPowerData.Blank_SAP.RH.gammaBandPower.All.mean_normS - ephysPowerData.Blank_SAP.RH.gammaBandPower.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(ephysPowerData.SSP_SAP.RH.gammaBandPower.All.mean_f,ephysPowerData.SSP_SAP.RH.gammaBandPower.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(ephysPowerData.SSP_SAP.RH.gammaBandPower.All.mean_f,ephysPowerData.SSP_SAP.RH.gammaBandPower.All.mean_normS + ephysPowerData.SSP_SAP.RH.gammaBandPower.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(ephysPowerData.SSP_SAP.RH.gammaBandPower.All.mean_f,ephysPowerData.SSP_SAP.RH.gammaBandPower.All.mean_normS - ephysPowerData.SSP_SAP.RH.gammaBandPower.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Ipsi Gamma Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% HbT power ('All')
subplot(2,6,3)
loglog(gcampPowerData.Blank_SAP.RH.HbT.All.mean_f,gcampPowerData.Blank_SAP.RH.HbT.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(gcampPowerData.Blank_SAP.RH.HbT.All.mean_f,gcampPowerData.Blank_SAP.RH.HbT.All.mean_normS + gcampPowerData.Blank_SAP.RH.HbT.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.Blank_SAP.RH.HbT.All.mean_f,gcampPowerData.Blank_SAP.RH.HbT.All.mean_normS - gcampPowerData.Blank_SAP.RH.HbT.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.RH.HbT.All.mean_f,gcampPowerData.SSP_SAP.RH.HbT.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(gcampPowerData.SSP_SAP.RH.HbT.All.mean_f,gcampPowerData.SSP_SAP.RH.HbT.All.mean_normS + gcampPowerData.SSP_SAP.RH.HbT.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.RH.HbT.All.mean_f,gcampPowerData.SSP_SAP.RH.HbT.All.mean_normS - gcampPowerData.SSP_SAP.RH.HbT.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Ipsi HbT Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% HbO power ('All')
subplot(2,6,4)
loglog(gcampPowerData.Blank_SAP.RH.HbO.All.mean_f,gcampPowerData.Blank_SAP.RH.HbO.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(gcampPowerData.Blank_SAP.RH.HbO.All.mean_f,gcampPowerData.Blank_SAP.RH.HbO.All.mean_normS + gcampPowerData.Blank_SAP.RH.HbO.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.Blank_SAP.RH.HbO.All.mean_f,gcampPowerData.Blank_SAP.RH.HbO.All.mean_normS - gcampPowerData.Blank_SAP.RH.HbO.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.RH.HbO.All.mean_f,gcampPowerData.SSP_SAP.RH.HbO.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(gcampPowerData.SSP_SAP.RH.HbO.All.mean_f,gcampPowerData.SSP_SAP.RH.HbO.All.mean_normS + gcampPowerData.SSP_SAP.RH.HbO.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.RH.HbO.All.mean_f,gcampPowerData.SSP_SAP.RH.HbO.All.mean_normS - gcampPowerData.SSP_SAP.RH.HbO.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Ipsi HbO Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% HbR power ('All')
subplot(2,6,5)
loglog(gcampPowerData.Blank_SAP.RH.HbR.All.mean_f,gcampPowerData.Blank_SAP.RH.HbR.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(gcampPowerData.Blank_SAP.RH.HbR.All.mean_f,gcampPowerData.Blank_SAP.RH.HbR.All.mean_normS + gcampPowerData.Blank_SAP.RH.HbR.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.Blank_SAP.RH.HbR.All.mean_f,gcampPowerData.Blank_SAP.RH.HbR.All.mean_normS - gcampPowerData.Blank_SAP.RH.HbR.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.RH.HbR.All.mean_f,gcampPowerData.SSP_SAP.RH.HbR.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(gcampPowerData.SSP_SAP.RH.HbR.All.mean_f,gcampPowerData.SSP_SAP.RH.HbR.All.mean_normS + gcampPowerData.SSP_SAP.RH.HbR.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.RH.HbR.All.mean_f,gcampPowerData.SSP_SAP.RH.HbR.All.mean_normS - gcampPowerData.SSP_SAP.RH.HbR.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Ipsi HbR Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% GCaMP power ('All')
subplot(2,6,6)
loglog(gcampPowerData.Blank_SAP.RH.GCaMP.All.mean_f,gcampPowerData.Blank_SAP.RH.GCaMP.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(gcampPowerData.Blank_SAP.RH.GCaMP.All.mean_f,gcampPowerData.Blank_SAP.RH.GCaMP.All.mean_normS + gcampPowerData.Blank_SAP.RH.GCaMP.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.Blank_SAP.RH.GCaMP.All.mean_f,gcampPowerData.Blank_SAP.RH.GCaMP.All.mean_normS - gcampPowerData.Blank_SAP.RH.GCaMP.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.RH.GCaMP.All.mean_f,gcampPowerData.SSP_SAP.RH.GCaMP.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(gcampPowerData.SSP_SAP.RH.GCaMP.All.mean_f,gcampPowerData.SSP_SAP.RH.GCaMP.All.mean_normS + gcampPowerData.SSP_SAP.RH.GCaMP.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.RH.GCaMP.All.mean_f,gcampPowerData.SSP_SAP.RH.GCaMP.All.mean_normS - gcampPowerData.SSP_SAP.RH.GCaMP.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Ipsi GCaMP Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% HbT power ('All')
subplot(2,6,7)
loglog(ephysPowerData.Blank_SAP.LH.HbT.All.mean_f,ephysPowerData.Blank_SAP.LH.HbT.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(ephysPowerData.Blank_SAP.LH.HbT.All.mean_f,ephysPowerData.Blank_SAP.LH.HbT.All.mean_normS + ephysPowerData.Blank_SAP.LH.HbT.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(ephysPowerData.Blank_SAP.LH.HbT.All.mean_f,ephysPowerData.Blank_SAP.LH.HbT.All.mean_normS - ephysPowerData.Blank_SAP.LH.HbT.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(ephysPowerData.SSP_SAP.LH.HbT.All.mean_f,ephysPowerData.SSP_SAP.LH.HbT.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(ephysPowerData.SSP_SAP.LH.HbT.All.mean_f,ephysPowerData.SSP_SAP.LH.HbT.All.mean_normS + ephysPowerData.SSP_SAP.LH.HbT.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(ephysPowerData.SSP_SAP.LH.HbT.All.mean_f,ephysPowerData.SSP_SAP.LH.HbT.All.mean_normS - ephysPowerData.SSP_SAP.LH.HbT.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Contra HbT Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% Gamma power ('All')
subplot(2,6,8)
loglog(ephysPowerData.Blank_SAP.LH.gammaBandPower.All.mean_f,ephysPowerData.Blank_SAP.LH.gammaBandPower.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(ephysPowerData.Blank_SAP.LH.gammaBandPower.All.mean_f,ephysPowerData.Blank_SAP.LH.gammaBandPower.All.mean_normS + ephysPowerData.Blank_SAP.LH.gammaBandPower.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(ephysPowerData.Blank_SAP.LH.gammaBandPower.All.mean_f,ephysPowerData.Blank_SAP.LH.gammaBandPower.All.mean_normS - ephysPowerData.Blank_SAP.LH.gammaBandPower.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(ephysPowerData.SSP_SAP.LH.gammaBandPower.All.mean_f,ephysPowerData.SSP_SAP.LH.gammaBandPower.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(ephysPowerData.SSP_SAP.LH.gammaBandPower.All.mean_f,ephysPowerData.SSP_SAP.LH.gammaBandPower.All.mean_normS + ephysPowerData.SSP_SAP.LH.gammaBandPower.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(ephysPowerData.SSP_SAP.LH.gammaBandPower.All.mean_f,ephysPowerData.SSP_SAP.LH.gammaBandPower.All.mean_normS - ephysPowerData.SSP_SAP.LH.gammaBandPower.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Contra Gamma Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% HbT power ('All')
subplot(2,6,9)
loglog(gcampPowerData.Blank_SAP.LH.HbT.All.mean_f,gcampPowerData.Blank_SAP.LH.HbT.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(gcampPowerData.Blank_SAP.LH.HbT.All.mean_f,gcampPowerData.Blank_SAP.LH.HbT.All.mean_normS + gcampPowerData.Blank_SAP.LH.HbT.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.Blank_SAP.LH.HbT.All.mean_f,gcampPowerData.Blank_SAP.LH.HbT.All.mean_normS - gcampPowerData.Blank_SAP.LH.HbT.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.LH.HbT.All.mean_f,gcampPowerData.SSP_SAP.LH.HbT.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(gcampPowerData.SSP_SAP.LH.HbT.All.mean_f,gcampPowerData.SSP_SAP.LH.HbT.All.mean_normS + gcampPowerData.SSP_SAP.LH.HbT.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.LH.HbT.All.mean_f,gcampPowerData.SSP_SAP.LH.HbT.All.mean_normS - gcampPowerData.SSP_SAP.LH.HbT.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Contra HbT Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% HbO power ('All')
subplot(2,6,10)
loglog(gcampPowerData.Blank_SAP.LH.HbO.All.mean_f,gcampPowerData.Blank_SAP.LH.HbO.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(gcampPowerData.Blank_SAP.LH.HbO.All.mean_f,gcampPowerData.Blank_SAP.LH.HbO.All.mean_normS + gcampPowerData.Blank_SAP.LH.HbO.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.Blank_SAP.LH.HbO.All.mean_f,gcampPowerData.Blank_SAP.LH.HbO.All.mean_normS - gcampPowerData.Blank_SAP.LH.HbO.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.LH.HbO.All.mean_f,gcampPowerData.SSP_SAP.LH.HbO.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(gcampPowerData.SSP_SAP.LH.HbO.All.mean_f,gcampPowerData.SSP_SAP.LH.HbO.All.mean_normS + gcampPowerData.SSP_SAP.LH.HbO.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.LH.HbO.All.mean_f,gcampPowerData.SSP_SAP.LH.HbO.All.mean_normS - gcampPowerData.SSP_SAP.LH.HbO.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Contra HbO Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% HbR power ('All')
subplot(2,6,11)
loglog(gcampPowerData.Blank_SAP.LH.HbR.All.mean_f,gcampPowerData.Blank_SAP.LH.HbR.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(gcampPowerData.Blank_SAP.LH.HbR.All.mean_f,gcampPowerData.Blank_SAP.LH.HbR.All.mean_normS + gcampPowerData.Blank_SAP.LH.HbR.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.Blank_SAP.LH.HbR.All.mean_f,gcampPowerData.Blank_SAP.LH.HbR.All.mean_normS - gcampPowerData.Blank_SAP.LH.HbR.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.LH.HbR.All.mean_f,gcampPowerData.SSP_SAP.LH.HbR.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(gcampPowerData.SSP_SAP.LH.HbR.All.mean_f,gcampPowerData.SSP_SAP.LH.HbR.All.mean_normS + gcampPowerData.SSP_SAP.LH.HbR.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.LH.HbR.All.mean_f,gcampPowerData.SSP_SAP.LH.HbR.All.mean_normS - gcampPowerData.SSP_SAP.LH.HbR.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Contra HbR Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

% GCaMP power ('All')
subplot(2,6,12)
loglog(gcampPowerData.Blank_SAP.LH.GCaMP.All.mean_f,gcampPowerData.Blank_SAP.LH.GCaMP.All.mean_normS,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(gcampPowerData.Blank_SAP.LH.GCaMP.All.mean_f,gcampPowerData.Blank_SAP.LH.GCaMP.All.mean_normS + gcampPowerData.Blank_SAP.LH.GCaMP.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.Blank_SAP.LH.GCaMP.All.mean_f,gcampPowerData.Blank_SAP.LH.GCaMP.All.mean_normS - gcampPowerData.Blank_SAP.LH.GCaMP.All.stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.LH.GCaMP.All.mean_f,gcampPowerData.SSP_SAP.LH.GCaMP.All.mean_normS,'color',colors('electric purple'),'LineWidth',2);
loglog(gcampPowerData.SSP_SAP.LH.GCaMP.All.mean_f,gcampPowerData.SSP_SAP.LH.GCaMP.All.mean_normS + gcampPowerData.SSP_SAP.LH.GCaMP.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
loglog(gcampPowerData.SSP_SAP.LH.GCaMP.All.mean_f,gcampPowerData.SSP_SAP.LH.GCaMP.All.mean_normS - gcampPowerData.SSP_SAP.LH.GCaMP.All.stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
title('Contra GCaMP Power')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.01,0.5])
set(gca,'box','off')
axis square

%% save figure and stats
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figs & Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS4,[dirpath 'FigS4']);
    % statistical diary
    diaryFile = [dirpath 'FigS4_Stats.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on

    % Ipsi HbT power spectra (ephys)
    disp('======================================================================================================================')
    disp('Ipsi HbT power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(EphysPowerStats.RH.HbT.All.Stats)

    % Ipsi gamma-band power spectra (ephys)
    disp('======================================================================================================================')
    disp('Ipsi HbT power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(EphysPowerStats.RH.gammaBandPower.All.Stats)

    % Ipsi HbT power spectra (GCaMP)
    disp('======================================================================================================================')
    disp('Ipsi HbT power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(GCaMPPowerStats.RH.HbT.All.Stats)

    % Ipsi HbO power spectra (GCaMP)
    disp('======================================================================================================================')
    disp('Ipsi HbO power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(GCaMPPowerStats.RH.HbO.All.Stats)

    % Ipsi HbR power spectra (GCaMP)
    disp('======================================================================================================================')
    disp('Ipsi HbR power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(GCaMPPowerStats.RH.HbR.All.Stats)

    % Ipsi GCaMP power spectra (GCaMP)
    disp('======================================================================================================================')
    disp('Ipsi GCaMP power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(GCaMPPowerStats.RH.GCaMP.All.Stats)

     % Contra HbT power spectra (ephys)
    disp('======================================================================================================================')
    disp('Contra HbT power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(EphysPowerStats.LH.HbT.All.Stats)

    % Contra gamma-band power spectra (ephys)
    disp('======================================================================================================================')
    disp('Contra HbT power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(EphysPowerStats.LH.gammaBandPower.All.Stats)

    % Contra HbT power spectra (GCaMP)
    disp('======================================================================================================================')
    disp('Contra HbT power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(GCaMPPowerStats.LH.HbT.All.Stats)

    % Contra HbO power spectra (GCaMP)
    disp('======================================================================================================================')
    disp('Contra HbO power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(GCaMPPowerStats.LH.HbO.All.Stats)

    % Contra HbR power spectra (GCaMP)
    disp('======================================================================================================================')
    disp('Contra HbR power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(GCaMPPowerStats.LH.HbR.All.Stats)

    % Contra GCaMP power spectra (GCaMP)
    disp('======================================================================================================================')
    disp('Contra GCaMP power spectra (ephys), n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(GCaMPPowerStats.LH.GCaMP.All.Stats)

    diary off
end
