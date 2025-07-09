function [] = FigS7_eLife2025(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% LFP power
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PowerSpec_LFP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Alert','Asleep','All'};
variables = {'S','f','deltaS'};
dimensions = [2,1,1];
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PowerSpec_LFP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                lfpData.(group).(hemisphere).(behavior).dummCheck = 1;
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    dimension = dimensions(ee);
                    if isfield(lfpData.(group).(hemisphere).(behavior),(variable)) == false
                        lfpData.(group).(hemisphere).(behavior).(variable) = [];
                        lfpData.(group).(hemisphere).(behavior).group = {};
                        lfpData.(group).(hemisphere).(behavior).animalID = {};
                        lfpData.(group).(hemisphere).(behavior).hemisphere = {};
                        lfpData.(group).(hemisphere).(behavior).behavior = {};
                    end
                    % pull data if field isn't empty
                    if isempty(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S) == false
                        if strcmp(variable,'deltaS') == true
                            index = find(round(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).f,2) == 4);
                            deltaIndex = index(end);
                            lfpData.(group).(hemisphere).(behavior).(variable) = cat(dimension,lfpData.(group).(hemisphere).(behavior).(variable),mean(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S(1:deltaIndex)));
                            % for stats
                            lfpData.(group).(hemisphere).(behavior).group = cat(1,lfpData.(group).(hemisphere).(behavior).group,group);
                            lfpData.(group).(hemisphere).(behavior).animalID = cat(1,lfpData.(group).(hemisphere).(behavior).animalID,animalID);
                            lfpData.(group).(hemisphere).(behavior).hemisphere = cat(1,lfpData.(group).(hemisphere).(behavior).hemisphere,hemisphere);
                            lfpData.(group).(hemisphere).(behavior).behavior = cat(1,lfpData.(group).(hemisphere).(behavior).behavior,behavior);
                        else
                            lfpData.(group).(hemisphere).(behavior).(variable) = cat(dimension,lfpData.(group).(hemisphere).(behavior).(variable),Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).(variable));
                        end
                    end
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
            for dd = 1:length(variables)
                variable = variables{1,dd};
                dimension = dimensions(dd);
                lfpData.(group).(hemisphere).(behavior).(['mean_' variable]) = mean(lfpData.(group).(hemisphere).(behavior).(variable),dimension,'omitnan');
                lfpData.(group).(hemisphere).(behavior).(['stdErr_' variable]) = std(lfpData.(group).(hemisphere).(behavior).(variable),0,dimension,'omitnan')./sqrt(size(lfpData.(group).(hemisphere).(behavior).(variable),dimension));
            end
        end
    end
end

%% statistics - generalized linear mixed effects model
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        lfpStats.(hemisphere).(behavior).tableSize = cat(1,lfpData.Blank_SAP.(hemisphere).(behavior).deltaS,lfpData.SSP_SAP.(hemisphere).(behavior).deltaS,lfpData.Naive.(hemisphere).(behavior).deltaS);
        lfpStats.(hemisphere).(behavior).Table = table('Size',[size(lfpStats.(hemisphere).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'group','animalID','behavior','deltaS'});
        lfpStats.(hemisphere).(behavior).Table.group = cat(1,lfpData.Blank_SAP.(hemisphere).(behavior).group,lfpData.SSP_SAP.(hemisphere).(behavior).group,lfpData.Naive.(hemisphere).(behavior).group);
        lfpStats.(hemisphere).(behavior).Table.animalID = cat(1,lfpData.Blank_SAP.(hemisphere).(behavior).animalID,lfpData.SSP_SAP.(hemisphere).(behavior).animalID,lfpData.Naive.(hemisphere).(behavior).animalID);
        lfpStats.(hemisphere).(behavior).Table.behavior = cat(1,lfpData.Blank_SAP.(hemisphere).(behavior).behavior,lfpData.SSP_SAP.(hemisphere).(behavior).behavior,lfpData.Naive.(hemisphere).(behavior).behavior);
        lfpStats.(hemisphere).(behavior).Table.deltaS = cat(1,lfpData.Blank_SAP.(hemisphere).(behavior).deltaS,lfpData.SSP_SAP.(hemisphere).(behavior).deltaS,lfpData.Naive.(hemisphere).(behavior).deltaS);
        lfpStats.(hemisphere).(behavior).FitFormula = 'deltaS ~ 1 + group + behavior + (1|animalID)';
        lfpStats.(hemisphere).(behavior).Stats = fitglme(lfpStats.(hemisphere).(behavior).Table,lfpStats.(hemisphere).(behavior).FitFormula);
    end
end

%% IOS variance signals
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','Whisk','Stim','NREM','REM','Iso'};
dataTypes = {'HbT'};
variables = {'avg','p2p','vari'};
fs = 30;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                iosSigData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(iosSigData.(group).(hemisphere).(dataType),behavior) == false
                        iosSigData.(group).(hemisphere).(dataType).(behavior).avg = [];
                        iosSigData.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        iosSigData.(group).(hemisphere).(dataType).(behavior).vari = [];
                        iosSigData.(group).(hemisphere).(dataType).(behavior).group = {};
                        iosSigData.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    animalVar = [];
                    animalP2P = [];
                    for ff = 1:length(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT)
                        if strcmp(behavior,'Rest') == true
                            dataArray = Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1}(2*fs:end);
                        elseif strcmp(behavior,'Stim') == true
                            dataArray = Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1} - mean(cell2mat(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT),1);
                        else
                            dataArray = Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1};
                        end
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    iosSigData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).HbT));
                    iosSigData.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    iosSigData.(group).(hemisphere).(dataType).(behavior).vari = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    iosSigData.(group).(hemisphere).(dataType).(behavior).group = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).group,group);
                    iosSigData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    iosSigData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(iosSigData.(group).(hemisphere).(dataType).(behavior).(variable),1,'omitnan');
                    iosSigData.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(iosSigData.(group).(hemisphere).(dataType).(behavior).(variable),1,'omitnan');
                end
            end
        end
    end
end

%% IOS pulse variance
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_Pulse';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
variables = {'avg','p2p','vari'};
fs = 30;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_Pulse.(group));
    pulseSigData.(group).dummCheck = 1;
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        if isfield(pulseSigData.(group),'avg') == false
            pulseSigData.(group).avg = [];
            pulseSigData.(group).p2p = [];
            pulseSigData.(group).vari = [];
            pulseSigData.(group).group = {};
            pulseSigData.(group).animalID = {};
        end
        animalVar = [];
        animalP2P = [];
        if isfield(Results_IntSig_Pulse.(group).(animalID),'Rest') == true
            for ff = 1:length(Results_IntSig_Pulse.(group).(animalID).Rest.indHbT)
                dataArray = Results_IntSig_Pulse.(group).(animalID).Rest.indHbT{ff,1}(2*fs:end);
                animalVar(ff,1) = var(dataArray);
                animalP2P(ff,1) = max(dataArray) - min(dataArray);
            end
            pulseSigData.(group).avg = cat(1,pulseSigData.(group).avg,mean(Results_IntSig_Pulse.(group).(animalID).Rest.HbT,'omitnan'));
            pulseSigData.(group).p2p = cat(1,pulseSigData.(group).p2p,mean(animalP2P,'omitnan'));
            pulseSigData.(group).vari = cat(1,pulseSigData.(group).vari,mean(animalVar,'omitnan'));
            pulseSigData.(group).group = cat(1,pulseSigData.(group).group,group);
            pulseSigData.(group).animalID = cat(1,pulseSigData.(group).animalID,animalID);
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:length(variables)
        variable = variables{1,ee};
        pulseSigData.(group).(['mean_' variable]) = mean(pulseSigData.(group).(variable),1,'omitnan');
        pulseSigData.(group).(['std_' variable]) = std(pulseSigData.(group).(variable),1,'omitnan');
    end
end

%% GCaMP variance signals
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','Whisk','Stim','NREM','REM'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
variables = {'avg','p2p','vari'};
fs = 10;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                gcampSigdata.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(gcampSigdata.(group).(hemisphere).(dataType),behavior) == false
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg = [];
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).vari = [];
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).group = {};
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    animalVar = [];
                    animalP2P = [];
                    for ff = 1:length(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData)
                        if strcmp(behavior,'Rest') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1}(2*fs:end);
                        elseif strcmp(behavior,'Stim') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1} - mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData),1);
                        else
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1};
                        end
                        if strcmp(dataType,'GCaMP') == true
                            dataArray = (dataArray - 1)*100;
                        end
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    try
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean));
                    catch
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg,mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean)));
                    end
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).vari = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).group,group);
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(gcampSigdata.(group).(hemisphere).(dataType).(behavior).(variable),1,'omitnan');
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(gcampSigdata.(group).(hemisphere).(dataType).(behavior).(variable),1,'omitnan');
                end
            end
        end
    end
end

%% HbT variance statistics
blankHbTVarData = cat(1,iosSigData.Blank_SAP.LH.HbT.Rest.vari,gcampSigdata.Blank_SAP.LH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari);
sspHbTVarData = cat(1,iosSigData.SSP_SAP.LH.HbT.Rest.vari,gcampSigdata.SSP_SAP.LH.HbT.Rest.vari,pulseSigData.SSP_SAP.vari);
restHbTVarStats.tableSize = cat(1,iosSigData.Blank_SAP.LH.HbT.Rest.vari,iosSigData.SSP_SAP.LH.HbT.Rest.vari,gcampSigdata.Blank_SAP.LH.HbT.Rest.vari,gcampSigdata.SSP_SAP.LH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari,pulseSigData.SSP_SAP.vari);
restHbTVarStats.Table = table('Size',[size(restHbTVarStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','Variance'});
restHbTVarStats.Table.Mouse = cat(1,iosSigData.Blank_SAP.LH.HbT.Rest.animalID,iosSigData.SSP_SAP.LH.HbT.Rest.animalID,gcampSigdata.Blank_SAP.LH.HbT.Rest.animalID,gcampSigdata.SSP_SAP.LH.HbT.Rest.animalID,pulseSigData.Blank_SAP.animalID,pulseSigData.SSP_SAP.animalID);
restHbTVarStats.Table.Group = cat(1,iosSigData.Blank_SAP.LH.HbT.Rest.group,iosSigData.SSP_SAP.LH.HbT.Rest.group,gcampSigdata.Blank_SAP.LH.HbT.Rest.group,gcampSigdata.SSP_SAP.LH.HbT.Rest.group,pulseSigData.Blank_SAP.group,pulseSigData.SSP_SAP.group);
restHbTVarStats.Table.Variance = cat(1,iosSigData.Blank_SAP.LH.HbT.Rest.vari,iosSigData.SSP_SAP.LH.HbT.Rest.vari,gcampSigdata.Blank_SAP.LH.HbT.Rest.vari,gcampSigdata.SSP_SAP.LH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari,pulseSigData.SSP_SAP.vari);
restHbTVarStats.FitFormula = 'Variance ~ 1 + Group + (1|Mouse)';
restHbTVarStats.Stats = fitglme(restHbTVarStats.Table,restHbTVarStats.FitFormula);

%% IOS variance power spectra
iosPowerSpectra.Blank_SAP.psd = [];
iosPowerSpectra.SSP_SAP.psd = [];
freqRes = 0.01;
commonFreqs = 0:freqRes:5; % Limit to lowest Nyquist (10 Hz group)
fs = 30;
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        dataCellArray = Results_IntSig_Ephys.(group).(animalID).LH.Rest.indHbT;
        nTrials = length(dataCellArray);
        psdMatrix = zeros(length(commonFreqs), nTrials);
        for ii = 1:nTrials
            data = detrend(dataCellArray{ii},'constant');
            [pxx,f] = pwelch(data, [],[],[],fs);
            psdMatrix(:,ii) = interp1(f,pxx,commonFreqs,'linear',0);
        end
        meanPSD = mean(psdMatrix, 2);
        iosPowerSpectra.(group).psd = cat(2,iosPowerSpectra.(group).psd,meanPSD);
    end
end

fs = 10;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        dataCellArray = Results_IntSig_GCaMP.(group).(animalID).LH.HbT.Rest.indData;
        nTrials = length(dataCellArray);
        psdMatrix = zeros(length(commonFreqs),nTrials);
        for ii = 1:nTrials
            data = detrend(dataCellArray{ii},'constant');
            [pxx,f] = pwelch(data,[],[],[],fs);
            psdMatrix(:,ii) = interp1(f,pxx,commonFreqs,'linear',0);
        end
        meanPSD = mean(psdMatrix, 2);
        iosPowerSpectra.(group).psd = cat(2,iosPowerSpectra.(group).psd,meanPSD);
    end
end
iosBlankPSD = mean(iosPowerSpectra.Blank_SAP.psd,2);
iosSPPSD = mean(iosPowerSpectra.SSP_SAP.psd,2);

%% Figure S7
FigS7 = figure('Name','Fig. S7');

% Ephys alert LFP
ax1 = subplot(1,5,1);
loglog(lfpData.Blank_SAP.LH.Alert.mean_f,lfpData.Blank_SAP.LH.Alert.mean_S,'color',colors('black'),'LineWidth',2);
hold on
loglog(lfpData.Blank_SAP.LH.Alert.mean_f,lfpData.Blank_SAP.LH.Alert.mean_S + lfpData.Blank_SAP.LH.Alert.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(lfpData.Blank_SAP.LH.Alert.mean_f,lfpData.Blank_SAP.LH.Alert.mean_S - lfpData.Blank_SAP.LH.Alert.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(lfpData.SSP_SAP.LH.Alert.mean_f,lfpData.SSP_SAP.LH.Alert.mean_S,'color',colors('sapphire'),'LineWidth',2);
loglog(lfpData.SSP_SAP.LH.Alert.mean_f,lfpData.SSP_SAP.LH.Alert.mean_S + lfpData.SSP_SAP.LH.Alert.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
loglog(lfpData.SSP_SAP.LH.Alert.mean_f,lfpData.SSP_SAP.LH.Alert.mean_S - lfpData.SSP_SAP.LH.Alert.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);

% Ephys asleep LFP
ax2 = subplot(1,5,2);
loglog(lfpData.Blank_SAP.LH.Asleep.mean_f,lfpData.Blank_SAP.LH.Asleep.mean_S,'color',colors('black'),'LineWidth',2);
hold on
loglog(lfpData.Blank_SAP.LH.Asleep.mean_f,lfpData.Blank_SAP.LH.Asleep.mean_S + lfpData.Blank_SAP.LH.Asleep.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(lfpData.Blank_SAP.LH.Asleep.mean_f,lfpData.Blank_SAP.LH.Asleep.mean_S - lfpData.Blank_SAP.LH.Asleep.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(lfpData.SSP_SAP.LH.Asleep.mean_f,lfpData.SSP_SAP.LH.Asleep.mean_S,'color',colors('sapphire'),'LineWidth',2);
loglog(lfpData.SSP_SAP.LH.Asleep.mean_f,lfpData.SSP_SAP.LH.Asleep.mean_S + lfpData.SSP_SAP.LH.Asleep.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
loglog(lfpData.SSP_SAP.LH.Asleep.mean_f,lfpData.SSP_SAP.LH.Asleep.mean_S - lfpData.SSP_SAP.LH.Asleep.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);

% Ephys all LFPS
ax3 = subplot(1,5,3);
loglog(lfpData.Blank_SAP.LH.All.mean_f,lfpData.Blank_SAP.LH.All.mean_S,'color',colors('black'),'LineWidth',2);
hold on
loglog(lfpData.Blank_SAP.LH.All.mean_f,lfpData.Blank_SAP.LH.All.mean_S + lfpData.Blank_SAP.LH.All.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(lfpData.Blank_SAP.LH.All.mean_f,lfpData.Blank_SAP.LH.All.mean_S - lfpData.Blank_SAP.LH.All.stdErr_S,'color',colors('black'),'LineWidth',0.25);
loglog(lfpData.SSP_SAP.LH.All.mean_f,lfpData.SSP_SAP.LH.All.mean_S,'color',colors('sapphire'),'LineWidth',2);
loglog(lfpData.SSP_SAP.LH.All.mean_f,lfpData.SSP_SAP.LH.All.mean_S + lfpData.SSP_SAP.LH.All.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
loglog(lfpData.SSP_SAP.LH.All.mean_f,lfpData.SSP_SAP.LH.All.mean_S - lfpData.SSP_SAP.LH.All.stdErr_S,'color',colors('sapphire'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);
linkaxes([ax1,ax2,ax3],'xy')

% HbT rest variance
subplot(1,5,4)
xInds = ones(1,length(blankHbTVarData));
scatter(xInds*1,blankHbTVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('black'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(blankHbTVarData,'omitnan'),std(blankHbTVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(sspHbTVarData));
scatter(xInds*2,sspHbTVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,mean(sspHbTVarData,'omitnan'),std(sspHbTVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[HbT]^2 (\muM)')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
axis tight
xlim([0,3]);
ylim([0,90])

% HbT power spectra
subplot(1,5,5)
semilogy(commonFreqs,iosBlankPSD,'color',colors('black'),'LineWidth',2);
hold on;
semilogy(commonFreqs,iosSPPSD,'color',colors('electric purple'),'LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');
set(gca,'box','off')
axis square
grid on
xlim([0.1,2.4])

%% save figure and stats
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figs & Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS7,[dirpath 'FigS7']);
    % statistical diary
    diaryFile = [dirpath 'FigS7_Stats.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on

    comparisons = 1;
    alphaA = 0.05/comparisons;
    alphaB = 0.01/comparisons;
    alphaC = 0.001/comparisons;

    % LFP delta power (alert)
    disp('======================================================================================================================')
    disp('LFP Alert delta power: Blank (N = 9, 4M/5F); SSP (N = 9, 5M/4F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(lfpData.Blank_SAP.LH.Alert.deltaS,'omitnan')) ' +/- ' num2str(std(lfpData.Blank_SAP.LH.Alert.deltaS,0,1,'omitnan')./sqrt(size(lfpData.Blank_SAP.LH.Alert.deltaS,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(lfpData.SSP_SAP.LH.Alert.deltaS,'omitnan')) ' +/- ' num2str(std(lfpData.SSP_SAP.LH.Alert.deltaS,0,1,'omitnan')./sqrt(size(lfpData.SSP_SAP.LH.Alert.deltaS,1)))]); disp(' ')
    disp('GLME statistics for LFP delta (1:4 Hz)')
    disp(lfpStats.LH.Alert.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % LFP delta power (asleep)
    disp('======================================================================================================================')
    disp('LFP Asleep delta power: Blank (N = 7, 3M/4F); SSP (N = 7, 4M/3F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(lfpData.Blank_SAP.LH.Asleep.deltaS,'omitnan')) ' +/- ' num2str(std(lfpData.Blank_SAP.LH.Asleep.deltaS,0,1,'omitnan')./sqrt(size(lfpData.Blank_SAP.LH.Asleep.deltaS,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(lfpData.SSP_SAP.LH.Asleep.deltaS,'omitnan')) ' +/- ' num2str(std(lfpData.SSP_SAP.LH.Asleep.deltaS,0,1,'omitnan')./sqrt(size(lfpData.SSP_SAP.LH.Asleep.deltaS,1)))]); disp(' ')
    disp('GLME statistics for LFP delta (1:4 Hz)')
    disp(lfpStats.LH.Asleep.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % LFP delta power (all)
    disp('======================================================================================================================')
    disp('LFP All delta power: Blank (N = 9, 4M/5F); SSP (N = 9, 5M/4F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(lfpData.Blank_SAP.LH.All.deltaS,'omitnan')) ' +/- ' num2str(std(lfpData.Blank_SAP.LH.All.deltaS,0,1,'omitnan')./sqrt(size(lfpData.Blank_SAP.LH.All.deltaS,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(lfpData.SSP_SAP.LH.All.deltaS,'omitnan')) ' +/- ' num2str(std(lfpData.SSP_SAP.LH.All.deltaS,0,1,'omitnan')./sqrt(size(lfpData.SSP_SAP.LH.All.deltaS,1)))]); disp(' ')
    disp('GLME statistics for LFP delta (1:4 Hz)')
    disp(lfpStats.LH.All.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

    % IOS [HbT] resting variance
    disp('======================================================================================================================')
    disp('IOS resting variance: Blank (N = 16, 7M/9F); SSP (N = 17, 9M/8F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(blankHbTVarData,'omitnan')) ' +/- ' num2str(std(blankHbTVarData,0,1,'omitnan')./sqrt(size(blankHbTVarData,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(sspHbTVarData,'omitnan')) ' +/- ' num2str(std(sspHbTVarData,0,1,'omitnan')./sqrt(size(sspHbTVarData,1)))]); disp(' ')
    disp('GLME statistics for resting HbT variance')
    disp(restHbTVarStats.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);
end
