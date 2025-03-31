function [] = Fig7_nNOS(rootFolder,saveState,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

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

%% 2P variance
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Diameter_2P';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
variables = {'avg','p2p','vari'};
fs = 5;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Diameter_2P.(group));
    diameterData.(group).dummCheck = 1;
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Diameter_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            if isfield(diameterData.(group),'avg') == false
                diameterData.(group).avg = [];
                diameterData.(group).p2p = [];
                diameterData.(group).vari = [];
                diameterData.(group).group = {};
                diameterData.(group).animalID = {};
            end
            animalVar = [];
            animalP2P = [];
            if isfield(Results_Diameter_2P.(group).(animalID).(vID),'Rest') == true
                for ff = 1:length(Results_Diameter_2P.(group).(animalID).(vID).Rest.indEvents)
                    dataArray = Results_Diameter_2P.(group).(animalID).(vID).Rest.indEvents{ff,1}(2*fs:end);
                    animalVar(ff,1) = var(dataArray);
                    animalP2P(ff,1) = max(dataArray) - min(dataArray);
                end
                diameterData.(group).avg = cat(1,diameterData.(group).avg,mean(Results_Diameter_2P.(group).(animalID).(vID).Rest.mean,'omitnan'));
                diameterData.(group).p2p = cat(1,diameterData.(group).p2p,mean(animalP2P,'omitnan'));
                diameterData.(group).vari = cat(1,diameterData.(group).vari,mean(animalVar,'omitnan'));
                diameterData.(group).group = cat(1,diameterData.(group).group,group);
                diameterData.(group).animalID = cat(1,diameterData.(group).animalID,animalID);
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:length(variables)
        variable = variables{1,ee};
        diameterData.(group).(['mean_' variable]) = mean(diameterData.(group).(variable),1,'omitnan');
        diameterData.(group).(['std_' variable]) = std(diameterData.(group).(variable),1,'omitnan');
    end
end

%% HbT variance statistics
blankHbTVarData = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari);
sspHbTVarData = cat(1,iosSigData.SSP_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,pulseSigData.SSP_SAP.vari);
restHbTVarStats.tableSize = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.vari,iosSigData.SSP_SAP.RH.HbT.Rest.vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari,pulseSigData.SSP_SAP.vari);
restHbTVarStats.Table = table('Size',[size(restHbTVarStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','Variance'});
restHbTVarStats.Table.Mouse = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.animalID,iosSigData.SSP_SAP.RH.HbT.Rest.animalID,gcampSigdata.Blank_SAP.RH.HbT.Rest.animalID,gcampSigdata.SSP_SAP.RH.HbT.Rest.animalID,pulseSigData.Blank_SAP.animalID,pulseSigData.SSP_SAP.animalID);
restHbTVarStats.Table.Group = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.group,iosSigData.SSP_SAP.RH.HbT.Rest.group,gcampSigdata.Blank_SAP.RH.HbT.Rest.group,gcampSigdata.SSP_SAP.RH.HbT.Rest.group,pulseSigData.Blank_SAP.group,pulseSigData.SSP_SAP.group);
restHbTVarStats.Table.Variance = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.vari,iosSigData.SSP_SAP.RH.HbT.Rest.vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari,pulseSigData.SSP_SAP.vari);
restHbTVarStats.FitFormula = 'Variance ~ 1 + Group + (1|Mouse)';
restHbTVarStats.Stats = fitglme(restHbTVarStats.Table,restHbTVarStats.FitFormula);

%% Diameter variance statistics
blankDiameterVarData = cat(1,diameterData.Blank_SAP.vari);
sspDiameterVarData = cat(1,diameterData.SSP_SAP.vari);
restDiameterVarStats.tableSize = cat(1,diameterData.Blank_SAP.vari,diameterData.SSP_SAP.vari);
restDiameterVarStats.Table = table('Size',[size(restDiameterVarStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','Variance'});
restDiameterVarStats.Table.Mouse = cat(1,diameterData.Blank_SAP.animalID,diameterData.SSP_SAP.animalID);
restDiameterVarStats.Table.Group = cat(1,diameterData.Blank_SAP.group,diameterData.SSP_SAP.group);
% restDiameterVarStats.Table.Variance = cat(1,diameterData.Blank_SAP.vari,diameterData.SSP_SAP.vari);
restDiameterVarStats.Table.Variance = cat(1,diameterData.Blank_SAP.vari,sspDiameterVarData);
restDiameterVarStats.FitFormula = 'Variance ~ 1 + Group + (1|Mouse)';
restDiameterVarStats.Stats = fitglme(restDiameterVarStats.Table,restDiameterVarStats.FitFormula);

%% GCaMP variance statistics
blankGCaMPVarData = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.vari);
sspGCaMPVarData = cat(1,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.vari);
restGCaMPVarStats.tableSize = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.vari,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.vari);
restGCaMPVarStats.Table = table('Size',[size(restGCaMPVarStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','Variance'});
restGCaMPVarStats.Table.Mouse = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.animalID,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.animalID);
restGCaMPVarStats.Table.Group = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.group,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.group);
restGCaMPVarStats.Table.Variance = cat(1,gcampSigdata.Blank_SAP.RH.GCaMP.Rest.vari,gcampSigdata.SSP_SAP.RH.GCaMP.Rest.vari);
restGCaMPVarStats.FitFormula = 'Variance ~ 1 + Group + (1|Mouse)';
restGCaMPVarStats.Stats = fitglme(restGCaMPVarStats.Table,restGCaMPVarStats.FitFormula);

%% arousal-state hemodnyamics [Ephys]
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
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
                ephysData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(ephysData.(group).(hemisphere).(dataType),behavior) == false
                        ephysData.(group).(hemisphere).(dataType).(behavior).avg = [];
                        ephysData.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        ephysData.(group).(hemisphere).(dataType).(behavior).vari = [];
                        ephysData.(group).(hemisphere).(dataType).(behavior).group = {};
                        ephysData.(group).(hemisphere).(dataType).(behavior).animalID = {};
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
                    ephysData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).HbT));
                    ephysData.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    ephysData.(group).(hemisphere).(dataType).(behavior).vari = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    ephysData.(group).(hemisphere).(dataType).(behavior).group = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).group,group);
                    ephysData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,ephysData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
                    ephysData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(ephysData.(group).(hemisphere).(dataType).(behavior).(variable),1,'omitnan');
                    ephysData.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(ephysData.(group).(hemisphere).(dataType).(behavior).(variable),1,'omitnan');
                end
            end
        end
    end
end
% statistics - ttest
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                % statistics - unpaired ttest
                [EphysStats.(hemisphere).(dataType).(behavior).(variable).h,EphysStats.(hemisphere).(dataType).(behavior).(variable).p] = ttest2(ephysData.Blank_SAP.(hemisphere).(dataType).(behavior).(variable),ephysData.SSP_SAP.(hemisphere).(dataType).(behavior).(variable));
            end
        end
    end
end

%% arousal-state hemodnyamics [GCaMP]
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
                gcampData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(gcampData.(group).(hemisphere).(dataType),behavior) == false
                        gcampData.(group).(hemisphere).(dataType).(behavior).avg = [];
                        gcampData.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        gcampData.(group).(hemisphere).(dataType).(behavior).vari = [];
                        gcampData.(group).(hemisphere).(dataType).(behavior).group = {};
                        gcampData.(group).(hemisphere).(dataType).(behavior).animalID = {};
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
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    try
                        gcampData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean));
                    catch
                        gcampData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).avg,mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean)));
                    end
                    gcampData.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    gcampData.(group).(hemisphere).(dataType).(behavior).vari = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    gcampData.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).group,group);
                    gcampData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
                    gcampData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(gcampData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    gcampData.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(gcampData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                end
            end
        end
    end
end
% statistics - ttest
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(variables)
                variable = variables{1,dd};
                % statistics - unpaired ttest
                [GCaMPStats.(hemisphere).(dataType).(behavior).(variable).h,GCaMPStats.(hemisphere).(dataType).(behavior).(variable).p] = ttest2(gcampData.Blank_SAP.(hemisphere).(dataType).(behavior).(variable),gcampData.SSP_SAP.(hemisphere).(dataType).(behavior).(variable));
            end
        end
    end
end

%% two photo diameter
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Diameter_2P';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
behaviors = {'Rest','Whisk'};
variables = {'data'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Diameter_2P.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Diameter_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                diameterData.(group).(behavior).dummCheck = 1;
                for ee = 1:length(variables)
                    if isfield(diameterData.(group).(behavior),(variables{1,ee})) == false
                        diameterData.(group).(behavior).(variables{1,ee}) = [];
                    end
                end
                if isfield(Results_Diameter_2P.(group).(animalID).(vID),behavior) == true
                    diameterData.(group).(behavior).data = cat(1,diameterData.(group).(behavior).data,mean(Results_Diameter_2P.(group).(animalID).(vID).(behavior).mean));
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        diameterData.(group).(behavior).meanData = mean(diameterData.(group).(behavior).data,1);
        diameterData.(group).(behavior).stdData = std(diameterData.(group).(behavior).data,0,1);
    end
end
% statistics - ttest
for cc = 1:length(behaviors)
    behavior = behaviors{1,cc};
    [DiameterStats.(behavior).h,DiameterStats.(behavior).p] = ttest2(diameterData.Blank_SAP.(behavior).data,diameterData.SSP_SAP.(behavior).data);
end

%% two photon isoflurane shift
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Baseline_2P';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
variables = {'diameter','baseline'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Baseline_2P.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Baseline_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            diameterShiftData.(group).dummCheck = 1;
            for dd = 1:length(variables)
                if isfield(diameterShiftData.(group),(variables{1,dd})) == false
                    diameterShiftData.(group).(variables{1,dd}) = [];
                end
            end
            diameterShiftData.(group).diameter = cat(1,diameterShiftData.(group).diameter,((Results_Baseline_2P.(group).(animalID).(vID).diameter - Results_Baseline_2P.(group).(animalID).(vID).baseline)/Results_Baseline_2P.(group).(animalID).(vID).baseline)*100);
            diameterShiftData.(group).baseline = cat(1,diameterShiftData.(group).baseline,Results_Baseline_2P.(group).(animalID).(vID).baseline);
        end
    end
end
% mean/std
for ee = 1:length(groups)
    group = groups{1,ee};
    diameterShiftData.(group).meanDiameter = mean(diameterShiftData.(group).diameter,1);
    diameterShiftData.(group).stdDiameter = std(diameterShiftData.(group).diameter,0,1);
    diameterShiftData.(group).meanBaseline = mean(diameterShiftData.(group).baseline,1);
    diameterShiftData.(group).stdBaseline = std(diameterShiftData.(group).baseline,0,1);
end
[TwoPIsoStats.h,TwoPIsoStats.p] = ttest2(diameterShiftData.Blank_SAP.diameter,diameterShiftData.SSP_SAP.diameter);

%% Figure 7
Fig7 = figure('Name','Fig. 7');

% HbT rest variance
subplot(2,2,1)
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

% Diameter rest variance
subplot(2,2,2)
xInds = ones(1,length(blankDiameterVarData));
scatter(xInds*1,blankDiameterVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('black'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(blankDiameterVarData,'omitnan'),std(blankDiameterVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(sspDiameterVarData));
scatter(xInds*2,sspDiameterVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('dark candy apple red'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,mean(sspDiameterVarData,'omitnan'),std(sspDiameterVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[Diameter]^2 (\muM)')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
axis tight
xlim([0,3]);
ylim([-5,60])


subplot(2,2,3);
hold on
xInds = ones(1,length(ephysData.Blank_SAP.RH.HbT.Rest.avg));
s2 = scatter(xInds*1,ephysData.Blank_SAP.RH.HbT.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('black'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(1,ephysData.Blank_SAP.RH.HbT.Rest.mean_avg,ephysData.Blank_SAP.RH.HbT.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysData.SSP_SAP.RH.HbT.Rest.avg));
s3 = scatter(xInds*2,ephysData.SSP_SAP.RH.HbT.Rest.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(2,ephysData.SSP_SAP.RH.HbT.Rest.mean_avg,ephysData.SSP_SAP.RH.HbT.Rest.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

xInds = ones(1,length(ephysData.Blank_SAP.RH.HbT.NREM.avg));
scatter(xInds*4,ephysData.Blank_SAP.RH.HbT.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('black'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(4,ephysData.Blank_SAP.RH.HbT.NREM.mean_avg,ephysData.Blank_SAP.RH.HbT.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysData.SSP_SAP.RH.HbT.NREM.avg));
scatter(xInds*5,ephysData.SSP_SAP.RH.HbT.NREM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(5,ephysData.SSP_SAP.RH.HbT.NREM.mean_avg,ephysData.SSP_SAP.RH.HbT.NREM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

xInds = ones(1,length(ephysData.Blank_SAP.RH.HbT.REM.avg));
scatter(xInds*7,ephysData.Blank_SAP.RH.HbT.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('black'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(7,ephysData.Blank_SAP.RH.HbT.REM.mean_avg,ephysData.Blank_SAP.RH.HbT.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysData.SSP_SAP.RH.HbT.REM.avg));
scatter(xInds*8,ephysData.SSP_SAP.RH.HbT.REM.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(8,ephysData.SSP_SAP.RH.HbT.REM.mean_avg,ephysData.SSP_SAP.RH.HbT.REM.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
xInds = ones(1,length(ephysData.Blank_SAP.RH.HbT.Iso.avg));
scatter(xInds*10,ephysData.Blank_SAP.RH.HbT.Iso.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('black'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(10,ephysData.Blank_SAP.RH.HbT.Iso.mean_avg,ephysData.Blank_SAP.RH.HbT.Iso.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(ephysData.SSP_SAP.RH.HbT.Iso.avg));
scatter(xInds*11,ephysData.SSP_SAP.RH.HbT.Iso.avg,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(11,ephysData.SSP_SAP.RH.HbT.Iso.mean_avg,ephysData.SSP_SAP.RH.HbT.Iso.std_avg,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('\Delta[HbT] (\muM)')
xlim([0,12])
legend([s2,s3],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
set(gca,'xtick',[])
axis square

subplot(2,2,4)
xInds = ones(1,length(diameterShiftData.Blank_SAP.diameter));
scatter(xInds*1,diameterShiftData.Blank_SAP.diameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('black'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,diameterShiftData.Blank_SAP.meanDiameter,diameterShiftData.Blank_SAP.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(diameterShiftData.SSP_SAP.diameter));
scatter(xInds*2,diameterShiftData.SSP_SAP.diameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('dark candy apple red'),'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,diameterShiftData.SSP_SAP.meanDiameter,diameterShiftData.SSP_SAP.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\DeltaD/D (%)')
set(gca,'xtick',[1,1.2,2,2.2])
set(gca,'xticklabel',{['Avg Baseline: ' num2str(round(diameterShiftData.Blank_SAP.meanBaseline,1)) ' \muM'],['n = ' num2str(length(diameterShiftData.Blank_SAP.diameter)) ' arterioles'],['Avg Baseline: ' num2str(round(diameterShiftData.SSP_SAP.meanBaseline,1)) ' \muM'],['n = ' num2str(length(diameterShiftData.SSP_SAP.diameter)) ' arterioles']})
xtickangle(45)
axis square
xlim([0.5,2.5])
set(gca,'box','off')
axis square

%% Save figure and stats
if saveState == true
    dirpath = [rootFolder delim 'MATLAB Figs/Stats' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig7,[dirpath 'Fig4']);
    diaryFile = [dirpath 'Fig4_Stats.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    % statistical diary
    diary(diaryFile)
    diary on

    comparisons = 1;
    alphaA = 0.05/comparisons;
    alphaB = 0.01/comparisons;
    alphaC = 0.001/comparisons;

    % IOS [HbT] resting variance
    disp('======================================================================================================================')
    disp('IOS resting variance: Blank (N = 24, 12M/12F); SSP (N = 24, 11M/13F); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(blankHbTVarData,'omitnan')) ' +/- ' num2str(std(blankHbTVarData,0,1,'omitnan')./sqrt(size(blankHbTVarData,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(sspHbTVarData,'omitnan')) ' +/- ' num2str(std(sspHbTVarData,0,1,'omitnan')./sqrt(size(sspHbTVarData,1)))]); disp(' ')
    disp('GLME statistics for resting HbT variance')
    disp(restHbTVarStats.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

     % [HbT] during each arousal-state
    disp('======================================================================================================================')
    disp('[HbT] during each arousal-state (ephys), n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(ephysData.Blank_SAP.RH.HbT.Rest.mean_avg) ' +/- ' num2str(ephysData.Blank_SAP.RH.HbT.Rest.std_avg)]); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(ephysData.SSP_SAP.RH.HbT.Rest.mean_avg) ' +/- ' num2str(ephysData.SSP_SAP.RH.HbT.Rest.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(EphysStats.RH.HbT.Rest.avg.p)]); disp(' ')
    disp(['Blank-SAP NREM: ' num2str(ephysData.Blank_SAP.RH.HbT.NREM.mean_avg) ' +/- ' num2str(ephysData.Blank_SAP.RH.HbT.NREM.std_avg)]); disp(' ')
    disp(['SSP-SAP NREM: ' num2str(ephysData.SSP_SAP.RH.HbT.NREM.mean_avg) ' +/- ' num2str(ephysData.SSP_SAP.RH.HbT.NREM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP NREM ttest p = ' num2str(EphysStats.RH.HbT.NREM.avg.p)]); disp(' ')
    disp(['Blank-SAP REM: ' num2str(ephysData.Blank_SAP.RH.HbT.REM.mean_avg) ' +/- ' num2str(ephysData.Blank_SAP.RH.HbT.REM.std_avg)]); disp(' ')
    disp(['SSP-SAP REM: ' num2str(ephysData.SSP_SAP.RH.HbT.REM.mean_avg) ' +/- ' num2str(ephysData.SSP_SAP.RH.HbT.REM.std_avg)]); disp(' ')
    disp(['Blank vs. SAP REM ttest p = ' num2str(EphysStats.RH.HbT.REM.avg.p)]); disp(' ')
    disp(['Blank-SAP Iso: ' num2str(ephysData.Blank_SAP.RH.HbT.Iso.mean_avg) ' +/- ' num2str(ephysData.Blank_SAP.RH.HbT.Iso.std_avg)]); disp(' ')
    disp(['SSP-SAP Iso: ' num2str(ephysData.SSP_SAP.RH.HbT.Iso.mean_avg) ' +/- ' num2str(ephysData.SSP_SAP.RH.HbT.Iso.std_avg)]); disp(' ')
    disp(['Blank vs. SAP Iso ttest p = ' num2str(EphysStats.RH.HbT.Iso.avg.p)]); disp(' ')

    % Diameter resting variance
    disp('======================================================================================================================')
    disp('Diameter resting variance: Blank (N = 9, 5M/4F, n = 70); SSP (N = 7, 2M/5F, n = 65); mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(blankDiameterVarData,'omitnan')) ' +/- ' num2str(std(blankDiameterVarData,0,1,'omitnan')./sqrt(size(blankDiameterVarData,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(sspDiameterVarData,'omitnan')) ' +/- ' num2str(std(sspDiameterVarData,0,1,'omitnan')./sqrt(size(sspDiameterVarData,1)))]); disp(' ')
    disp('GLME statistics for resting HbT variance')
    disp(restDiameterVarStats.Stats)
    disp(['*p < ' num2str(alphaA) ' **p < ' num2str(alphaB) ' ***p < ' num2str(alphaC)]);

        % isoflurane shift 2P diameter
    disp('======================================================================================================================')
    disp('isoflurane diameter shift, n = 9 mice per group, mean +/- StD'); disp(' ')
    disp(['Blank-SAP Rest: ' num2str(diameterShiftData.Blank_SAP.meanDiameter) ' +/- ' num2str(diameterShiftData.Blank_SAP.stdDiameter) ' baseline diameter: ' num2str(diameterShiftData.Blank_SAP.meanBaseline) '+/-' num2str(diameterShiftData.Blank_SAP.stdBaseline) ' \muM']); disp(' ')
    disp(['SSP-SAP Rest: ' num2str(diameterShiftData.SSP_SAP.meanDiameter) ' +/- ' num2str(diameterShiftData.SSP_SAP.stdDiameter) ' baseline diameter: ' num2str(diameterShiftData.SSP_SAP.meanBaseline) '+/-' num2str(diameterShiftData.SSP_SAP.stdBaseline) ' /muM']); disp(' ')
    disp(['Blank vs. SAP Rest ttest p = ' num2str(TwoPIsoStats.p)]); disp(' ')


    diary off
end