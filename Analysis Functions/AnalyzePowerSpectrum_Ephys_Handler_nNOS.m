function [] = AnalyzePowerSpectrum_Ephys_Handler(rootFolder,delim,runFromStart)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
if runFromStart == true
    Results_PowerSpec_Ephys.Naive = [];
    Results_PowerSpec_Ephys.Blank_SAP = [];
    Results_PowerSpec_Ephys.SSP_SAP = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_PowerSpec_Ephys.mat','file') == 2
        load('Results_PowerSpec_Ephys.mat','-mat')
    else
        Results_PowerSpec_Ephys.Naive = [];
        Results_PowerSpec_Ephys.Blank_SAP = [];
        Results_PowerSpec_Ephys.SSP_SAP = [];
    end
end
cd([rootFolder delim 'Data']);
groups = {'Naive','Blank_SAP','SSP_SAP'};
set = 'Ephys';
% determine waitbar length
waitBarLength = 0;
for aa = 1:length(groups)
    folderList = dir([groups{1,aa} delim set]);
    folderList = folderList(~startsWith({folderList.name}, '.'));
    folderAnimalIDs = {folderList.name};
    waitBarLength = waitBarLength + length(folderAnimalIDs);
end
% run analysis for each animal in the group
cc = 1;
multiWaitbar('Analyzing power spectrum for Ephys',0,'Color','P');
for aa = 1:length(groups)
    folderList = dir([groups{1,aa} delim set]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_PowerSpec_Ephys.(groups{1,aa}),animalIDs{1,bb}) == false
            [Results_PowerSpec_Ephys] = AnalyzePowerSpectrum_Ephys(animalIDs{1,bb},groups{1,aa},set,rootFolder,delim,Results_PowerSpec_Ephys);
        end
        multiWaitbar('Analyzing power spectrum for Ephys','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end
multiWaitbar('close all');