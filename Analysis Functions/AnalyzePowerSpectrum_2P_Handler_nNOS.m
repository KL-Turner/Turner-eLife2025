function [] = AnalyzePowerSpectrum_2P_Handler(rootFolder,delim,runFromStart)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
if runFromStart == true
        Results_PowerSpec_2P.Blank_SAP = [];
        Results_PowerSpec_2P.SSP_SAP = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_PowerSpec_2P.mat','file') == 2
        load('Results_PowerSpec_2P.mat','-mat')
    else
        Results_PowerSpec_2P.Blank_SAP = [];
        Results_PowerSpec_2P.SSP_SAP = [];
    end
end
cd([rootFolder delim 'Data']);
groups = {'Blank_SAP','SSP_SAP'};
set = '2P';
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
multiWaitbar('Analyzing power spectrum for 2P',0,'Color','P');
for aa = 1:length(groups)
    folderList = dir([groups{1,aa} delim set]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_PowerSpec_2P.(groups{1,aa}),animalIDs{1,bb}) == false
            [Results_PowerSpec_2P] = AnalyzePowerSpectrum_2P(animalIDs{1,bb},groups{1,aa},set,rootFolder,delim,Results_PowerSpec_2P);
        end
        multiWaitbar('Analyzing power spectrum for 2P',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end
multiWaitbar('close all');