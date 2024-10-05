function [] = AnalyzeModelAccuracy_GCaMP_Handler(rootFolder,delim,runFromStart)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
if runFromStart == true
    Results_ModelAccuracy_GCaMP.Blank_SAP = [];
    Results_ModelAccuracy_GCaMP.SSP_SAP = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_ModelAccuracy_GCaMP.mat','file') == 2
        load('Results_ModelAccuracy_GCaMP.mat','-mat')
    else
        Results_ModelAccuracy_GCaMP.Blank_SAP = [];
        Results_ModelAccuracy_GCaMP.SSP_SAP = [];
    end
end
cd([rootFolder delim 'Data']);
groups = {'Blank_SAP','SSP_SAP'};
set = 'GCaMP';
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
multiWaitbar('Analyzing sleep model accuracy for GCaMP',0,'Color','P');
for aa = 1:length(groups)
    folderList = dir([groups{1,aa} delim set]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_ModelAccuracy_GCaMP.(groups{1,aa}),animalIDs{1,bb}) == false
            [Results_ModelAccuracy_GCaMP] = AnalyzeModelAccuracy_GCaMP(animalIDs{1,bb},groups{1,aa},set,rootFolder,delim,Results_ModelAccuracy_GCaMP);
        end
        multiWaitbar('Analyzing sleep model accuracy for GCaMP','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end
multiWaitbar('close all');