function [] = AnalyzeNeuralHemoCoherence_GCaMP_Handler(rootFolder,delim,runFromStart)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
if runFromStart == true
    Results_NeuralHemoCoher_GCaMP.Blank_SAP = [];
    Results_NeuralHemoCoher_GCaMP.SSP_SAP = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_NeuralHemoCoher_GCaMP.mat','file') == 2
        load('Results_NeuralHemoCoher_GCaMP.mat','-mat')
    else
        Results_NeuralHemoCoher_GCaMP.Blank_SAP = [];
        Results_NeuralHemoCoher_GCaMP.SSP_SAP = [];
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
multiWaitbar('Analyzing neural-hemo coherence for GCaMP',0,'Color','P');
for aa = 1:length(groups)
    folderList = dir([groups{1,aa} delim set]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_NeuralHemoCoher_GCaMP.(groups{1,aa}),(animalIDs{1,bb})) == false
            [Results_NeuralHemoCoher_GCaMP] = AnalyzeNeuralHemoCoherence_GCaMP(animalIDs{1,bb},groups{1,aa},set,rootFolder,delim,Results_NeuralHemoCoher_GCaMP);
        end
        multiWaitbar('Analyzing neural-hemo coherence for GCaMP','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end
multiWaitbar('close all');