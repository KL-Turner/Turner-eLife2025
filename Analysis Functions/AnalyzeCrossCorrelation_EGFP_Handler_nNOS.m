function [] = AnalyzeCrossCorrelation_EGFP_Handler(rootFolder,delim,runFromStart)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
if runFromStart == true
    Results_CrossCorr_EGFP.EGFP = [];
elseif runFromStart == false
    cd([rootFolder delim 'Results_Turner']);
    % load existing results structure, if it exists
    if exist('Results_CrossCorr_EGFP.mat','file') == 2
        load('Results_CrossCorr_EGFP.mat','-mat')
    else
        Results_CrossCorr_EGFP.EGFP = [];
    end
end
cd([rootFolder delim 'Data']);
groups = {'EGFP'};
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
multiWaitbar('Analyzing cross correlation for EGFP',0,'Color','P');
for aa = 1:length(groups)
    folderList = dir([groups{1,aa} delim set]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    for bb = 1:length(animalIDs)
        if isfield(Results_CrossCorr_EGFP.(groups{1,aa}),(animalIDs{1,bb})) == false
            [Results_CrossCorr_EGFP] = AnalyzeCrossCorrelation_EGFP(animalIDs{1,bb},groups{1,aa},set,rootFolder,delim,Results_CrossCorr_EGFP);
        end
        multiWaitbar('Analyzing cross correlation for EGFP','Value',cc/waitBarLength); pause(0.5);
        cc = cc + 1;
    end
end
multiWaitbar('close all');