%% this code was written to organize all the files in different folders for different mouse and analyze them and list in a group folder. 
function BatchRun_code_OFT()
%{
clc; clear all; close all;

cd('I:\ResizedVideos\SSP_SAP_Blank')
dirList = dir(); 
rmfolders = isfolder({dirList.name});
dirList(~rmfolders) = [];
excludenames = ["Control_c57","Experiement_group"];
dirList(~ismember({dirList.name},excludenames(:))) = [];
for exptypeN = 1:1:length(dirList)
    experimentPath = fullfile(dirList(exptypeN).folder,dirList(exptypeN).name);
    sessionList = dir(experimentPath);
    excludesession = [".",".."];
    sessionList(ismember({sessionList.name},excludesession(:))) = [];
    for sessiontypeN = 1:1:length(sessionList)
        sessionPath = fullfile(sessionList(sessiontypeN).folder,sessionList(sessiontypeN).name);
        subjectList = dir(sessionPath);
        excludeFolders = [".",".."];
        subjectList(ismember({subjectList.name},excludeFolders(:))) = [];
        for subjectN = 1:1:length(subjectList)
            subjectPath = fullfile(subjectList(subjectN).folder,subjectList(subjectN).name);
            KalmanFilterCorrection_OFT(subjectPath);
        end
    end
end
fprintf('Done KalmanFilterCorrection \n')
%}

% clear all; clc
% cd('I:\ResizedVideos\SSP_SAP_Blank\Recalculated_Raw_Population')
% subjectList = dir(); 
% rmfolders = isfolder({subjectList.name});
% excludeFolders = [".",".."];
% subjectList(~rmfolders) = [];
% subjectList(ismember({subjectList.name},excludeFolders(:))) = [];
% for subjectN = 1:1:length(subjectList)
%      subjectPath = fullfile(subjectList(subjectN).folder,subjectList(subjectN).name);
%      CorrectedOFTAnalysis(subjectPath)
% end


% clear all; clc
% cd('I:\ResizedVideos\SSP_SAP_Blank\Recalculated_Raw_Population')
% subjectList = dir(); 
% rmfolders = isfolder({subjectList.name});
% excludeFolders = [".",".."];
% subjectList(~rmfolders) = [];
% subjectList(ismember({subjectList.name},excludeFolders(:))) = [];
% for subjectN = 1:1:length(subjectList)
%      subjectPath = fullfile(subjectList(subjectN).folder,subjectList(subjectN).name);
%      copyDatatocompare(subjectPath)
% end


fprintf('Done CorrectedOFTAnalysis \n')
% FinalOFTAnalysis()
fprintf('Done FinalOFTAnalysis \n')

% don't need to run this code now. 