function [] = MainScript_HRF2020()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
% Purpose: Generates KLT's main and supplemental figs for Turner et al. Manuscript TBD
%
% Scripts used to pre-process the original data are located in the folder "Pre-Processing Scripts".
% Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________

clear; clc; close all;
%% make sure the code repository and data are present in the current directory
currentFolder = pwd;
addpath(genpath(currentFolder));
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootFolder = fullfile(filesep,fileparts{1:end});
    delim = '/';
else
    rootFolder = fullfile(fileparts{1:end});
    delim = '\';
end
% add root folder to Matlab's working directory
addpath(genpath(rootFolder))
%% run the data analysis. The progress bars will show the analysis progress
rerunAnalysis = 'y';
saveFigs = 'n';
if exist('AnalysisResults.mat','file') ~= 2 || strcmp(rerunAnalysis,'y') == true
    multiWaitbar_HRF2020('Analyzing HRF deconvolution',0,'Color','R'); pause(0.5);
    multiWaitbar_HRF2020('Analyzing HRF prediction accuracy',0,'Color','Y'); pause(0.5);
    multiWaitbar_HRF2020('Analyzing full trial HRF predictions',0,'Color','R'); pause(0.5);
    multiWaitbar_HRF2020('Analyzing power spectra',0,'Color','Y'); pause(0.5);
    multiWaitbar_HRF2020('Analyzing coherence',0,'Color','R'); pause(0.5);
    multiWaitbar_HRF2020('Analyzing evoked responses',0,'Color','Y'); pause(0.5);
    multiWaitbar_HRF2020('Analyzing cross correlation',0,'Color','R'); pause(0.5);
    % run analysis and output a structure containing all the analyzed data
    [AnalysisResults] = AnalyzeData_HRF2020(rootFolder,delim);
    multiWaitbar_HRF2020('CloseAll');
else
    disp('Loading analysis results and generating figures...'); disp(' ')
    load('AnalysisResults.mat')
end
%% main figure panels
Fig1_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
Fig2_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
Fig3_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
Fig3_HipTest_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
Fig4_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
Fig5_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
%% fin.
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults] = AnalyzeData_HRF2020(rootFolder,delim)
% animal IDs
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
saveFigs = 'y';
if exist('AnalysisResults.mat','file') == 2
    load('AnalysisResults.mat')
else
    AnalysisResults = [];
end
%% Block [1] Analyze the HRF deconvolution.
runFromStart = 'n';
for aa = 1:length(animalIDs)
    if isfield(AnalysisResults,'HRFs') == false || isfield(AnalysisResults.HRFs,animalIDs{1,aa}) == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeHRF_HRF2020(animalIDs{1,aa},saveFigs,rootFolder,delim,AnalysisResults);
    end
    multiWaitbar_HRF2020('Analyzing HRF deconvolution','Value',aa/length(animalIDs));
end
%% Block [2] Analyze the HRF predictions.
runFromStart = 'n';
for bb = 1:length(animalIDs)
    if isfield(AnalysisResults,'Predictions') == false || isfield(AnalysisResults.Predictions,animalIDs{1,bb}) == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzePredictions_HRF2020(animalIDs{1,bb},rootFolder,AnalysisResults);
    end
    multiWaitbar_HRF2020('Analyzing HRF prediction accuracy','Value',bb/length(animalIDs));
end
%% Block [3] Analyze HRF full trial predictions.
runFromStart = 'n';
for cc = 1:length(animalIDs)
    if isfield(AnalysisResults,'FullTrials') == false || isfield(AnalysisResults.FullTrials,animalIDs{1,cc}) == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeFullTrialPredictions_HRF2020(animalIDs{1,cc},rootFolder,delim,AnalysisResults);
    end
    multiWaitbar_HRF2020('Analyzing full trial HRF predictions','Value',cc/length(animalIDs));
end
%% Block [4] Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
runFromStart = 'n';
for dd = 1:length(animalIDs)
    if isfield(AnalysisResults,'Power') == false || isfield(AnalysisResults.Power,animalIDs{1,dd}) == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzePredictionPower_HRF2020(animalIDs{1,dd},rootFolder,AnalysisResults);
    end
    multiWaitbar_HRF2020('Analyzing power spectra','Value',dd/length(animalIDs));
end
%% Block [5] Analyze coherence.
runFromStart = 'n';
for ee = 1:length(animalIDs)
    if isfield(AnalysisResults,'Coherence') == false || isfield(AnalysisResults.Coherence,animalIDs{1,ee}) == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzePredictionCoherence_HRF2020(animalIDs{1,ee},rootFolder,AnalysisResults);
    end
    multiWaitbar_HRF2020('Analyzing coherence','Value',ee/length(animalIDs));
end
%% Block [6] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
runFromStart = 'n';
for ff = 1:length(animalIDs)
    if isfield(AnalysisResults,'Evoked') == false || isfield(AnalysisResults.Evoked,animalIDs{1,ff}) == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeEvokedResponses_HRF2020(animalIDs{1,ff},rootFolder,AnalysisResults);
    end
    multiWaitbar_HRF2020('Analyzing evoked responses','Value',ff/length(animalIDs));
end
%% Block [7] Analyze the cross-correlation between neural activity and hemodynamics [HbT] (IOS)
runFromStart = 'n';
for gg = 1:length(animalIDs)
    if isfield(AnalysisResults,'XCorr') == false || isfield(AnalysisResults.XCorr,animalIDs{1,gg}) == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeXCorr_HRF2020(animalIDs{1,gg},rootFolder,AnalysisResults);
    end
    multiWaitbar_HRF2020('Analyzing cross correlation','Value',gg/length(animalIDs));
end
%% fin.
disp('Loading analysis results and generating figures...'); disp(' ')

end
