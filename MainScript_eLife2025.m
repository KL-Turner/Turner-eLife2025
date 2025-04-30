function [] = MainScript_eLife2025()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
currentFolder = pwd;
addpath(genpath(currentFolder));
zap;
saveState = true;
runAnalysis = false;
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootFolder = fullfile(filesep,fileparts{1:end});
    delim = '/';
else
    rootFolder = fullfile(fileparts{1:end});
    delim = '\';
end
addpath(genpath(rootFolder))
if runAnalysis == true
    RunAnalysis_eLife2025(rootFolder,delim)
end
% analysis subfunctions
disp('Loading analysis results and generating figures...'); disp(' ')

% main figure panels
Fig1_eLife2025(rootFolder,saveState,delim)
Fig2_eLife2025(rootFolder,saveState,delim)
Fig3_eLife2025(rootFolder,saveState,delim)
Fig4_eLife2025(rootFolder,saveState,delim)
Fig5_eLife2025(rootFolder,saveState,delim)
Fig6_eLife2025(rootFolder,saveState,delim)
Fig7_eLife2025(rootFolder,saveState,delim)

% supplemental figure panels
FigS1_eLife2025(rootFolder,saveState,delim)
FigS2_eLife2025(rootFolder,saveState,delim)
FigS3_eLife2025(rootFolder,saveState,delim)
FigS4_eLife2025(rootFolder,saveState,delim)
FigS5_eLife2025(rootFolder,saveState,delim)
FigS6_eLife2025(rootFolder,saveState,delim)
FigS7_eLife2025(rootFolder,saveState,delim)

% go back to original directory
cd(currentFolder)

function [] = RunAnalysis_eLife2025(rootFolder,delim)
% evoked responses
AnalyzeEvokedResponses_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzeEvokedResponses_GCaMP_Handler_eLife2025(rootFolder,delim,false)
AnalyzeEvokedResponses_2P_Handler_eLife2025(rootFolder,delim,false)
AnalyzeEvokedResponses_Pulse_Handler_eLife2025(rootFolder,delim,false)

% HbT, HbO, HbR, GCaMP
AnalyzeIntrinsicSignals_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzeIntrinsicSignals_GCaMP_Handler_eLife2025(rootFolder,delim,false)
AnalyzeIntrinsicSignals_Pulse_Handler_eLife2025(rootFolder,delim,false)
% AnalyzeHemoGFPRelationship_EGFP_Handler_eLife2025(rootFolder,delim,false)
% AnalyzeArousalDerivative_Ephys_Handler_eLife2025(rootFolder,delim,false)

% arteriole diameter and baseline shift
AnalyzeArterioleBaseline_2P_Handler_eLife2025(rootFolder,delim,false)
AnalyzeArterioleBaselineShift_2P_Handler_eLife2025(rootFolder,delim,true)
AnalyzeArterioleDiameter_2P_Handler_eLife2025(rootFolder,delim,false)

% power spectral density
AnalyzePowerSpectrum_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzePowerSpectrum_GCaMP_Handler_eLife2025(rootFolder,delim,false)
AnalyzePowerSpectrum_LFP_Handler_eLife2025(rootFolder,delim,false)
AnalyzePowerSpectrum_2P_Handler_eLife2025(rootFolder,delim,false)

% pre-whitened power spectral density
% AnalyzePreWhitenedPowerSpectrum_Ephys_Handler_eLife2025(rootFolder,delim,false)
% AnalyzePreWhitenedPowerSpectrum_GCaMP_Handler_eLife2025(rootFolder,delim,false)
% AnalyzePreWhitenedPowerSpectrum_2P_Handler_eLife2025(rootFolder,delim,false)

% coherence between hemipheres
AnalyzeBilateralCoherence_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzeBilateralCoherence_GCaMP_Handler_eLife2025(rootFolder,delim,false)

% coherence between neural-hemo
AnalyzeNeuralHemoCoherence_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzeNeuralHemoCoherence_GCaMP_Handler_eLife2025(rootFolder,delim,false)

% Pearson's correlation between hemispheres
AnalyzePearsonCorrelation_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzePearsonCorrelation_GCaMP_Handler_eLife2025(rootFolder,delim,false)

% cross correlation between neural-hemo
AnalyzeCrossCorrelation_Ephys_Handler_eLife2025(rootFolder,delim,false)
% AnalyzeHbTCrossCorrelation_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzeCrossCorrelation_GCaMP_Handler_eLife2025(rootFolder,delim,false)
% AnalyzeCrossCorrelation_EGFP_Handler_eLife2025(rootFolder,delim,false)

% pupil analysis
AnalyzePupilArea_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzePupilEvokedResponses_Ephys_Handler_eLife2025(rootFolder,delim,false)
% AnalyzePupilCoherence_Ephys_Handler_eLife2025(rootFolder,delim,false)
% AnalyzePupilCrossCorrelation_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzePupilInterBlinkInterval_Ephys_Handler_eLife2025(rootFolder,delim,false)

% state probability
AnalyzeArousalStateProbability_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzeArousalStateProbability_GCaMP_Handler_eLife2025(rootFolder,delim,false)

% sleep model accuracy
AnalyzeModelAccuracy_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzeModelAccuracy_GCaMP_Handler_eLife2025(rootFolder,delim,false)

% quanitification of whisking
AnalyzeWhiskingBehavior_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzeWhiskingBehavior_GCaMP_Handler_eLife2025(rootFolder,delim,false)

% arousal state transitions
AnalyzeArousalTransitions_Ephys_Handler_eLife2025(rootFolder,delim,false)
AnalyzeArousalTransitions_GCaMP_Handler_eLife2025(rootFolder,delim,false)

% hemodynamic response function
AnalyzeHRF_Ephys_Handler_eLife2025(rootFolder,delim,false)