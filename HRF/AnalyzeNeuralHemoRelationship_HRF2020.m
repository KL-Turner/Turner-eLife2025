function [AnalysisResults] = AnalyzeNeuralHemoRelationship_HRF2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the HRF deconvolution of neural activity and HbT
%________________________________________________________________________________________________________________________

%% only run analysis for valid animal IDs
dataLocation = [rootFolder '\' animalID '\Bilateral Imaging\'];
cd(dataLocation)
% find and load manual baseline event information
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% find and load Manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID)
% find and load EventData.mat structure
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)
%% left hemisphere gamma-band and HbT after a contralateral stimulus
% extract neural and hemodynamic data from event structure
[LH_NeuralDataStruct,LH_NeuralFiltArray] = SelectConvolutionBehavioralEvents_HRF2020(EventData.cortical_LH.gammaBandPower,'Contra','adjLH');
[LH_HemoDataStruct,LH_HemoFiltArray] = SelectConvolutionBehavioralEvents_HRF2020(EventData.CBV_HbT.adjLH,'Contra','adjLH');
% remove events that don't meet criteria
[LH_NeuralData,~,~,~] = RemoveInvalidData_HRF2020(LH_NeuralDataStruct.NormData(LH_NeuralFiltArray,:),LH_NeuralDataStruct.fileIDs(LH_NeuralFiltArray,:),LH_NeuralDataStruct.duration(LH_NeuralFiltArray,:),LH_NeuralDataStruct.eventTime(LH_NeuralFiltArray,:),ManualDecisions);
[LH_HemoData,~,~,~] = RemoveInvalidData_HRF2020(LH_HemoDataStruct.data(LH_HemoFiltArray,:),LH_HemoDataStruct.fileIDs(LH_HemoFiltArray,:),LH_HemoDataStruct.duration(LH_HemoFiltArray,:),LH_HemoDataStruct.eventTime(LH_HemoFiltArray,:),ManualDecisions);
%% right hemisphere gamma-band and HbT after a contralateral stimulus
% extract neural and hemodynamic data from event structure
[RH_NeuralDataStruct,RH_NeuralFiltArray] = SelectConvolutionBehavioralEvents_HRF2020(EventData.cortical_RH.gammaBandPower,'Contra','adjRH');
[RH_HemoDataStruct,RH_HemoFiltArray] = SelectConvolutionBehavioralEvents_HRF2020(EventData.CBV_HbT.adjRH,'Contra','adjRH');
% remove events that don't meet criteria
[RH_NeuralData,~,~,~] = RemoveInvalidData_HRF2020(RH_NeuralDataStruct.NormData(RH_NeuralFiltArray,:),RH_NeuralDataStruct.fileIDs(RH_NeuralFiltArray,:),RH_NeuralDataStruct.duration(RH_NeuralFiltArray,:),RH_NeuralDataStruct.eventTime(RH_NeuralFiltArray,:),ManualDecisions);
[RH_HemoData,~,~,~] = RemoveInvalidData_HRF2020(RH_HemoDataStruct.data(RH_HemoFiltArray,:),RH_HemoDataStruct.fileIDs(RH_HemoFiltArray,:),RH_HemoDataStruct.duration(RH_HemoFiltArray,:),RH_HemoDataStruct.eventTime(RH_HemoFiltArray,:),ManualDecisions);
%% take the max of each
% pre-allocate
LH_gamma = nan(size(LH_NeuralData,1),1);
LH_HbT = nan(size(LH_HemoData,1),1);
RH_gamma = nan(size(RH_NeuralData,1),1);
RH_HbT = nan(size(RH_HemoData,1),1);
offset = 2; 
samplingRate = 30;
% left hem
for aa = 1:size(LH_NeuralData,1)
    LH_gammaArray = LH_NeuralData(aa,:);
    LH_gamma(aa,1) = max(LH_gammaArray(offset*samplingRate:(offset + 3)*samplingRate));
    LH_hemoArray = LH_HemoData(aa,:);
    LH_HbT(aa,1) = max(LH_hemoArray(offset*samplingRate:(offset + 3)*samplingRate));
end
% right hem
for bb = 1:size(RH_NeuralData,1)
    RH_gammaArray = RH_NeuralData(bb,:);
    RH_gamma(bb,1) = max(RH_gammaArray(offset*samplingRate:(offset + 3)*samplingRate));
    RH_hemoArray = RH_HemoData(bb,:);
    RH_HbT(bb,1) = max(RH_hemoArray(offset*samplingRate:(offset + 3)*samplingRate));
end
% save results
AnalysisResults.GammaHbT.(animalID).LH_gamma = LH_gamma;
AnalysisResults.GammaHbT.(animalID).LH_HbT = LH_HbT;
AnalysisResults.GammaHbT.(animalID).RH_gamma = RH_gamma;
AnalysisResults.GammaHbT.(animalID).RH_HbT = RH_HbT;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end
