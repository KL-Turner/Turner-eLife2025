function [AnalysisResults] = AnalyzeFullTrialPredictions_HRF2020(animalID,rootFolder,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
HRF_hemDataTypes = {'adjLH','adjRH'};
%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % load the Resting baselines structure
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    for aa = 1:length(HRF_hemDataTypes)
        hemDataType = HRF_hemDataTypes{1,aa};
        % character list of all ProcData files
        procDataFileStruct = dir('*_ProcData.mat');
        procDataFiles = {procDataFileStruct.name}';
        procDataFileIDs = char(procDataFiles);
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [AnalysisResults] = GenerateSingleFigures_HRF2020(procDataFileID,hemDataType,RestingBaselines,rootFolder,delim,AnalysisResults);
        end
    end
    AnalysisResults.FullTrials.(animalID) = true;
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
