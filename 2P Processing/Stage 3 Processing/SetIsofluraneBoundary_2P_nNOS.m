function [] = SetIsofluraneBoundary_2P(RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________

mergedDataFileIDs = uigetfile('*_MergedData.mat','Multiselect','on');
for aa = 1:size(mergedDataFileIDs,1)
    mergedDataFileID = mergedDataFileIDs{1,aa};
    [figHandle] = GenerateSingleFigures_2P(mergedDataFileID,'manualSelection','n',RestingBaselines);
    MergedData.notes.isoStart = input('Input isoflurane start time: '); disp(' ')
    MergedData.notes.isoEnd = input('Input isoflurane end time: '); disp(' ')
    save(mergedDataFileID,'MergedData')
    close(figHandle)
end