function [animalID,fileDate,fileID] = GetFileInfo_IOS(fileName)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% identify the extension
extInd = strfind(fileName(1,:),'.');
extension = fileName(1,extInd + 1:end);
% identify the underscores
fileBreaks = strfind(fileName(1,:),'_');
switch extension
    case 'bin'
        animalID = [];
        fileDate = fileName(:,1:fileBreaks(1) - 1);
        fileID = fileName(:,1:fileBreaks(4) - 1);
    case 'mat'
        % use the known format to parse
        animalID = fileName(:,1:fileBreaks(1) - 1);
        if numel(fileBreaks) > 3
            fileDate = fileName(:,fileBreaks(1) + 1:fileBreaks(2) - 1);
            fileID = fileName(:,fileBreaks(1) + 1:fileBreaks(5) - 1);
        else
            fileDate = [];
            fileID = [];
        end
    case 'pcoraw'
        animalID = [];
        fileDate = fileName(:,1:fileBreaks(1) - 1);
        fileID = fileName(:,1:fileBreaks(4) - 1);
end