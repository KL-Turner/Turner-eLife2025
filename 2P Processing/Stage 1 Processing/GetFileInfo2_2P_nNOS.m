function [animalID,hem,fileDate,fileID,imageID,vesselID] = GetFileInfo2_2P(fileName)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% Identify the extension
extInd = strfind(fileName(1,:),'.');
extension = fileName(1,extInd + 1:end);
% Identify the underscores
fileBreaks = strfind(fileName(1,:),'_');
switch extension
    case 'bin'
        animalID = [];
        vesselID = [];
        fileDate = fileName(:,1:fileBreaks(1) - 1);
        fileID = fileName(:,1:fileBreaks(4) - 1);
    case 'mat'
        % Use the known format to parse
        animalID = fileName(:,1:fileBreaks(1) - 1);
        fileDate = fileName(:,fileBreaks(2) + 1:fileBreaks(3) - 1);
        if numel(fileBreaks) > 3
            fileID = fileName(:,fileBreaks(2) + 1:fileBreaks(6) - 1);
            hem = fileName(:,fileBreaks(1) + 1:fileBreaks(2) - 1);
            imageID = fileName(:,fileBreaks(6) + 1: fileBreaks(7) - 1);
            vesselID = fileName(:,fileBreaks(7) + 1:fileBreaks(8) - 1);
        else
            fileDate = [];
            fileID = [];
        end
end

end
