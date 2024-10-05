function [animalID,hem,fileDate,fileID] = GetFileInfo_2P_nNOS(fileName)
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
        hem = [];
        fileDate = fileName(:,1:fileBreaks(1) - 1);
        fileID = fileName(:,1:fileBreaks(4) - 1);
    case 'mat'
        % Use the known format to parse
        animalID = fileName(:,1:fileBreaks(1) - 1);
        hem = fileName(:, fileBreaks(1) + 1:fileBreaks(2) - 1);
        if numel(fileBreaks) > 3
            fileDate = fileName(:,fileBreaks(2) + 1:fileBreaks(3) - 1);
            fileID = fileName(:,fileBreaks(2) + 1:fileBreaks(6) - 1);
        else
            fileDate = [];
            fileID = [];
        end
end

end
