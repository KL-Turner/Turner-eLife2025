function [] = GenerateTrialFigures_IOS_eLife2025(procDataFileIDs,RestingBaselines)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    % load file and gather information
    load(procDataFileID)
    [animalID,~,fileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
    % imaging type
    imagingType = ProcData.notes.imagingType;
    if strcmpi(imagingType,{'Single ROI (SI)'}) == true
        [figHandle] = GenerateTrialFigures_SI_IOS_eLife2025(procDataFileID);
    elseif strcmpi(imagingType,'Single ROI (SSS)') == true
        [figHandle] = GenerateTrialFigures_SSS_IOS_eLife2025(procDataFileID,RestingBaselines);
    elseif any(strcmpi(imagingType,{'Bilateral ROI (SI)','Bilateral ROI (SI,FC)'})) == true
        [figHandle] = GenerateBilateralTrialFigures_IOS_eLife2025(procDataFileID);
    end
    % save the file to directory.
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(figHandle,[dirpath animalID '_' fileID '_SingleTrialFig']);
    close(figHandle)
end