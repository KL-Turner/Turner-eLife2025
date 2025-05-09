function [ROIs] = CheckROIDates_IOS_eLife2025(animalID,ROIs,ROInames,lensMag,imagingType,imagingWavelengths,imagingCamera)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% character list of all WindowCam files
windowCamFilesDir = dir('*_PCO_Cam01.pcoraw');
windowCamDataFiles = {windowCamFilesDir.name}';
windowCamDataFileIDs = char(windowCamDataFiles);
% establish the number of unique days based on file IDs
[~,fileDates,~] = GetFileInfo_IOS_eLife2025(windowCamDataFileIDs);
[uniqueDays,~,DayID] = GetUniqueDays_IOS_eLife2025(fileDates);
firstsFileOfDay = cell(1,length(uniqueDays));
for a = 1:length(uniqueDays)
    FileInd = DayID == a;
    dayFilenames = windowCamDataFileIDs(FileInd,:);
    firstsFileOfDay(a) = {dayFilenames(1,:)};
end
% load existing ROI structure if it exists
ROIFileDir = dir('*_ROIs.mat');
ROIFileName = {ROIFileDir.name}';
ROIFileID = char(ROIFileName);
if exist(ROIFileID,'file')
    load(ROIFileID);
else
    ROIs = [];
end
% create the desired window ROI for each day if it doesn't yet exist
for b = 1:length(firstsFileOfDay)
    fileID = firstsFileOfDay{1,b};
    strDay = ConvertDate_IOS_eLife2025(fileID);
    for c = 1:length(ROInames)
        ROIname = ROInames{1,c};
        if ~isfield(ROIs,strDay) == true || ~isfield(ROIs.(strDay),ROIname) == true
            if any(strcmp(ROInames{1,c},{'LH','RH','fLH','fRH','barrels'})) == true
                % circular ROIs based on lens mag size and imaging type
                [ROIs] = PlaceROIs_IOS_eLife2025(animalID,fileID,ROIs,lensMag,imagingType,imagingWavelengths,imagingCamera);
            elseif strcmp(ROInames{1,c},{'SSS'}) == true
                % ROIS specifically for SSS
                [ROIs] = DrawSagSinusROIs_IOS_eLife2025(animalID,fileID,ROIs,imagingWavelengths,imagingCamera);
            else
                % any additional free hand ROIs
                [ROIs] = CreateFreeHandROIs_IOS_eLife2025(animalID,fileID,ROIs,imagingWavelengths,imagingCamera);
            end
            save([animalID '_ROIs.mat'],'ROIs');
        end
    end
end