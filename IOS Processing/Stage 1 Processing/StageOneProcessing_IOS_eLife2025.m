%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
zap;
% select wavelengths and ROI options
[imagingParameters,animalID,fileList] = SelectImagingParameters_IOS_nNOS();
% process TDMS file and run whisker/pupil tracking
for aa = 1:size(fileList,1)
    [~,~,fileID] = GetFileInfo_IOS_nNOS(fileList(aa,:));
    % determine if file has already been processed
    if ~exist([animalID '_' fileID '_RawData.mat'],'file') == true
        disp(['Analyzing file (' num2str(aa) ' of ' num2str(size(fileList,1)) ')']); disp(' ')
        % import .tdms data
        trialData = ReadTDMSData_IOS_nNOS(animalID,[fileID '.tdms']);
        % start pupil tracker
        [PupilStruct] = TrackPupilDiameter_IOS_nNOS([fileID '_PupilCam.bin'],aa,trialData,imagingParameters);
        % start whisker tracker
        [whiskerAngle] = WhiskerTrackerParallel_IOS_nNOS([fileID '_WhiskerCam.bin'],trialData);
        % create RawData file
        CreateRawDataFile_IOS_nNOS(trialData,fileID,imagingParameters,whiskerAngle,PupilStruct)
    else
        disp(['File (' num2str(aa) '/' num2str(size(fileList,1)) ') already analyzed.']); disp(' ')
    end
end