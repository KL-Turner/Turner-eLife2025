function [imagingParameters,animalID,fileList] = SelectImagingParameters_IOS()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% select imaging type
cameraOptions = {'PCO Edge 5.5','Dalsa Pantera TF 1M60'};
disp('Please select the IOS camera'); disp(' ')
imagingParameters.camera  = SelectImagingType_IOS(cameraOptions);
% select imaging type
imagingOptions = {'Single ROI (SI)','Single ROI (SSS)','Bilateral ROI (SI)','Bilateral ROI (SI,FC)'};
disp('Please select the window type'); disp(' ')
imagingParameters.type = SelectImagingType_IOS(imagingOptions);
% select imaging type
wavelengthOptions = {'Green','Lime','Blue','Green & Blue','Lime & Blue','Red, Green, & Blue','Red, Lime, & Blue'};
disp('Please select the imaging wavelengths'); disp(' ')
imagingParameters.wavelengths = SelectImagingType_IOS(wavelengthOptions);
% select imaging type
pupilOptions = {'Yes','No'};
disp('Please select whether to run the pupil tracker'); disp(' ')
imagingParameters.pupilOptions = SelectImagingType_IOS(pupilOptions);
% animal ID
animalID = input('Enter animal ID: ','s'); disp(' ')
% file IDs
whiskerFileDir = dir('*_WhiskerCam.bin');
whiskerFileNames = {whiskerFileDir.name}';
fileList = char(whiskerFileNames);