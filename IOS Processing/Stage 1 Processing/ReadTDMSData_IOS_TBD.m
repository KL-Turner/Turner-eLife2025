function [TDMSFile] = ReadTDMSData_IOS(animalID,fileName)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% convert the .tdms file into something that Matlab understands
[tempStruct,~] = ConvertTDMS_IOS(0,fileName);
% extract whisker camera info and transfer from tempStruct
TDMSFile.experimenter = tempStruct.Data.Root.Experimenter;
TDMSFile.animal_ID = tempStruct.Data.Root.Animal_ID;
TDMSFile.animalID = animalID;
TDMSFile.hemisphere = tempStruct.Data.Root.Hemisphere;
TDMSFile.solenoidPSI = tempStruct.Data.Root.Solenoid_PSI;
TDMSFile.isofluraneTime = tempStruct.Data.Root.Isoflurane_time;
TDMSFile.sessionID = tempStruct.Data.Root.Session_ID;
TDMSFile.amplifierGain = tempStruct.Data.Root.Amplifier_Gain;
TDMSFile.Opto_LED_mW = tempStruct.Data.Root.Opto_LED_mW;
TDMSFile.wavelengths = tempStruct.Data.Root.Wavelengths;
TDMSFile.lensMag = tempStruct.Data.Root.Lens_mag;
TDMSFile.CBVCamSamplingRate = tempStruct.Data.Root.CBV_Cam_Fs;
TDMSFile.CBVCameraID = tempStruct.Data.Root.CBV_CameraID;
TDMSFile.CBVCamTriggerMode = tempStruct.Data.Root.PCO_TriggerMode;
TDMSFile.CBVCamExposureTime_microsec = num2str(str2double(tempStruct.Data.Root.PCO_ExposureTime_Sec)*1000*1000);
TDMSFile.CBVCamTimeStampMode = tempStruct.Data.Root.PCO_TimeStampMode;
TDMSFile.CBVCamBitDepth = '16';
TDMSFile.CBVCamPixelWidthx0 = tempStruct.Data.Root.PCO_ROI_x0;
TDMSFile.CBVCamPixelHeighty0 = tempStruct.Data.Root.PCO_ROI_y0;
TDMSFile.CBVCamPixelWidth = tempStruct.Data.Root.PCO_ROI_x1;
TDMSFile.CBVCamPixelHeight = tempStruct.Data.Root.PCO_ROI_y1;
TDMSFile.CBVCamBinning = tempStruct.Data.Root.PCO_Binning_vert;
TDMSFile.CBVCamBinningHorz = tempStruct.Data.Root.PCO_Binning_horz;
TDMSFile.CBVCamBinningVert = tempStruct.Data.Root.PCO_Binning_vert;
TDMSFile.whiskCamSamplingRate = tempStruct.Data.Root.Whisk_Cam_Fs;
TDMSFile.webCamSamplingRate = tempStruct.Data.Root.Web_Cam_Fs;
TDMSFile.pupilCamSamplingRate = tempStruct.Data.Root.Pupil_Cam_Fs;
TDMSFile.analogSamplingRate = tempStruct.Data.Root.Analog_Fs;
TDMSFile.trialDuration_sec = tempStruct.Data.Root.TrialDuration_sec;
TDMSFile.pupilCamPixelWidth = tempStruct.Data.Root.PupilCam_Width_pix;
TDMSFile.pupilCamPixelHeight = tempStruct.Data.Root.PupilCam_Height_pix;
TDMSFile.whiskCamPixelWidth = tempStruct.Data.Root.WhiskCam_Width_pix;
TDMSFile.whiskCamPixelHeight = tempStruct.Data.Root.WhiskCam_Height_pix;
TDMSFile.droppedPupilCamFrameIndex = tempStruct.Data.Root.PupilCam_DroppedFrameIndex;
TDMSFile.droppedWhiskCamFrameIndex = tempStruct.Data.Root.WhiskCam_DroppedFrameIndex;
TDMSFile.Sol_DutyCycle = tempStruct.Data.Root.Sol_DutyCycle;
TDMSFile.Sol_Freq = tempStruct.Data.Root.Sol_Freq;
TDMSFile.Sol_Duration_sec = tempStruct.Data.Root.Sol_Duration_sec;
TDMSFile.LED_DutyCycle = tempStruct.Data.Root.LED_DutyCycle;
TDMSFile.LED_Freq = tempStruct.Data.Root.LED_Freq;
TDMSFile.LED_Duration_sec = tempStruct.Data.Root.LED_Duration_sec;
TDMSFile.Interstim_sec = tempStruct.Data.Root.Interstim_sec;
TDMSFile.Stim_Offset_sec = tempStruct.Data.Root.Stim_Offset_sec;
% pre-allocate - data is contained in .vals folder in rows with corresponding labels in .names
TDMSFile.data.vals = NaN*ones(length(tempStruct.Data.MeasuredData),length(tempStruct.Data.MeasuredData(1).Data));
TDMSFile.data.names = cell(length(tempStruct.Data.MeasuredData),1);
% pull data from tempStruct and allocate it in the proper areas
for k = 1:length(tempStruct.Data.MeasuredData)
    TDMSFile.data.vals(k,:) = tempStruct.Data.MeasuredData(k).Data;
    TDMSFile.data.names{k} = strrep(tempStruct.Data.MeasuredData(k).Name,'Analog_Data','');
end