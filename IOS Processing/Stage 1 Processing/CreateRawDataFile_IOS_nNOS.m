function [] = CreateRawDataFile_IOS_nNOS(trialData,fileID,imagingParameters,whiskerAngle,PupilStruct)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% left, right, and hippocampal electrodes
dataRow = strcmp(trialData.data.names,'Cortical_LH');
cortical_LH = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
dataRow = strcmp(trialData.data.names,'Cortical_RH');
cortical_RH = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
dataRow = strcmp(trialData.data.names,'Hippocampus');
hippocampus = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
% left, right, auditory solenoids. combine the arrays together.
dataRow = strcmp(trialData.data.names,'LPadSol');
LPadSol = gt(trialData.data.vals(dataRow,:),0.5)*1; % ID amplitude is 1
dataRow = strcmp(trialData.data.names,'RPadSol');
RPadSol = gt(trialData.data.vals(dataRow,:),0.5)*2; % ID amplitude is 2
dataRow = strcmp(trialData.data.names,'AudSol');
AudSol = gt(trialData.data.vals(dataRow,:),0.5)*3;  % ID amplitude is 3
dataRow = strcmp(trialData.data.names,'OptoLED');
OptoLED = gt(trialData.data.vals(dataRow,:),0.5)*4; % ID amplitude is 4
stimulations = LPadSol + RPadSol + AudSol + OptoLED;
% force sensor and EMG
dataRow = strcmp(trialData.data.names,'Force_Sensor');
forceSensor = trialData.data.vals(dataRow,:);
dataRow = strcmp(trialData.data.names,'EMG');
EMG = trialData.data.vals(dataRow,:)/str2double(trialData.amplifierGain);
% save the notes
RawData.notes.iosCamera = imagingParameters.camera;
RawData.notes.imagingType = imagingParameters.type;
RawData.notes.imagingWavelengths = imagingParameters.wavelengths;
RawData.notes.experimenter = trialData.experimenter;
RawData.notes.animalID = trialData.animalID;
RawData.notes.hemisphere = trialData.hemisphere;
RawData.notes.solenoidPSI = str2double(trialData.solenoidPSI);
RawData.notes.isofluraneTime = str2double(trialData.isofluraneTime);
RawData.notes.sessionID = trialData.sessionID;
RawData.notes.amplifierGain = str2double(trialData.amplifierGain);
RawData.notes.Opto_LED_mW = trialData.Opto_LED_mW;
RawData.notes.wavelengths = trialData.wavelengths;
RawData.notes.lensMag = trialData.lensMag;
RawData.notes.CBVCamSamplingRate = str2double(trialData.CBVCamSamplingRate);
RawData.notes.whiskCamSamplingRate = str2double(trialData.whiskCamSamplingRate);
RawData.notes.webCamSamplingRate = str2double(trialData.webCamSamplingRate);
RawData.notes.pupilCamSamplingRate = str2double(trialData.pupilCamSamplingRate);
RawData.notes.analogSamplingRate = str2double(trialData.analogSamplingRate);
RawData.notes.trialDuration_sec = str2double(trialData.trialDuration_sec);
RawData.notes.CBVCameraID = trialData.CBVCameraID;
RawData.notes.CBVCamTriggerMode = trialData.CBVCamTriggerMode;
RawData.notes.CBVCamExposureTime_microsec = str2double(trialData.CBVCamExposureTime_microsec);
RawData.notes.CBVCamTimeStampMode = trialData.CBVCamTimeStampMode;
RawData.notes.CBVCamBitDepth = str2double(trialData.CBVCamBitDepth);
RawData.notes.CBVCamPixelWidthx0 = trialData.CBVCamPixelWidthx0;
RawData.notes.CBVCamPixelHeighty0 = trialData.CBVCamPixelHeighty0;
RawData.notes.CBVCamPixelWidth = str2double(trialData.CBVCamPixelWidth);
RawData.notes.CBVCamPixelHeight = str2double(trialData.CBVCamPixelHeight);
RawData.notes.CBVCamBinning = trialData.CBVCamBinning;
RawData.notes.CBVCamBinningHorz = trialData.CBVCamBinningHorz;
RawData.notes.CBVCamBinningVert = trialData.CBVCamBinningVert;
RawData.notes.pupilCamPixelWidth = str2double(trialData.pupilCamPixelWidth);
RawData.notes.pupilCamPixelHeight = str2double(trialData.pupilCamPixelHeight);
RawData.notes.whiskCamPixelHeight = str2double(trialData.whiskCamPixelHeight);
RawData.notes.whiskCamPixelWidth = str2double(trialData.whiskCamPixelWidth);
RawData.notes.droppedPupilCamFrameIndex = trialData.droppedPupilCamFrameIndex;
RawData.notes.droppedWhiskCamFrameIndex = trialData.droppedWhiskCamFrameIndex;
RawData.notes.solenoidDutyCycle = str2double(trialData.Sol_DutyCycle);
RawData.notes.solenoidFreq = str2double(trialData.Sol_Freq);
RawData.notes.solenoidDuration_sec = str2double(trialData.Sol_Duration_sec);
RawData.notes.LEDdutyCycle = str2double(trialData.LED_DutyCycle);
RawData.notes.LEDfreq = str2double(trialData.LED_Freq);
RawData.notes.LEDduration_sec = str2double(trialData.LED_Duration_sec);
RawData.notes.interstim_sec = str2double(trialData.Interstim_sec);
RawData.notes.stimOffset_sec = str2double(trialData.Stim_Offset_sec);
% save the data
RawData.data.cortical_LH = cortical_LH;
RawData.data.cortical_RH = cortical_RH;
RawData.data.hippocampus = hippocampus;
RawData.data.forceSensor = forceSensor;
RawData.data.EMG = EMG;
RawData.data.whiskerAngle = whiskerAngle;
RawData.data.pupil = PupilStruct;
RawData.data.stimulations = stimulations;
disp('File Created. Saving RawData File...'); disp(' ')
save([trialData.animalID '_' fileID '_RawData'],'RawData')