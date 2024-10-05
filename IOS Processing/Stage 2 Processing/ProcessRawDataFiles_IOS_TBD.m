function [procDataFileIDs] = ProcessRawDataFiles_IOS()
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
for a = 1:size(rawDataFileIDs,1)
    rawDataFile = rawDataFileIDs(a,:);
    [animalID,fileDate,fileID] = GetFileInfo_IOS(rawDataFile);
    strDay = ConvertDate_IOS(fileDate);
    procDataFile = ([animalID '_' fileID '_ProcData.mat']);
    if ~exist(procDataFile,'file') == true
        disp(['Creating ProcData file (' num2str(a) '/' num2str(size(rawDataFileIDs,1)) ')']); disp(' ')
        load(rawDataFile);
        % transfer RawData notes to ProcData structure.
        ProcData.notes = RawData.notes;
        % expected durations
        analogExpectedLength = ProcData.notes.trialDuration_sec*ProcData.notes.analogSamplingRate;
        % save solenoid times (in seconds). identify the solenoids by amplitude.
        ProcData.data.stimulations.LPadSol = find(diff(RawData.data.stimulations) == 1)/RawData.notes.analogSamplingRate;
        ProcData.data.stimulations.RPadSol = find(diff(RawData.data.stimulations) == 2)/RawData.notes.analogSamplingRate;
        ProcData.data.stimulations.AudSol = find(diff(RawData.data.stimulations) == 3)/RawData.notes.analogSamplingRate;
        ProcData.data.stimulations.OptoLED = find(diff(RawData.data.stimulations) == 4)/RawData.notes.analogSamplingRate;
        % process neural data into its various forms.
        ProcData.notes.dsFs = 30; % downsampled Fs
        neuralDataTypes = {'cortical_LH','cortical_RH','hippocampus'};
        for c = 1:length(neuralDataTypes)
            neuralDataType = neuralDataTypes{1,c};
            % MUA [300 - 3000]
            [muaPower,~] = ProcessNeuro_IOS(RawData,analogExpectedLength,'MUA',neuralDataType);
            ProcData.data.(neuralDataType).muaPower = muaPower;
            % gamma [30 - 100]
            [gammaBandPower,~] = ProcessNeuro_IOS(RawData,analogExpectedLength,'Gam',neuralDataType);
            ProcData.data.(neuralDataType).gammaBandPower = gammaBandPower;
            % beta [13 - 30 Hz]
            [betaBandPower,~] = ProcessNeuro_IOS(RawData,analogExpectedLength,'Beta',neuralDataType);
            ProcData.data.(neuralDataType).betaBandPower = betaBandPower;
            % alpha [10 - 13 Hz]
            [alphaBandPower,~] = ProcessNeuro_IOS(RawData,analogExpectedLength,'Alpha',neuralDataType);
            ProcData.data.(neuralDataType).alphaBandPower = alphaBandPower;
            % theta [4 - 10 Hz]
            [thetaBandPower,~] = ProcessNeuro_IOS(RawData,analogExpectedLength,'Theta',neuralDataType);
            ProcData.data.(neuralDataType).thetaBandPower = thetaBandPower;
            % delta [1 - 4 Hz]
            [deltaBandPower,~] = ProcessNeuro_IOS(RawData,analogExpectedLength,'Delta',neuralDataType);
            ProcData.data.(neuralDataType).deltaBandPower = deltaBandPower;
        end
        % patch and binarize the whisker angle and set the resting angle to zero degrees.
        [patchedWhisk] = PatchWhiskerAngle_IOS(RawData.data.whiskerAngle,RawData.notes.whiskCamSamplingRate,RawData.notes.trialDuration_sec,RawData.notes.droppedWhiskCamFrameIndex);
        % create filter for whisking/movement
        filtThreshold = 20;
        filtOrder = 2;
        [z,p,k] = butter(filtOrder,filtThreshold/(RawData.notes.whiskCamSamplingRate/2),'low');
        [sos,g] = zp2sos(z,p,k);
        filteredWhiskers = filtfilt(sos,g,patchedWhisk - mean(patchedWhisk));
        resampledWhisk = resample(filteredWhiskers,ProcData.notes.dsFs,RawData.notes.whiskCamSamplingRate);
        % binarize the whisker waveform (wwf)
        threshfile = dir('*_Thresholds.mat');
        if ~isempty(threshfile)
            load(threshfile.name)
        end
        [ok] = CheckForThreshold_IOS(['binarizedWhiskersLower_' strDay],animalID);
        if ok == 0
            [whiskersThresh1,whiskersThresh2] = CreateWhiskThreshold_IOS(resampledWhisk,ProcData.notes.dsFs,strDay);
            Thresholds.(['binarizedWhiskersLower_' strDay]) = whiskersThresh1;
            Thresholds.(['binarizedWhiskersUpper_' strDay]) = whiskersThresh2;
            save([animalID '_Thresholds.mat'],'Thresholds');
        end
        load([animalID '_Thresholds.mat']);
        binWhisk = BinarizeWhiskers_IOS(resampledWhisk,ProcData.notes.dsFs,Thresholds.(['binarizedWhiskersLower_' strDay]),Thresholds.(['binarizedWhiskersUpper_' strDay]));
        [linkedBinarizedWhiskers] = LinkBinaryEvents_IOS(gt(binWhisk,0),[round(ProcData.notes.dsFs/3),0]);
        inds = linkedBinarizedWhiskers == 0;
        restAngle = mean(resampledWhisk(inds));
        ProcData.data.whiskerAngle.angle = resampledWhisk - restAngle;
        ProcData.data.whiskerAngle.binarization = [0,binWhisk,0];
        % downsample and binarize the force sensor.
        trimmedForce = RawData.data.forceSensor(1:min(analogExpectedLength,length(RawData.data.forceSensor)));
        % filter then downsample the Force Sensor waveform to desired frequency
        filtThreshold = 20;
        filtOrder = 2;
        [z,p,k] = butter(filtOrder,filtThreshold/(ProcData.notes.analogSamplingRate/2),'low');
        [sos,g] = zp2sos(z,p,k);
        filtForceSensor = filtfilt(sos,g,trimmedForce);
        ProcData.data.forceSensor.force = resample(filtForceSensor,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
        % binarize the force sensor waveform
        threshfile = dir('*_Thresholds.mat');
        if ~isempty(threshfile)
            load(threshfile.name)
        end
        [ok] = CheckForThreshold_IOS(['binarizedForceSensor_' strDay],animalID);
        if ok == 0
            [forceSensorThreshold] = CreateForceSensorThreshold_IOS(ProcData.data.forceSensor.force,ProcData.notes.dsFs,strDay);
            Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
            save([animalID '_Thresholds.mat'],'Thresholds');
        end
        ProcData.data.forceSensor.binarization = [BinarizeForceSensor_IOS(ProcData.data.forceSensor.force,Thresholds.(['binarizedForceSensor_' strDay])),0];
        % EMG filtering and processing
        fpass = [300,3000];
        trimmedEMG = RawData.data.EMG(1:min(analogExpectedLength,length(RawData.data.EMG)));
        [z,p,k] = butter(3,fpass/(ProcData.notes.analogSamplingRate/2));
        [sos,g] = zp2sos(z,p,k);
        filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
        kernelWidth = 0.5;
        smoothingKernel = gausswin(kernelWidth*ProcData.notes.analogSamplingRate)/sum(gausswin(kernelWidth*ProcData.notes.analogSamplingRate));
        EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
        resampEMG = resample(EMGPwr,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
        ProcData.data.EMG.power = resampEMG;
        % pupil data
        [ProcData] = PatchPupilData_IOS(RawData,ProcData);
        diameter = sqrt(ProcData.data.pupil.area./pi)*2;
        ProcData.data.pupil.diameter = diameter;
        ProcData.notes.pupilCam_mmPerPixel = 0.018;
        ProcData.data.pupil.mmDiameter = diameter.*ProcData.notes.pupilCam_mmPerPixel;
        ProcData.data.pupil.mmArea = ProcData.data.pupil.area.*(ProcData.notes.pupilCam_mmPerPixel^2);
        % save the processed data
        save(procDataFile,'ProcData')
    else
        disp(['ProcData file (' num2str(a) '/' num2str(size(rawDataFileIDs,1)) ') already analyzed.']); disp(' ')
    end
end
% character list of all ProcData files in the directory from ProcessRawDataFiles_IOS.m
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);