function [] = Process2PDataFiles_2P_nNOS(labviewDataFiles,mscanDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the force sensor, whisker angle, EMG, and neural bands. Create a threshold for binarized movement and
%            whisking if one does not already exist.
%________________________________________________________________________________________________________________________

%% MScan data file analysis
for a = 1:size(mscanDataFiles,1)
    mscanDataFile = mscanDataFiles(a,:);
    load(mscanDataFile);
    % Skip the file if it has already been processed
    if MScanData.notes.checklist.processData == false
        disp(['Analyzing MScan neural bands and analog signals for file number ' num2str(a) ' of ' num2str(size(mscanDataFiles,1)) '...']); disp(' ');
        animalID = MScanData.notes.animalID;
        imageID = MScanData.notes.imageID;
        vesselID = MScanData.notes.vesselID;
        date = MScanData.notes.date;
        strDay = ConvertDate_2P_nNOS(date);
        downSampledFs = 30;
        MScanData.notes.dsFs = downSampledFs;
        expectedLength = floor((MScanData.notes.numberOfFrames/MScanData.notes.frameRate)*MScanData.notes.analogSamplingRate);
        %% Process neural data into its various forms.
        % MUA Band [300 - 3000]
        neuralTypes = {'corticalNeural','hippocampalNeural'};
        neuralFields = {'cortical','hippocampal'};
        for q = 1:length(neuralTypes)
            neuralType = neuralTypes{1,q};
            neuralField = neuralFields{1,q};
            % MUA [300-3000 Hz]
            [MScanData.data.(neuralField).muaPower] = ProcessNeuro_2P_nNOS(MScanData,expectedLength,'MUA',neuralType);
            % Gamma [30-100 Hz]
            [MScanData.data.(neuralField).gammaBandPower] = ProcessNeuro_2P_nNOS(MScanData,expectedLength,'Gam',neuralType);
            % Beta [13-30 Hz]
            [MScanData.data.(neuralField).betaBandPower] = ProcessNeuro_2P_nNOS(MScanData,expectedLength,'Beta',neuralType);
            % Alpha [10-13 Hz]
            [MScanData.data.(neuralField).alphaBandPower] = ProcessNeuro_2P_nNOS(MScanData,expectedLength,'Alpha',neuralType);
            % Theta [4-10 Hz]
            [MScanData.data.(neuralField).thetaBandPower] = ProcessNeuro_2P_nNOS(MScanData,expectedLength,'Theta',neuralType);
            % Delta [1-4 Hz]
            [MScanData.data.(neuralField).deltaBandPower] = ProcessNeuro_2P_nNOS(MScanData,expectedLength,'Delta',neuralType);
        end
        
        %% Downsample and binarize the force sensor.
        trimmedForceM = MScanData.data.forceSensor(1:min(expectedLength,length(MScanData.data.forceSensor)));
        % Filter then downsample the Force Sensor waveform to desired frequency
        filtThreshold = 20;
        filtOrder = 2;
        [z,p,k] = butter(filtOrder,filtThreshold/(MScanData.notes.analogSamplingRate/2),'low');
        [sos,g] = zp2sos(z,p,k);
        filtForceSensorM = filtfilt(sos,g,trimmedForceM);
        MScanData.data.dsForceSensorM = resample(filtForceSensorM,downSampledFs,MScanData.notes.analogSamplingRate);
        % Binarize the force sensor waveform
        threshfile = dir('*_Thresholds.mat');
        if ~isempty(threshfile)
            load(threshfile.name)
        end
        [ok] = CheckForThreshold_2P_nNOS(['binarizedForceSensor_' strDay],animalID);
        if ok == 0
            [forceSensorThreshold] = CreateForceSensorThreshold_2P_nNOS(MScanData.data.dsForceSensorM);
            Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
            save([animalID '_Thresholds.mat'],'Thresholds');
        end
        MScanData.data.binForceSensorM = BinarizeForceSensor_2P_nNOS(MScanData.data.dsForceSensorM,Thresholds.(['binarizedForceSensor_' strDay]));
        
        %% EMG
        fpass = [300,3000];
        trimmedEMG = MScanData.data.EMG(1:min(expectedLength,length(MScanData.data.EMG)));
        [z,p,k] = butter(3,fpass/(MScanData.notes.analogSamplingRate/2));
        [sos,g] = zp2sos(z,p,k);
        filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
        kernelWidth = 0.5;
        smoothingKernel = gausswin(kernelWidth*MScanData.notes.analogSamplingRate)/sum(gausswin(kernelWidth*MScanData.notes.analogSamplingRate));
        EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
        resampEMG = resample(EMGPwr,MScanData.notes.dsFs,MScanData.notes.analogSamplingRate);
        MScanData.data.filtEMG = resampEMG;  
        % nancheck
        if sum(isnan(MScanData.data.filtEMG)) > 0
            keyboard
        end
        % inf check
        if sum(isinf(MScanData.data.filtEMG)) > 0
            keyboard
        end

        %% Save the data, set checklist to true
        MScanData.notes.checklist.processData = true;
        save([animalID '_' date '_'  imageID '_' vesselID '_MScanData'],'MScanData')
    else
        disp([mscanDataFile ' has already been processed. Continuing...']); disp(' ');
    end
end

%% LabVIEW data file analysis
for b = 1:size(labviewDataFiles,1)
    labviewDataFile = labviewDataFiles(b,:);
    load(labviewDataFile);
    if LabVIEWData.notes.checklist.processData == false
        disp(['Analyzing LabVIEW analog signals and whisker angle for file number ' num2str(b) ' of ' num2str(size(labviewDataFiles, 1)) '...']); disp(' ');
        [animalID,hem,fileDate,fileID] = GetFileInfo_2P_nNOS(labviewDataFile);
        strDay = ConvertDate_2P_nNOS(fileDate);
        expectedLength = LabVIEWData.notes.trialDuration_Seconds*LabVIEWData.notes.analogSamplingRate_Hz;

        %% Patch and binarize the whisker angle and set the resting angle to zero degrees.
        [patchedWhisk] = PatchWhiskerAngle_2P_nNOS(LabVIEWData.data.whiskerAngle,LabVIEWData.notes.whiskerCamSamplingRate_Hz,LabVIEWData.notes.trialDuration_Seconds,LabVIEWData.notes.droppedWhiskerCamFrameIndex);
        % Create filter for whisking/movement
        downSampledFs = 30;
        filtThreshold = 20;
        filtOrder = 2;
        [z,p,k] = butter(filtOrder,filtThreshold/(LabVIEWData.notes.whiskerCamSamplingRate_Hz/2),'low');
        [sos,g] = zp2sos(z,p,k);
        filteredWhiskers = filtfilt(sos, g, patchedWhisk - mean(patchedWhisk));
        resampledWhisk = resample(filteredWhiskers, downSampledFs,LabVIEWData.notes.whiskerCamSamplingRate_Hz);
        % Binarize the whisker waveform (wwf)
        threshfile = dir('*_Thresholds.mat');
        if ~isempty(threshfile)
            load(threshfile.name)
        end 
        [ok] = CheckForThreshold_2P_nNOS(['binarizedWhiskersLower_' strDay],animalID);
        if ok == 0
            [whiskersThresh1, whiskersThresh2] = CreateWhiskThreshold_2P_nNOS(resampledWhisk,downSampledFs);
            Thresholds.(['binarizedWhiskersLower_' strDay]) = whiskersThresh1;
            Thresholds.(['binarizedWhiskersUpper_' strDay]) = whiskersThresh2;
            save([animalID '_Thresholds.mat'], 'Thresholds');
        end
        load([animalID '_Thresholds.mat']);
        binWhisk = BinarizeWhiskers_2P_nNOS(resampledWhisk,downSampledFs,Thresholds.(['binarizedWhiskersLower_' strDay]),Thresholds.(['binarizedWhiskersUpper_' strDay]));
        [linkedBinarizedWhiskers] = LinkBinaryEvents_2P_nNOS(gt(binWhisk,0),[round(downSampledFs/3),0]);
        inds = linkedBinarizedWhiskers == 0;
        restAngle = mean(resampledWhisk(inds));
        LabVIEWData.data.dsWhiskerAngle = resampledWhisk - restAngle;
        LabVIEWData.data.binWhiskerAngle = binWhisk;
       
        %% Downsample and binarize the force sensor.
        trimmedForceL = LabVIEWData.data.forceSensor(1:min(expectedLength, length(LabVIEWData.data.forceSensor)));
        % Filter then downsample the Force Sensor waveform to desired frequency
        [z,p,k] = butter(filtOrder,filtThreshold/(LabVIEWData.notes.analogSamplingRate_Hz/2),'low');
        [sos,g] = zp2sos(z,p,k);
        filtForceSensorL = filtfilt(sos,g,trimmedForceL);
        LabVIEWData.data.dsForceSensorL = resample(filtForceSensorL,downSampledFs,LabVIEWData.notes.analogSamplingRate_Hz);
        % Binarize the force sensor waveform
        [ok] = CheckForThreshold_2P_nNOS(['binarizedForceSensor_' strDay],animalID);   
        if ok == 0
            [forceSensorThreshold] = CreateForceSensorThreshold_2P_nNOS(LabVIEWData.data.dsForceSensorL);
            Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
            save([animalID '_Thresholds.mat'],'Thresholds');
        end 
        LabVIEWData.data.binForceSensorL = BinarizeForceSensor_2P_nNOS(LabVIEWData.data.dsForceSensorL,Thresholds.(['binarizedForceSensor_' strDay]));
        
        %% Save the data, set checklist to true
        LabVIEWData.notes.checklist.processData = true;
        save([animalID '_' hem '_' fileID '_LabVIEWData'],'LabVIEWData')
    else
        disp([labviewDataFile ' has already been processed. Continuing...']); disp(' ');     
    end
end

end

