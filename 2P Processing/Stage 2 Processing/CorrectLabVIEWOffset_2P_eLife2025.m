function [] = CorrectLabVIEWOffset_2P_eLife2025(mscanDataFiles,trimTime)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: MScan triggers the LabVIEW acquisition program to start recording, but there is a slight (~ 1 second) lag
%            associated with the MScan data. The force sensor is duplicated, and this function serves to correct that
%            offset by finding the peak in the cross correlation, and shifting the LabVIEW signals based on the number of
%            lags. The beginning/end of all signals are then snipped appropriately after shifting.
%________________________________________________________________________________________________________________________

for a = 1:size(mscanDataFiles,1)
    %% Find offset between the two force sensor signals using the cross correlation
    mscanDataFile = mscanDataFiles(a,:);
    load(mscanDataFile);
    labviewDataFile = MScanData.notes.labviewFileID;
    load(labviewDataFile)
%     if MScanData.notes.checklist.offsetCorrect == false
        disp(['Correcting offset in file number ' num2str(a) ' of ' num2str(size(mscanDataFiles, 1)) '...']); disp(' ');
        [animalID,hem,fileDate,fileID] = GetFileInfo_2P_eLife2025(labviewDataFile);
        imageID = MScanData.notes.imageID;
        vesselID = MScanData.notes.vesselID;
        analogSamplingRate = LabVIEWData.notes.analogSamplingRate_Hz;
        whiskerCamSamplingRate = LabVIEWData.notes.whiskerCamSamplingRate_Hz;
        dsFs = MScanData.notes.dsFs;
        if strcmp(MScanData.notes.movieType,'MC') == true
            lineScanFs = MScanData.data.bloodFlow.Fs;
            desiredFs = 5;
            originalVelocity = abs(MScanData.data.bloodFlow.fixedVelocity);
            % lowpass filter signal below Nyquist frequency of downsample rate;
            [z,p,k] = butter(3,(0.5*floor(desiredFs))/(0.5*lineScanFs),'low');
            [sos,g] = zp2sos(z,p,k);
            filtOriginalVelocity = filtfilt(sos,g,originalVelocity);
            % generate sine wave for resampling
            nonIntegerSineWave = dsp.SineWave(1,desiredFs,1,'SampleRate',lineScanFs,'SamplesPerFrame',length(originalVelocity));
            triggerWave = nonIntegerSineWave();
            waveDiff = diff(triggerWave);
            downSlope = waveDiff < 0;
            edgeFind = diff(downSlope);
            downsampleInds = edgeFind == 1;
            MScanData.data.bloodFlow.dsVelocity = filtOriginalVelocity(:,downsampleInds);
            MScanData.notes.p2Fs = desiredFs;
            vesselSamplingRate = floor(MScanData.notes.p2Fs);
            % check figure
            dsVelocity = figure;
            p1 = plot((1:length(originalVelocity))/lineScanFs,originalVelocity,'r');
            hold on;
            p2 = plot((1:length(MScanData.data.bloodFlow.dsVelocity))/desiredFs,MScanData.data.bloodFlow.dsVelocity,'b');
            xlabel('Time (s)')
            ylabel('Velocity (\muM/sec)')
            legend([p1,p2],'Original signal','Downsampled signal')
        else
            vesselSamplingRate = floor(MScanData.notes.frameRate);
        end
        trialDuration_MScan = MScanData.notes.numberOfFrames/MScanData.notes.frameRate;
        trialDuration_LabVIEW = LabVIEWData.notes.trialDuration_Seconds;
        labviewForce = detrend(LabVIEWData.data.dsForceSensorL,'constant');
        mscanForce = detrend(MScanData.data.dsForceSensorM,'constant');
        sampleDiff = length(labviewForce) - length(mscanForce);
        mscanForce = vertcat(mscanForce,zeros(sampleDiff,1));
        analog_labviewForce = detrend(LabVIEWData.data.forceSensor,'constant');
        analog_mscanForce = detrend(MScanData.data.forceSensor,'constant');
        sampleDiff = length(analog_labviewForce) - length(analog_mscanForce);
        analog_mscanForce = vertcat(analog_mscanForce,zeros(sampleDiff,1));
        analog_labviewSolenoids = LabVIEWData.data.solenoids;
        % analog xcorr
        analog_MaxLag = 30*analogSamplingRate;
        [analog_r,analog_lags] = xcorr(analog_labviewForce,analog_mscanForce,analog_MaxLag);
        [~,analog_index] = max(analog_r);
        analog_offset = analog_lags(analog_index);
        % dsFs xcorr
        maxLag = 30*dsFs;
        [r,lags] = xcorr(labviewForce,mscanForce,maxLag);
        [~,index] = max(r);
        offset = lags(index);
        % offsets
        analog_forceOffset = abs(analog_offset);
        analog_SolenoidOffset = round(abs(analog_offset)/analogSamplingRate);
        analog_whiskerOffset = round(abs(analog_offset)/whiskerCamSamplingRate);
        dsOffset = round(dsFs*(abs(offset)/dsFs));
        disp(['LabVIEW trailed MScan by ' num2str(-offset/dsFs) ' seconds.']); disp(' ')
        if offset > 0
            analog_forceShift = analog_labviewForce(analog_forceOffset:end);
            analog_solenoidShift = analog_labviewSolenoids(analog_SolenoidOffset:end);
            analog_whiskerShift = LabVIEWData.data.whiskerAngle(analog_whiskerOffset:end);
            dsForceShift = labviewForce(offset:end);
            dsWhiskShift = LabVIEWData.data.dsWhiskerAngle(offset:end);
            binForceShift = LabVIEWData.data.binForceSensorL(offset:end);
            binWhiskShift = LabVIEWData.data.binWhiskerAngle(offset:end);
        elseif offset <= 0
            analog_fpad = zeros(1,abs(analog_forceOffset));
            analog_wpad = zeros(1,abs(analog_whiskerOffset));
            dsFs_pad = zeros(1,abs(dsOffset));
            analog_forceShift = horzcat(analog_fpad,analog_labviewForce);
            analog_solenoidShift = horzcat(analog_fpad,analog_labviewSolenoids);
            analog_whiskerShift = horzcat(analog_wpad,LabVIEWData.data.whiskerAngle);
            dsForceShift = horzcat(dsFs_pad,labviewForce);
            dsWhiskShift = horzcat(dsFs_pad,LabVIEWData.data.dsWhiskerAngle);
            binForceShift = horzcat(dsFs_pad,LabVIEWData.data.binForceSensorL);
            binWhiskShift = horzcat(dsFs_pad,LabVIEWData.data.binWhiskerAngle);
        end
        corrOffset = figure;
        ax1 = subplot(3,1,1);
        plot((1:length(mscanForce))/dsFs,mscanForce,'k')
        hold on;
        plot((1:length(labviewForce))/dsFs,labviewForce,'r')
        title({[animalID ' ' fileID ' ' vesselID ' ' imageID ' force sensor data'],'Offset correction between MScan and LabVIEW DAQ'})
        legend('Original MScan','Original LabVIEW')
        ylabel('A.U.')
        xlabel('Time (sec)')
        set(gca,'Ticklength',[0,0])
        axis tight
        ax2 = subplot(3,1,2); %#ok<NASGU>
        plot(analog_lags/dsFs,analog_r,'k')
        title('Cross Correlation between the two signals')
        ylabel('Correlation (A.U.)')
        xlabel('Lag (sec)')
        set(gca,'Ticklength',[0,0])
        axis tight
        ax3 = subplot(3,1,3);
        plot((1:length(mscanForce))/dsFs,mscanForce,'k')
        hold on;
        plot((1:length(dsForceShift))/dsFs,dsForceShift,'b')
        title({'Shifted correction between MScan and LabVIEW DAQ',['Offset value: ' num2str(offset) ' samples or ~' num2str(offset/dsFs) ' seconds']})
        legend('Original MScan','Shifted LabVIEW')
        ylabel('A.U.')
        xlabel('Time (sec)')
        set(gca,'Ticklength',[0,0])
        axis tight
        linkaxes([ax1,ax3],'x')
        %% Apply correction to the data, and trim excess time
        MScan_frontCut = trimTime;
        MScan_endCut = trimTime - (trialDuration_LabVIEW - trialDuration_MScan);
        LabVIEW_frontCut = trimTime;
        LabVIEW_endCut = trimTime;
        mscanAnalogSampleDiff = analogSamplingRate*trialDuration_MScan - length(MScanData.data.forceSensor);
        mscanAnalogCut = floor(MScan_endCut*analogSamplingRate - mscanAnalogSampleDiff);
        mscan_dsAnalogSampleDiff = dsFs*trialDuration_MScan - length(MScanData.data.dsForceSensorM);
        mscan_dsAnalogCut = floor(MScan_endCut*dsFs - mscan_dsAnalogSampleDiff);
        mscan_binForceSampleDiff = dsFs*trialDuration_MScan - length(MScanData.data.binForceSensorM);
        mscan_binForceCut = floor(MScan_endCut*dsFs - mscan_binForceSampleDiff);
        labview_AnalogSampleDiff = analogSamplingRate*trialDuration_LabVIEW - length(analog_forceShift);
        labview_AnalogCut = floor(LabVIEW_endCut*analogSamplingRate - labview_AnalogSampleDiff);
        labview_WhiskerSampleDiff = whiskerCamSamplingRate*trialDuration_LabVIEW - length(analog_whiskerShift);
        labview_WhiskerCut = floor(LabVIEW_endCut*whiskerCamSamplingRate - labview_WhiskerSampleDiff);
        labview_dsWhiskSamplingDiff = dsFs*trialDuration_LabVIEW - length(dsWhiskShift);
        labview_dsWhiskCut = floor(LabVIEW_endCut*dsFs - labview_dsWhiskSamplingDiff);
        labview_dsForceSamplingDiff = dsFs*trialDuration_LabVIEW - length(dsForceShift);
        labview_dsForceCut = floor(LabVIEW_endCut*dsFs - labview_dsForceSamplingDiff);
        labview_binForceSampleDiff = dsFs*trialDuration_LabVIEW - length(binForceShift);
        labview_binForceCut = floor(LabVIEW_endCut*dsFs - labview_binForceSampleDiff);
        labview_binWhiskSamplingDiff = dsFs*trialDuration_LabVIEW - length(binWhiskShift);
        labview_binWhiskCut = floor(LabVIEW_endCut*dsFs - labview_binWhiskSamplingDiff);
        MScanData.data.forceSensor_trim = MScanData.data.forceSensor(floor(MScan_frontCut*analogSamplingRate):end - (mscanAnalogCut + 1))';
        MScanData.data.corticalNeural_trim = MScanData.data.corticalNeural(floor(MScan_frontCut*analogSamplingRate):end - (mscanAnalogCut + 1))';
        MScanData.data.hippocampalNeural_trim = MScanData.data.hippocampalNeural(floor(MScan_frontCut*analogSamplingRate):end - (mscanAnalogCut + 1))';
        MScanData.data.EMG_trim = MScanData.data.EMG(floor(MScan_frontCut*analogSamplingRate):end - (mscanAnalogCut + 1))';
        if strcmp(MScanData.notes.movieType,'MC') == true
            MScanData.data.vesselDiameter_trim = MScanData.data.bloodFlow.dsVelocity(floor(MScan_frontCut*vesselSamplingRate) - 1:end - (floor(MScan_endCut*vesselSamplingRate)));
        else
            MScanData.data.vesselDiameter_trim = MScanData.data.vesselDiameter(floor(MScan_frontCut*vesselSamplingRate):end - (MScan_endCut*vesselSamplingRate + 1));
        end
        MScanData.data.cortical.muaPower_trim = MScanData.data.cortical.muaPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.cortical.gammaBandPower_trim = MScanData.data.cortical.gammaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.cortical.betaBandPower_trim = MScanData.data.cortical.betaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.cortical.alphaBandPower_trim = MScanData.data.cortical.alphaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.cortical.thetaBandPower_trim = MScanData.data.cortical.thetaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.cortical.deltaBandPower_trim = MScanData.data.cortical.deltaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.hippocampal.muaPower_trim = MScanData.data.hippocampal.muaPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.hippocampal.gammaBandPower_trim = MScanData.data.hippocampal.gammaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.hippocampal.betaBandPower_trim = MScanData.data.hippocampal.betaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.hippocampal.alphaBandPower_trim = MScanData.data.hippocampal.alphaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.hippocampal.thetaBandPower_trim = MScanData.data.hippocampal.thetaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.hippocampal.deltaBandPower_trim = MScanData.data.hippocampal.deltaBandPower(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.dsForceSensorM_trim = MScanData.data.dsForceSensorM(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.data.binForceSensorM_trim = MScanData.data.binForceSensorM(floor(MScan_frontCut*dsFs):end - (mscan_binForceCut + 1))';
        MScanData.data.filtEMG_trim = MScanData.data.filtEMG(floor(MScan_frontCut*dsFs):end - (mscan_dsAnalogCut + 1))';
        MScanData.notes.shiftLags = analog_lags;
        MScanData.notes.shiftXCorr = analog_r;
        MScanData.notes.checklist.offsetCorrect = true;
        LabVIEWData.data.solenoids_trim = analog_solenoidShift(floor(LabVIEW_frontCut*analogSamplingRate):end - (labview_AnalogCut + 1));
        LabVIEWData.data.forceSensor_trim = analog_forceShift(floor(LabVIEW_frontCut*analogSamplingRate):end - (labview_AnalogCut + 1));
        LabVIEWData.data.whiskerAngle_trim = analog_whiskerShift(floor(LabVIEW_frontCut*whiskerCamSamplingRate):end - (labview_WhiskerCut + 1));
        LabVIEWData.data.dsWhiskerAngle_trim = dsWhiskShift(floor(LabVIEW_frontCut*dsFs):end - (labview_dsWhiskCut + 1));
        LabVIEWData.data.binWhiskerAngle_trim = binWhiskShift(floor(LabVIEW_frontCut*dsFs):end - (labview_binWhiskCut + 1));
        LabVIEWData.data.dsForceSensorL_trim = dsForceShift(floor(LabVIEW_frontCut*dsFs):end - (labview_dsForceCut + 1));
        LabVIEWData.data.binForceSensorL_trim = binForceShift(floor(LabVIEW_frontCut*dsFs):end - (labview_binForceCut + 1));
        LabVIEWData.notes.checklist.offsetCorrect = true;
        LabVIEWData.notes.trimTime = trimTime;
        LabVIEWData.notes.trialDuration_Seconds_trim = LabVIEWData.notes.trialDuration_Seconds - 2*trimTime;
        save([animalID '_' fileDate '_' imageID '_' vesselID '_MScanData'],'MScanData')
        save([animalID '_' hem '_' fileID '_LabVIEWData'],'LabVIEWData')
        %% Save the file to directory.
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/XCorr Shift/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(corrOffset,[dirpath animalID '_' fileID '_' imageID '_' vesselID '_XCorrShift']);
        close(corrOffset)
        if strcmp(MScanData.notes.movieType,'MC') == true
            dirpath = [pathstr '/Figures/Velocity Downsample/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(dsVelocity,[dirpath animalID '_' fileID '_' imageID '_' vesselID '_dsVelocityShift']);
            close(dsVelocity)
        end
%     else
%         disp(['Offset in ' mscanDataFile ' and ' labviewDataFile ' has already been corrected. Continuing...']); disp(' ');
%     end
end

end
