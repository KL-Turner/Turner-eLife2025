function [] = ExtractHeartRate_IOS(procDataFileIDs)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
load(procDataFileIDs(1,:));
imagingType = ProcData.notes.imagingType; %#ok<NODEF>
imagingWavelengths = ProcData.notes.imagingWavelengths;
if any(strcmp(imagingWavelengths,{'Green','Blue','Green & Blue','Red, Green, & Blue'})) == true
    wavelength = 'green';
elseif any(strcmp(imagingWavelengths,{'Lime','Lime & Blue','Red, Lime, & Blue'})) == true
    wavelength = 'lime';
end
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    load(procDataFileID)
    if isfield(ProcData.data,'heartRate') == false
        disp(['Extracting heart rate from file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
        if ProcData.notes.wavelengthSamplingRate >= 30
            if strcmpi(imagingType,'Single ROI (SI)') == true
                [~,~,~,HR] = FindHeartRate_IOS(ProcData.data.(wavelength).barrels,ProcData.notes.SamplingRate);
            elseif strcmpi(imagingType,'Single ROI (SSS)') == true
                [~,~,~,HR] = FindHeartRate_IOS(ProcData.data.(wavelength).SSS,ProcData.notes.wavelengthSamplingRate);
            elseif strcmpi(imagingType,'Bilateral ROI (SI)') == true || strcmpi(imagingType,'Bilateral ROI (SI,FC)') == true
                % pull out the left and right window heart rate. they should be essentiall6 identical
                [~,~,~,LH_HR] = FindHeartRate_IOS(ProcData.data.(wavelength).LH,ProcData.notes.wavelengthSamplingRate);
                [~,~,~,RH_HR] = FindHeartRate_IOS(ProcData.data.(wavelength).RH,ProcData.notes.wavelengthSamplingRate);
                % average the two signals from the left and right windows
                HR = (LH_HR + RH_HR)/2;
            end
            % patch the missing data at the beginning and end of the signal
            patchedHR = horzcat(HR(1),HR,HR(end),HR(end));
            % smooth the signal with a 2 Hz low pass third-order butterworth filter
            [B,A] = butter(3,2/(ProcData.notes.wavelengthSamplingRate/2),'low');
            heartRate = filtfilt(B,A,patchedHR); % filtered heart rate signal
            ProcData.data.heartRate.frequency = heartRate;
            save(procDataFileID,'ProcData');
        else
            ProcData.data.heartRate.frequency = zeros(1,ProcData.notes.trialDuration_sec);
            save(procDataFileID,'ProcData');
        end
    else
        disp(['Heart rate for file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ') already analyzed.']); disp(' ')
    end
end