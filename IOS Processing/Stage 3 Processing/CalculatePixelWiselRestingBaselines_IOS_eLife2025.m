function [RestingBaselines] = CalculatePixelWiselRestingBaselines_IOS_eLife2025(procDataFileIDs,RestingBaselines)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
warning('off','imageio:tifftagsread:nextIfdPointerOutOfRange')
if isfield(RestingBaselines,'Pixel') == false
    [animalID,~,~] = GetFileInfo_IOS_eLife2025(procDataFileIDs(1,:));
    restFileList = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs); % obtain the list of unique fileIDs
    % obtain the information from all the resting files
    for bb = 1:length(restFileList)
        fileID = restFileList{bb,1}; % FileID of currently loaded file
        load([animalID '_' fileID '_ProcData.mat']);
        samplingRate = ProcData.notes.wavelengthSamplingRate;
        imagingWavelengths = ProcData.notes.imagingWavelengths;
        imagingCamera = ProcData.notes.iosCamera;
        strDay = ConvertDate_IOS_eLife2025(fileID(1:6));
        % load in neural data from current file
        disp(['Reading frames from file ' num2str(bb)]); disp(' ');
        if strcmp(imagingCamera,'PCO Edge 5.5') == true
            info = imfinfo([fileID '_PCO_Cam01.pcoraw']);
            numberOfPages = length(info);
            for k = 1:numberOfPages
                imageStack.full(:,:,k) = double(imread([fileID '_PCO_Cam01.pcoraw'],k));
            end
        elseif strcmp(imagingCamera,'Dalsa Pantera TF 1M60') == true
            imageWidth = ProcData.notes.CBVCamPixelWidth;
            imageHeight = ProcData.notes.CBVCamPixelHeight;
            pixelsPerFrame = imageWidth*imageHeight;
            % open the file, get file size, back to the begining
            fid = fopen([fileID '_WindowCam.bin']);
            fseek(fid,0,'eof');
            fseek(fid,0,'bof');
            % identify the number of frames to read. Each frame has a previously defined width and height (as inputs), along with a grayscale "depth" of 2"
            nFramesToRead = 10;
            % pre-allocate memory
            for n = 1:nFramesToRead
                z = fread(fid,pixelsPerFrame,'*int16','b');
                img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
                imageStack.full(:,:,n) = double(rot90(img',2));
            end
        end
        % separate image stack by wavelength
        if strcmp(imagingWavelengths,{'Red, Green, & Blue'}) == true
            wavelengths = {'red','green','blue'};
            imageStack.red = imageStack.full(:,:,2:3:end);
            imageStack.green = imageStack.full(:,:,3:3:end);
            imageStack.blue = imageStack.full(:,:,1:3:end - 1);
        elseif strcmp(imagingWavelengths,{'Red, Lime, & Blue'}) == true
            wavelengths = {'red','lime','blue'};
            imageStack.red = imageStack.full(:,:,2:3:end);
            imageStack.lime = imageStack.full(:,:,3:3:end);
            imageStack.blue = imageStack.full(:,:,1:3:end - 1);
        elseif strcmp(imagingWavelengths,{'Green & Blue'}) == true
            wavelengths = {'green','blue'};
            imageStack.green = imageStack.full(:,:,2:2:end);
            imageStack.blue = imageStack.full(:,:,1:2:end - 1);
        elseif strcmp(imagingWavelengths,{'Lime & Blue'}) == true
            wavelengths = {'lime','blue'};
            imageStack.lime = imageStack.full(:,:,2:2:end);
            imageStack.blue = imageStack.full(:,:,1:2:end - 1);
        elseif strcmp(imagingWavelengths,'Green') == true
            wavelengths = {'green'};
            imageStack.green = imageStack.full(:,:,1:end);
        elseif strcmp(imagingWavelengths,'Lime') == true
            wavelengths = {'lime'};
            imageStack.lime = imageStack.full(:,:,1:end);
        elseif strcmp(imagingWavelengths,'Blue') == true
            wavelengths = {'blue'};
            imageStack.blue = imageStack.full(:,:,1:end);
        end
        for zz = 1:length(wavelengths)
            wavelength = wavelengths{1,zz};
            restFileData = [];
            for cc = 1:length(RestingBaselines.manualSelection.baselineFileInfo.fileIDs)
                restFileID = RestingBaselines.manualSelection.baselineFileInfo.fileIDs{cc,1};
                if strcmp(fileID,restFileID)
                    restDuration = round(RestingBaselines.manualSelection.baselineFileInfo.durations(cc,1),1);
                    startTime = round(RestingBaselines.manualSelection.baselineFileInfo.eventTimes(cc,1),1);
                    % conditions and indexing
                    startTimeIndex = floor(startTime*samplingRate);
                    restDurationIndex = floor(restDuration*samplingRate - 1);
                    restEventData = mean(imageStack.(wavelength)(:,:,(startTimeIndex:(startTimeIndex + restDurationIndex))),3);
                    if sum(sum(isnan(restEventData))) == 0
                        restFileData = cat(3,restFileData,restEventData);
                    end
                end
            end
            trialRestData.([strDay '_' fileID]).(wavelength) = restFileData;
            clear restFileData restEventData
        end
        clear imageStack
    end
    fields = fieldnames(trialRestData);
    uniqueDays = GetUniqueDays_IOS_eLife2025(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
    for aa = 1:length(wavelengths)
        wavelength = wavelengths{1,aa};
        clear frameAvgs
        for f = 1:length(uniqueDays)
            g = 1;
            stringDay = ConvertDate_IOS_eLife2025(uniqueDays{f});
            frameAvgs.(wavelength).(stringDay) = [];
            for field = 1:length(fields)
                if strcmp(fields{field}(7:12),uniqueDays{f})
                    frameAvgs.(wavelength).(stringDay) = cat(3,frameAvgs.(wavelength).(stringDay),trialRestData.(fields{field}).(wavelength));
                    g = g + 1;
                end
            end
        end
        dayFields = fieldnames(frameAvgs.(wavelength));
        for h = 1:length(dayFields)
            disp(['Adding ' wavelength ' pixel-wise baseline to baseline file for ' dayFields{h} '...']); disp(' ')
            % for bb = 1:length(wavelengths)
            RestingBaselines.Pixel.(wavelength).(dayFields{h}) = mean(frameAvgs.(wavelength).(dayFields{h}),3);
        end
    end
    save([animalID '_RestingBaselines.mat'],'RestingBaselines')
end