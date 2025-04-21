function [] = CalculateFullSpectroscopy_IOS_eLife2025(procDataFileIDs,RestingBaselines)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner, adapted from code written by Qingguang Zhang
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
warning('off','imageio:tifftagsread:nextIfdPointerOutOfRange')
% load structure with ROI locations
ROIFileDir = dir('*_ROIs.mat');
ROIFileName = {ROIFileDir.name}';
ROIFileID = char(ROIFileName);
load(ROIFileID);
% pathlengths for each wavelength (cm)
X.red = 0.3846483;    % 630 nm
X.green = 0.0371713;  % 530 nm
X.lime = 0.0324347;   % 570 nm
X.blue = 0.064898;    % 480 nm
% extinction coefficients (cm^-1*M^-1)
E_HbO.red = 610;      % 630 nm
E_HbR.red = 5148.8;
E_HbO.green = 39500;  % 530 nm
E_HbR.green = 39500;
E_HbO.lime = 45000;   % 570 nm
E_HbR.lime = 45000;
E_HbO.blue = 14550;   % 480 nm
E_HbR.blue = 26629.2;
%% go through each file and do pixel-wise analysis
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    imagingCamera = ProcData.notes.iosCamera;
    imagingWavelengths = ProcData.notes.imagingWavelengths;
    if isfield(ProcData.data,'HbT') == false
        disp(['Adding pixel-wise spectroscopy to file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')...']); disp(' ')
        [~,fileDate,fileID] = GetFileInfo_IOS_eLife2025(procDataFileID);
        strDay = ConvertDate_IOS_eLife2025(fileDate);
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
        conv2um = 1e6;
        % split images R-G-B
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
        %% HbT analysis
        if any(strcmp(imagingWavelengths,{'Green','Green & Blue','Red, Green, & Blue'})) == true
            ledType = 'M530L3';
            bandfilterType = 'FB530-10';
            cutfilterType = 'MF525-39';
            [~,~,weightedcoeffHbT] = GetHbcoeffs_IOS_eLife2025(ledType,bandfilterType,cutfilterType);
            HbT_ImageStack = (log(imageStack.green./RestingBaselines.Pixel.green.(strDay))).*weightedcoeffHbT*conv2um;
            roiImage = imageStack.green(:,:,1);
        elseif any(strcmp(imagingWavelengths,{'Lime','Lime & Blue','Red, Lime, & Blue'})) == true
            ledType = 'M565L3';
            bandfilterType = 'FB570-10';
            cutfilterType = 'EO65160';
            [~,~,weightedcoeffHbT] = GetHbcoeffs_IOS_eLife2025(ledType,bandfilterType,cutfilterType);
            HbT_ImageStack = (log(imageStack.lime./RestingBaselines.Pixel.lime.(strDay))).*weightedcoeffHbT*conv2um;
            roiImage = imageStack.lime(:,:,1);
        else
            weightedcoeffHbT = NaN;
            HbT_ImageStack = (log(imageStack.blue./RestingBaselines.Pixel.blue.(strDay))).*weightedcoeffHbT*conv2um;
            roiImage = imageStack.blue(:,:,1);
        end
        % extract pixel-wise ROIs
        ROInames = fieldnames(ROIs.(strDay));
        for dd = 1:length(ROInames)
            roiName = ROInames{dd,1};
            maskFig = figure;
            imagesc(roiImage);
            axis image;
            colormap gray
            circROI = drawcircle('Center',ROIs.(strDay).(roiName).circPosition,'Radius',ROIs.(strDay).(roiName).circRadius);
            mask = createMask(circROI,roiImage);
            close(maskFig)
            for cc = 1:size(HbT_ImageStack,3)
                HbT_roiMask = mask.*double(HbT_ImageStack(:,:,cc));
                ProcData.data.HbT.(roiName)(1,cc) = mean(nonzeros(HbT_roiMask),'omitnan');
            end
        end
        %% HbO/HbR analysis
        if any(strcmp(wavelengths,'red')) == true
            for bb = 1:length(wavelengths)
                wavelength = wavelengths{1,bb};
                % normalize each stack by resting baseline frame
                normImageStack.(wavelength) = ((imageStack.(wavelength) - RestingBaselines.Pixel.(wavelength).(strDay))./RestingBaselines.Pixel.(wavelength).(strDay)) + 1;
                % calculate absorption coefficient for each pixel = -1/pathlength*ln(deltaR/R+1)
                Mu.(wavelength) = -1/X.(wavelength)*log(normImageStack.(wavelength)); % cm^-1, natural logarithm
            end
            % Calculate concentrations (uM) using green/lime and red light
            if any(strcmp(wavelengths,'green')) == true
                HbR_ImageStack = (E_HbO.green*Mu.red - E_HbO.red*Mu.green)./(E_HbO.green*E_HbR.red - E_HbO.red*E_HbR.green)*conv2um; % in uM
                HbO_ImageStack = (E_HbR.green*Mu.red - E_HbR.red*Mu.green)./(E_HbR.green*E_HbO.red - E_HbR.red*E_HbO.green)*conv2um; % in uM
                % extract pixel-wise ROIs
                ROInames = fieldnames(ROIs.(strDay));
                for dd = 1:length(ROInames)
                    roiName = ROInames{dd,1};
                    maskFig = figure;
                    imagesc(imageStack.green(:,:,1));
                    axis image;
                    colormap gray
                    circROI = drawcircle('Center',ROIs.(strDay).(roiName).circPosition,'Radius',ROIs.(strDay).(roiName).circRadius);
                    mask = createMask(circROI,imageStack.green(:,:,1));
                    close(maskFig)
                    for cc = 1:size(HbO_ImageStack,3)
                        HbO_roiMask = mask.*double(HbO_ImageStack(:,:,cc));
                        ProcData.data.HbO.(roiName)(1,cc) = mean(nonzeros(HbO_roiMask),'omitnan');
                        HbR_roiMask = mask.*double(HbR_ImageStack(:,:,cc));
                        ProcData.data.HbR.(roiName)(1,cc) = mean(nonzeros(HbR_roiMask),'omitnan');
                    end
                end
            elseif any(strcmp(wavelengths,'lime')) == true
                HbR_ImageStack = (E_HbO.lime*Mu.red - E_HbO.red*Mu.lime)./(E_HbO.lime*E_HbR.red - E_HbO.red*E_HbR.lime)*conv2um; % in uM
                HbO_ImageStack = (E_HbR.lime*Mu.red - E_HbR.red*Mu.lime)./(E_HbR.lime*E_HbO.red - E_HbR.red*E_HbO.lime)*conv2um; % in uM
                % extract pixel-wise ROIs
                ROInames = fieldnames(ROIs.(strDay));
                for dd = 1:length(ROInames)
                    roiName = ROInames{dd,1};
                    maskFig = figure;
                    imagesc(imageStack.lime(:,:,1));
                    axis image;
                    colormap gray
                    circROI = drawcircle('Center',ROIs.(strDay).(roiName).circPosition,'Radius',ROIs.(strDay).(roiName).circRadius);
                    mask = createMask(circROI,imageStack.lime(:,:,1));
                    close(maskFig)
                    for cc = 1:size(HbO_ImageStack,3)
                        HbO_roiMask = mask.*double(HbO_ImageStack(:,:,cc));
                        ProcData.data.HbO.(roiName)(1,cc) = mean(nonzeros(HbO_roiMask),'omitnan');
                        HbR_roiMask = mask.*double(HbR_ImageStack(:,:,cc));
                        ProcData.data.HbR.(roiName)(1,cc) = mean(nonzeros(HbR_roiMask),'omitnan');
                    end
                end
            end
        end
        %% correct hemodynamic attenuation of GCaMP fluorescence
        if any(strcmp(imagingWavelengths,{'Green & Blue','Lime & Blue','Red, Green, & Blue','Red, Lime, & Blue'})) == true
            if any(strcmp(wavelengths,'green')) == true
                scale = RestingBaselines.Pixel.green.(strDay)./RestingBaselines.Pixel.blue.(strDay);
                correctedGCaMP_ImageStack = (imageStack.blue./imageStack.green).*scale;
                % extract pixel-wise ROIs
                ROInames = fieldnames(ROIs.(strDay));
                for dd = 1:length(ROInames)
                    roiName = ROInames{dd,1};
                    maskFig = figure;
                    imagesc(imageStack.green(:,:,1));
                    axis image;
                    colormap gray
                    circROI = drawcircle('Center',ROIs.(strDay).(roiName).circPosition,'Radius',ROIs.(strDay).(roiName).circRadius);
                    mask = createMask(circROI,imageStack.green(:,:,1));
                    close(maskFig)
                    for cc = 1:size(correctedGCaMP_ImageStack,3)
                        GCaMP_roiMask = mask.*double(correctedGCaMP_ImageStack(:,:,cc));
                        ProcData.data.GCaMP.(roiName)(1,cc) = mean(nonzeros(GCaMP_roiMask),'omitnan');
                    end
                end
            elseif any(strcmp(wavelengths,'lime')) == true
                scale = RestingBaselines.Pixel.lime.(strDay)./RestingBaselines.Pixel.blue.(strDay);
                correctedGCaMP_ImageStack = (imageStack.blue./imageStack.lime).*scale;
                % extract pixel-wise ROIs
                ROInames = fieldnames(ROIs.(strDay));
                for dd = 1:length(ROInames)
                    roiName = ROInames{dd,1};
                    maskFig = figure;
                    imagesc(imageStack.lime(:,:,1));
                    axis image;
                    colormap gray
                    circROI = drawcircle('Center',ROIs.(strDay).(roiName).circPosition,'Radius',ROIs.(strDay).(roiName).circRadius);
                    mask = createMask(circROI,imageStack.lime(:,:,1));
                    close(maskFig)
                    for cc = 1:size(correctedGCaMP_ImageStack,3)
                        GCaMP_roiMask = mask.*double(correctedGCaMP_ImageStack(:,:,cc));
                        ProcData.data.GCaMP.(roiName)(1,cc) = mean(nonzeros(GCaMP_roiMask),'omitnan');
                    end
                end
            end
        end
        % save data
        save(procDataFileID,'ProcData')
    end
end