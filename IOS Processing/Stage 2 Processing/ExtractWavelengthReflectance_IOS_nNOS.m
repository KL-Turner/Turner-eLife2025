function [] = ExtractWavelengthReflectance_IOS_nNOS(ROIs,ROInames,procDataFileIDs,imagingWavelengths,imagingCamera)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
warning('off','imageio:tifftagsread:nextIfdPointerOutOfRange')
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    [~,fileDate,fileID] = GetFileInfo_IOS_nNOS(procDataFileID);
    strDay = ConvertDate_IOS_nNOS(fileDate);
    load(procDataFileID)
    if isfield(ProcData.data,'green') == false
        disp(['Analyzing IOS reflectance from file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
        if strcmp(imagingCamera,'PCO Edge 5.5') == true
            info = imfinfo([fileID '_PCO_Cam01.pcoraw']);
            numberOfPages = length(info);
            for k = 1:numberOfPages
                imageStack(:,:,k) = double(imread([fileID '_PCO_Cam01.pcoraw'],k));
            end
        elseif strcmp(imagingCamera,'Dalsa Pantera TF 1M60') == true
            imageHeight = ProcData.notes.CBVCamPixelHeight;
            imageWidth = ProcData.notes.CBVCamPixelWidth;
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
                imageStack(:,:,n) = double(rot90(img',2));
            end
        end
        % separate image stack by wavelength
        if strcmp(imagingWavelengths,{'Red, Green, & Blue'}) == true
            wavelengths = {'red','green','blue'};
            ProcData.notes.wavelengthSamplingRate = ProcData.notes.CBVCamSamplingRate/3;
            data.red = imageStack(:,:,2:3:end);
            data.green = imageStack(:,:,3:3:end);
            data.blue = imageStack(:,:,1:3:end - 1);
        elseif strcmp(imagingWavelengths,{'Red, Lime, & Blue'}) == true
            wavelengths = {'red','lime','blue'};
            ProcData.notes.wavelengthSamplingRate = ProcData.notes.CBVCamSamplingRate/3;
            data.red = imageStack(:,:,2:3:end);
            data.lime = imageStack(:,:,3:3:end);
            data.blue = imageStack(:,:,1:3:end - 1);
        elseif strcmp(imagingWavelengths,{'Green & Blue'}) == true
            ProcData.notes.wavelengthSamplingRate = ProcData.notes.CBVCamSamplingRate/2;
            wavelengths = {'green','blue'};
            data.green = imageStack(:,:,2:2:end);
            data.blue = imageStack(:,:,1:2:end - 1);
        elseif strcmp(imagingWavelengths,{'Lime & Blue'}) == true
            ProcData.notes.wavelengthSamplingRate = ProcData.notes.CBVCamSamplingRate/2;
            wavelengths = {'lime','blue'};
            data.lime = imageStack(:,:,2:2:end);
            data.blue = imageStack(:,:,1:2:end - 1);
        elseif strcmp(imagingWavelengths,'Green') == true
            ProcData.notes.wavelengthSamplingRate = ProcData.notes.CBVCamSamplingRate;
            wavelengths = {'green'};
            data.green = imageStack(:,:,1:end);
        elseif strcmp(imagingWavelengths,'Lime') == true
            ProcData.notes.wavelengthSamplingRate = ProcData.notes.CBVCamSamplingRate;
            wavelengths = {'lime'};
            data.lime = imageStack(:,:,1:end);
        elseif strcmp(imagingWavelengths,'Blue') == true
            ProcData.notes.wavelengthSamplingRate = ProcData.notes.CBVCamSamplingRate;
            wavelengths = {'blue'};
            data.blue = imageStack(:,:,1:end);
        end
        for b = 1:length(ROInames)
            ROIshortName = ROInames{1,b};
            ROIname = ROInames{1,b};
            for cc = 1:length(wavelengths)
                reflField = wavelengths{1,cc};
                if any(strcmp(ROIshortName,{'LH','RH','fLH','fRH','barrels'})) == true
                    maskFig = figure;
                    imagesc(data.(reflField)(:,:,1));
                    axis image;
                    colormap gray
                    circROI = drawcircle('Center',ROIs.(strDay).(ROIname).circPosition,'Radius',ROIs.(strDay).(ROIname).circRadius);
                    ROImask = createMask(circROI,data.(reflField)(:,:,1));
                    close(maskFig)
                    for aa = 1:size(data.(reflField),3)
                        mask = ROImask.*double(data.(reflField)(:,:,aa));
                        refl(aa,1) = mean(nonzeros(mask));
                    end
                    ProcData.data.(reflField).(ROIname) = refl';
                else
                    maskFig = figure;
                    imagesc(data.(reflField)(:,:,1));
                    axis image;
                    colormap gray
                    ROImask = roipoly(double(data.(reflField)(:,:,1)),ROIs.(ROIname).xi,ROIs.(ROIname).yi);
                    close(maskFig)
                    for aa = 1:size(data.(reflField),3)
                        mask = ROImask.*double(data.(reflField)(:,:,aa));
                        refl(aa,1) = mean(nonzeros(mask));
                    end
                    ProcData.data.(reflField).(ROIname) = refl';
                end
            end
        end
        save(procDataFileID,'ProcData')
    else
        disp(['IOS reflectance for file (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ') already analyzed.']); disp(' ')
    end
end