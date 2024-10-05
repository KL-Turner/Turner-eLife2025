function [ROIs] = PlaceROIs_IOS(animalID,fileID,ROIs,lensMag,imagingType,imagingWavelengths,imagingCamera)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
strDay = ConvertDate_IOS(fileID);
fileDate = fileID(1:6);
% determine which ROIs to draw based on imaging type
if strcmpi(imagingType,'Single ROI (SI)') == true
    ROInames = {'barrels'};
elseif strcmpi(imagingType,'Single ROI (SSS)') == true
    ROInames = {'SSS','lSSS','rSSS'};
elseif strcmpi(imagingType,'Bilateral ROI (SI)') == true
    ROInames = {'LH','RH'};
elseif strcmpi(imagingType,'Bilateral ROI (SI,FC)') == true
    ROInames = {'LH','RH','fLH','fRH'};
end
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
if strcmp(imagingCamera,'PCO Edge 5.5') == true
    % character list of all WindowCam files
    windowDataFileStruct = dir('*_PCO_Cam01.pcoraw');
    windowDataFiles = {windowDataFileStruct.name}';
    windowDataFileIDs = char(windowDataFiles);
    bb = 1;
    % double check file lists
    for aa = 1:size(procDataFileIDs)
        procDataFileID = procDataFileIDs(aa,:);
        [~,procfileDate,~] = GetFileInfo_IOS(procDataFileID);
        windowDataFileID = windowDataFileIDs(aa,:);
        if strcmp(procfileDate,fileDate) == true
            windowDataFileList(bb,:) = windowDataFileID;
            bb = bb + 1;
        end
    end
    nFramesToRead = 10;
    % pre-allocate memory
    frames = cell(1,nFramesToRead);
    for n = 1:nFramesToRead
        frames{n} = double(imread(windowDataFileList(1,:),n));
    end
elseif strcmp(imagingCamera,'Dalsa Pantera TF 1M60')
    % character list of all WindowCam files
    windowDataFileStruct = dir('*_WindowCam.bin');
    windowDataFiles = {windowDataFileStruct.name}';
    windowDataFileIDs = char(windowDataFiles);
    bb = 1;
    % double check file lists
    for aa = 1:size(procDataFileIDs)
        procDataFileID = procDataFileIDs(aa,:);
        [~,procfileDate,~] = GetFileInfo_IOS(procDataFileID);
        windowDataFileID = windowDataFileIDs(aa,:);
        if strcmp(procfileDate,fileDate) == true
            windowDataFileList(bb,:) = windowDataFileID;
            bb = bb + 1;
        end
    end
    imageHeight = ProcData.notes.CBVCamPixelHeight;
    imageWidth = ProcData.notes.CBVCamPixelWidth;
    pixelsPerFrame = imageWidth*imageHeight;
    % open the file, get file size, back to the begining
    fid = fopen(windowDataFileList(1,:));
    fseek(fid,0,'eof');
    fseek(fid,0,'bof');
    % identify the number of frames to read. Each frame has a previously defined width and height (as inputs), along with a grayscale "depth" of 2"
    nFramesToRead = 10;
    % pre-allocate memory
    frames = cell(1,nFramesToRead);
    for n = 1:nFramesToRead
        z = fread(fid,pixelsPerFrame,'*int16','b');
        img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
        frames{n} = double(rot90(img',2));
    end
end
if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Red, Lime, & Blue'})) == true
    roiFrame = frames{3};
elseif any(strcmp(imagingWavelengths,{'Green & Blue','Lime & Blue'})) == true
    roiFrame = frames{2};
else
    roiFrame = frames{1};
end
% determine the proper size of the ROI based on camera/lens magnification
if strcmpi(lensMag,'0.75X') == true
    circRadius = 37/2; % pixels to be 1 mm in diameter
elseif strcmpi(lensMag,'1.0X') == true
    circRadius = 45/2;
elseif strcmpi(lensMag,'1.5X') == true
    circRadius = 60/2;
elseif strcmpi(lensMag,'2.0X') == true
    circRadius = 75/2;
elseif strcmpi(lensMag,'2.5X') == true
    circRadius = 90/2;
elseif strcmpi(lensMag,'3.0X') == true
    circRadius = 105/2;
end
% place circle along the most relevant region of each hemisphere
for ff = 1:length(ROInames)
    % generate image
    isok = false;
    while isok == false
        windowFig = figure;
        imagesc(roiFrame)
        title([animalID ' ' ROInames{1,ff} ' ROI'])
        xlabel('Image size (pixels)')
        ylabel('Image size (pixels)')
        colormap gray
        axis image
        disp(['Move the ROI over the desired region for ' ROInames{1,ff}]); disp(' ')
        drawnow
        circ = drawcircle('Center',[0,0],'Radius',circRadius,'Color','r');
        checkCircle = input('Is the ROI okay? (y/n): ','s'); disp(' ')
        circPosition = round(circ.Center);
        if strcmpi(checkCircle,'y') == true
            isok = true;
            ROIs.(strDay).(ROInames{1,ff}).circPosition = circPosition;
            ROIs.(strDay).(ROInames{1,ff}).circRadius = circRadius;
        end
        delete(windowFig);
    end
end
% check final image
fig = figure;
imagesc(roiFrame)
hold on;
for aa = 1:length(ROInames)
    drawcircle('Center',ROIs.(strDay).(ROInames{1,aa}).circPosition,'Radius',ROIs.(strDay).(ROInames{1,aa}).circRadius,'Color','r');
end
title([animalID ' final ROI placement'])
xlabel('Image size (pixels)')
ylabel('Image size (pixels)')
colormap gray
axis image
% save the file to directory.
[pathstr,~,~] = fileparts(cd);
dirpath = [pathstr '/Figures/IOS ROIs/'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(fig,[dirpath animalID '_' strDay '_ROIs.fig'])