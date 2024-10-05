function [PupilData] = SelectPupilFiles_IOS(trialData)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
pupilFileDir = dir('*_PupilData.mat');
if isempty(pupilFileDir) == true
    PupilData = [];
    PupilData.EyeROI = [];
else
    pupilFileName = {pupilFileDir.name}';
    pupilFileID = char(pupilFileName);
    load(pupilFileID);
end
% pupil file IDs
pupilFileDir = dir('*_PupilCam.bin');
pupilFileNames = {pupilFileDir.name}';
pupilFileIDs = char(pupilFileNames);
% check first file of each day
if isfield(PupilData,'firstFileOfDay') == false
    % establish the number of unique days based on file IDs
    [~,fileDates,~] = GetFileInfo_IOS(pupilFileIDs);
    [uniqueDays,~,DayID] = GetUniqueDays_IOS(fileDates);
    for aa = 1:length(uniqueDays)
        strDay = ConvertDate_IOS(uniqueDays{aa,1});
        FileInd = DayID == aa;
        dayFilenames.(strDay) = pupilFileIDs(FileInd,:);
    end
    for bb = 1:length(uniqueDays)
        strDay = ConvertDate_IOS(uniqueDays{bb,1});
        for cc = 1:size(dayFilenames.(strDay),1)
            pupilCamFileID = dayFilenames.(strDay)(cc,:);
            fid = fopen(pupilCamFileID); % reads the binary file in to the work space
            fseek(fid,0,'eof'); % find the end of the video frame
            imageHeight = str2double(trialData.pupilCamPixelHeight);
            imageWidth = str2double(trialData.pupilCamPixelWidth);
            pixelsPerFrame = imageWidth*imageHeight;
            skippedPixels = pixelsPerFrame;
            roiImage = zeros(imageHeight,imageWidth,1);
            fseek(fid,1*skippedPixels,'bof'); % read .bin File to roiImage
            z = fread(fid,pixelsPerFrame,'*uint8','b');
            img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
            roiImage(:,:,1) = flip(imrotate(img,-90),2);
            roiImage = uint8(roiImage); % convert double floating point data to unsignned 8bit integers
            workingImg{cc,1} = imcomplement(roiImage); % grab frame from image stack
        end
        desiredFig = figure;
        if length(workingImg) > 25
            workingImg = workingImg(1:25,1);
        end
        for dd = 1:length(workingImg)
            subplot(5,5,dd)
            imagesc(workingImg{dd,1})
            axis image
            axis off
            colormap gray
            title(['Session ' num2str(dd)])
        end
        % choose desired file
        drawnow()
        desiredFile = input('Which file looks best for ROI drawing: '); disp(' ')
        PupilData.firstFileOfDay{1,bb} = dayFilenames.(strDay)(desiredFile,:);
        % resize image if larger than 200x200
        if imageHeight > 200
            newROI = [1,1,199,199];
            newROIFig = figure;
            imagesc(workingImg{desiredFile,1})
            axis image
            axis off
            colormap gray
            h = imrect(gca,newROI); %#ok<IMRECT>
            addNewPositionCallback(h,@(p) title(mat2str(p,3)));
            fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
            setPositionConstraintFcn(h,fcn)
            position = wait(h);
            newImage = imcrop(workingImg{desiredFile,1},position);
            newImageFig = figure;
            imagesc(newImage)
            axis image
            axis off
            colormap gray
            PupilData.resizePosition{1,bb} = position;
            close(newROIFig)
            close(newImageFig)
        end
        clear workingImg
        close(desiredFig)
    end
    firstFileOfDay = PupilData.firstFileOfDay;
    save([trialData.animalID '_PupilData.mat'],'PupilData');
else
    firstFileOfDay = PupilData.firstFileOfDay;
end
% Create the desired window ROI for each day if it doesn't yet exist
for bb = 1:length(firstFileOfDay)
    firstFile = firstFileOfDay{1,bb};
    load(firstFile)
    [~,fileDate,fileID] = GetFileInfo_IOS(firstFile);
    strDay = ConvertDate_IOS(fileDate);
    if ~isfield(PupilData.EyeROI,(strDay))
        pupilCamFileID = [fileID '_PupilCam.bin'];
        fid = fopen(pupilCamFileID); % reads the binary file in to the work space
        fseek(fid,0,'eof'); % find the end of the video frame
        imageHeight = str2double(trialData.pupilCamPixelHeight);
        imageWidth = str2double(trialData.pupilCamPixelWidth);
        pixelsPerFrame = imageWidth*imageHeight;
        skippedPixels = pixelsPerFrame;
        roiImage = zeros(imageHeight,imageWidth,1);
        fseek(fid,1*skippedPixels,'bof'); % read .bin File to roiImage
        z = fread(fid,pixelsPerFrame,'*uint8','b');
        img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
        roiImage(:,:,1) = flip(imrotate(img,-90),2);
        roiImage = uint8(roiImage); % convert double floating point data to unsignned 8bit integers
        if isfield(PupilData,'resizePosition') == true
            roiImage = imcrop(roiImage,PupilData.resizePosition{1,bb});
        end
        workingImg = imcomplement(roiImage); % grab frame from image stack
        disp('Draw roi around eye'); disp(' ')
        eyeFigure = figure;
        title('Draw ROI around eye')
        [eyeROI] = roipoly(workingImg);
        % model the distribution of pixel intensities as a gaussian to estimate/isolate the population of pupil pixels
        pupilHistEdges = 1:1:256; % camera data is unsigned 8bit integers. Ignore 0 values
        threshSet = 1.5; % StD beyond mean intensity to binarize image for pupil tracking
        medFiltParams = [5,5]; % [x,y] dimensions for 2d median filter of images
        filtImg = medfilt2(workingImg,medFiltParams); % median filter image
        threshImg = double(filtImg).*eyeROI; % only look at pixel values in ROI
        [phat,~] = mle(reshape(threshImg(threshImg ~= 0),1,numel(threshImg(threshImg ~= 0))),'distribution','Normal');
        % figure for verifying pupil threshold
        pupilROIFig = figure;
        subplot(1,3,1)
        pupilHist = histogram(threshImg((threshImg ~= 0)),'BinEdges',pupilHistEdges,'Normalization','Probability');
        xlabel('Pixel intensities');
        ylabel('Bin Counts');
        title('Histogram of image pixel intensities')
        normCounts = pupilHist.BinCounts./sum(pupilHist.BinCounts); % normalizes bin count to total bin counts
        theFit = pdf('normal',pupilHist.BinEdges,phat(1),phat(2)); % generate distribution from mle fit of data
        normFit = theFit./sum(theFit); % normalize fit so sum of gaussian ==1
        intensityThresh = phat(1) + (threshSet*phat(2)); % set threshold as 4.5 sigma above population mean estimated from MLE
        testImg = threshImg;
        testImg(threshImg >= intensityThresh) = 1;
        testImg(threshImg < intensityThresh) = 0;
        testThresh = labeloverlay(roiImage(:,:,1),testImg);
        axis square
        subplot(1,3,2)
        plot(pupilHist.BinEdges(2:end),normCounts,'k','LineWidth',1);
        xlabel('Pixel intensities');
        ylabel('Normalized bin counts');
        title('Normalized histogram and MLE fit of histogram');
        hold on;
        plot(pupilHist.BinEdges,normFit,'r','LineWidth',2);
        xline(intensityThresh,'--c','LineWidth',1);
        legend({'Normalized Bin Counts','MLE fit of data','Pixel intensity threshold'},'Location','northwest');
        xlim([0,256]);
        axis square
        subplot(1,3,3)
        imshow(testThresh);
        title('Pixels above threshold');
        % check threshold
        threshOK = false;
        manualPupilThreshold = [];
        while threshOK == false
            disp(['Intensity threshold: ' num2str(intensityThresh)]); disp (' ')
            threshCheck = input('Is pupil threshold value ok? (y/n): ','s'); disp(' ')
            if strcmp(threshCheck,'y') == true
                threshOK = true;
            else
                intensityThresh = input('Manually set pupil intensity threshold: '); disp(' ')
                testImg(threshImg >= intensityThresh) = 1;
                testImg(threshImg < intensityThresh) = 0;
                testThresh = labeloverlay(roiImage(:,:,1),testImg);
                manualPupilThreshold = figure;
                imshow(testThresh);
                title('Pixels above threshold');
            end
        end
        % save the file to directory.
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Pupil ROI/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(eyeFigure,[dirpath fileID '_PupilROI'])
        savefig(pupilROIFig,[dirpath fileID '_PupilThreshold'])
        if isempty(manualPupilThreshold) == false
            savefig(manualPupilThreshold,[dirpath fileID '_ManualPupilThreshold'])
        end
        close all
        PupilData.EyeROI.(strDay) = eyeROI;
        PupilData.Threshold.(strDay) = intensityThresh;
        save([trialData.animalID '_PupilData.mat'],'PupilData');
    end
end