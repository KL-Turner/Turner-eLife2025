function [PupilStruct] = TrackPupilDiameter_IOS_eLife2025(fileName,iteration,trialData,imagingParameters)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% check pupil frames and establish first file of each day
if strcmp(imagingParameters.pupilOptions,'Yes') == true
    if iteration == 1
        [PupilData] = SelectPupilFiles_IOS_eLife2025(trialData);
    else
        load([trialData.animalID '_PupilData.mat'])
    end
    disp(['Running pupil tracker for: ' fileName]); disp(' ')
    pupilStart = tic;
    [~,fileDate,~] = GetFileInfo_IOS_eLife2025(fileName);
    strDay = ConvertDate_IOS_eLife2025(fileDate);
    idx = 0;
    for qq = 1:size(PupilData.firstFileOfDay,2)
        if strfind(PupilData.firstFileOfDay{1,qq},fileDate) >= 5
            idx = qq;
        end
    end
    % run pupil/blink tracking on all data files
    theAngles = 1:1:180; % projection angles measured during radon transform of pupil
    radonThresh = 0.05; % arbitrary threshold used to clean up radon transform above values == 1 below == 0
    pupilThresh = 0.25; % arbitrary threshold used to clean up inverse radon transform above values == 1 below == 0
    blinkThresh = 0.35; % arbitrary threshold used to binarize data for blink detection above values == 1 below == 0
    medFiltParams = [5,5]; % [x,y] dimensions for 2d median filter of images
    fid = fopen(fileName); % reads the binary file in to the work space
    fseek(fid,0,'eof'); % find the end of the video frame
    fileSize = ftell(fid); % calculate file size
    fseek(fid,0,'bof'); % find the begining of video frames
    imageHeight = str2double(trialData.pupilCamPixelHeight); % how many pixels tall is the frame
    imageWidth = str2double(trialData.pupilCamPixelWidth); % how many pixels wide is the frame
    pixelsPerFrame = imageWidth*imageHeight;
    skippedPixels = pixelsPerFrame;
    nFramesToRead = floor(fileSize/(pixelsPerFrame));
    imageStack = zeros(imageHeight,imageWidth,nFramesToRead); % empty variable for entire pupil tracking video
    % read .bin file to imageStack
    for dd = 1:nFramesToRead
        fseek(fid,(dd - 1)*skippedPixels,'bof');
        z = fread(fid,pixelsPerFrame,'*uint8','b'); % read in next movie frame
        img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight); % reshape to correct image dimension
        if isfield(PupilData,'resizePosition') == true
            imageStack(:,:,dd) = imcrop(flip(imrotate(img,-90),2),PupilData.resizePosition{1,idx}); % orient movie frame in proper dorsal-ventral view
        else
            imageStack(:,:,dd) = flip(imrotate(img,-90),2);
        end
    end
    % convert double floating point data to unsignned 8bit integers
    imageStack = uint8(imageStack); % change to 8-bit to conserve memory
    % grab frame from image stack
    workingImg = imcomplement(imageStack(:,:,2)); % invert image pixel intensities so pupil becomes "bright" and sclera "dark"
    PupilStruct.firstFrame = imageStack(:,:,2);
    % pre-allocate empty structures
    pupilArea(1:size(imageStack,3)) = NaN; % area of pupil
    pupilMajor(1:size(imageStack,3)) = NaN; % length of major axis of pupil
    pupilMinor(1:size(imageStack,3)) = NaN; % length of minor axis of pupil
    pupilCentroid(1:size(imageStack,3),2) = NaN; % center of pupil
    pupilBoundary(1:size(imageStack,1),1:size(imageStack,2),1:size(imageStack,3)) = NaN;
    % select ROI containing pupil
    eyeROI = PupilData.EyeROI.(strDay);
    % select pupil intensity threshold
    intensityThresh = PupilData.Threshold.(strDay); % this is the threshold used to binarize inverted movie frame
    imageFrames = gpuArray(imageStack);
    roiInt(1:size(imageFrames,3)) = NaN;
    roiInt = gpuArray(roiInt);
    correctedFlag = false;
    lastFrameOK = 1;
    for frameNum = 1:size(imageStack,3)
        filtImg = medfilt2(imcomplement(imageFrames(:,:,frameNum)),medFiltParams); % invert and median filter raw movie frame
        % only look at pixel values in ROI
        threshImg = uint8(double(filtImg).*eyeROI); % apply mask to remove pixel data outside of eye
        roiIntTemp = sum(threshImg,1);
        roiInt(frameNum) = sum(roiIntTemp,2); % sum pixel intensities inside eye used for blink detection
        isoPupil = threshImg;
        isoPupil(isoPupil < intensityThresh) = 0; % set pixel values below threshold to 0
        isoPupil(isoPupil >= intensityThresh) = 1; % set pixel values above threshold to 1
        radPupil = radon(isoPupil); % transform image to radon space
        minPupil = min(radPupil,[],1);
        minMat = repmat(minPupil,size(radPupil,1),1);
        maxMat = repmat(max((radPupil - minMat),[],1),size(radPupil,1),1);
        % normalize each projection angle to its min and max values. Each value should now be between [0 1]
        normPupil = (radPupil - minMat)./maxMat; % normalize radon transform to account for inter-trial luminance changes
        threshPupil = normPupil;
        % binarize radon projection
        threshPupil(normPupil >= radonThresh) = 1; % threshold in normalized radon space
        threshPupil(normPupil < radonThresh) = 0;
        % transform back to image space
        radonPupil = gather(iradon(double(threshPupil),theAngles,'linear','Hamming',size(workingImg,2))); % transform back to image space.
        % find area corresponding to pupil on binary image
        [~,pupilBoundaries,objNum,~] = bwboundaries(radonPupil >pupilThresh*max(radonPupil(:)),8,'noholes'); % identify objects within image
        inspectFlag = false; % flag to indicate if first frame has multiple ROI defined
        if objNum == 1
            fillPupil = pupilBoundaries;
            % fill any subthreshold pixels inside the pupil boundary
            fillPupil = imfill(fillPupil,8,'holes'); % fill any holes within identified objects
            % get properties of detected objects in image
            areaFilled = regionprops(fillPupil,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength'); %identify binarized image object properties. This should be the Pupil!
        else
            if frameNum == 1
                fillPupil = pupilBoundaries;
                % fill any subthreshold pixels inside the pupil boundary
                fillPupil = imfill(fillPupil,8,'holes'); % fill any holes within identified objects
                areaFilled = regionprops(fillPupil,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength');%identify binarized image object properties. This should be the Pupil!
                for num = 1:size(areaFilled,1)
                    theArea(num) = areaFilled(num).FilledArea; %#ok<*SAGROW>
                end
                if isempty(areaFilled) == false
                    maxArea = max(theArea); % find the ROI with the largest area
                    areaLogical = theArea == maxArea;
                    areaFilled = areaFilled(areaLogical); % keep ROI with largest area, assumed to be pupil
                    areaFilled = areaFilled(1); % if two ROI have same area take first one, this will need to be corrected
                    inspectFlag = true; % use this to flag files that have multiple ROI for first frame
                else
                    areaFilled = [];
                    areaFilled.Centroid = NaN;
                    areaFilled.MajorAxisLength = NaN;
                    areaFilled.MinorAxisLength = NaN;
                    areaFilled.FilledArea = NaN;
                    inspectFlag = true; % use this to flag files that have multiple ROI for first frame
                end
            else
                areaFilled = regionprops(pupilBoundaries,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength'); % identify binarized image object properties.
            end
        end
        if frameNum > 1
            if abs((roiInt(frameNum) - roiInt(frameNum - 1))/roiInt(frameNum)) >= blinkThresh % Exclude fitting frame during blinks
                areaFilled = [];
                fillPupil(:) = 0;
            end
        end
        if ~isempty(areaFilled) % is an pupil identified
            if size(areaFilled,1) > 1 % is the pupil fragmented in to multiple ROI
                clear theArea areaLogical
                for num = 1:size(areaFilled,1)
                    theArea(num) = areaFilled(num).FilledArea;
                end
                maxArea = max(theArea); % find the ROI with the largest area
                areaLogical = theArea == maxArea;
                areaFilled = areaFilled(areaLogical);
                areaFilled = areaFilled(1);
                % check for aberrant pupil diameter changes
                if frameNum > 1
                    fracChange = (maxArea - pupilArea(frameNum-1))/pupilArea(frameNum-1); % frame-wise fractional change
                    volFlag = fracChange <- 0.1; % does the change exceed a 10% reduction in pupil size
                    if ~isnan(pupilArea(frameNum - 1)) % does the current frame follow a blink
                        if volFlag == true
                            if correctedFlag == false
                                lastFrameOK = frameNum - 1; % find last good pupil frame for size and location comparison
                                correctedFlag = true;
                            end
                            % correct aberrant diameters by altering radon threshold
                            pupilSweep = intensityThresh - (1:100); % adjust threshold of binary image instead of radon image KWG
                            for sweepNum = 1:size(pupilSweep,2)
                                if volFlag == true
                                    isoSweep = threshImg; % get video frame
                                    isoSweep(isoSweep < pupilSweep(sweepNum)) = 0; % set all pixel intensities below new thresh to zero
                                    isoSweep(isoSweep >= pupilSweep(sweepNum)) = 1; % set all pixel intensities above new thresh to one
                                    inPupil = [];
                                    [~,correctBoundaries,~] = bwboundaries(gather(isoSweep),8,'noholes'); % Identify above threshold regions
                                    fillCorrection = correctBoundaries;
                                    fillCorrection =imfill(fillCorrection,8,'holes'); % fill holes in regions
                                    areaCorrect = regionprops(fillCorrection,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength','PixelList'); %find region properties
                                    for areaNum = 1:size(areaCorrect,1)
                                        areaCentroid = round(areaCorrect(areaNum).Centroid,0);
                                        inPupil(areaNum) = pupilBoundary(areaCentroid(2),areaCentroid(1),lastFrameOK); % determine if centroid of region was in previous pupil volume
                                    end
                                    theInds = find(inPupil == 1);
                                    keepRegions = areaCorrect(theInds); %#ok<FNDSB> % keep only regions within previous pupil measurements
                                    keepMask = zeros(size(fillCorrection,1),size(fillCorrection,2));
                                    for keepNum = 1:size(keepRegions,1)
                                        for pixNum = 1:size(keepRegions(keepNum).PixelList)
                                            keepMask(keepRegions(keepNum).PixelList(pixNum,2),keepRegions(keepNum).PixelList(pixNum,1)) = 1; % remap kept regions in to image frame
                                        end
                                    end
                                    fuseMask = bwconvhull(keepMask); % use convex hull operation to enclose regions previously within pupil
                                    fusedCorrect = regionprops(fuseMask,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength','PixelList'); % measure new corrected pupil volume properties
                                    if length(fusedCorrect) > 1
                                        fusedCorrect = fusedCorrect(1);
                                    end
                                    if ~isempty(fusedCorrect)
                                        fracChange = (fusedCorrect.FilledArea - pupilArea(lastFrameOK))/pupilArea(lastFrameOK); %compare pupil size to last known pupil size
                                        volFlag = fracChange < -0.1;
                                    else
                                        fracChange = 100; % if no data present force frac change outside threshold to keep frame
                                    end
                                end
                            end
                            % this can be used to insert NaN if the change is > 10%
                            if  abs(fracChange) < 0.1 % changed to only fill data withing a +/- 10% change in area
                                fillPupil = fuseMask;
                                areaFilled = fusedCorrect;
                            end
                            if ~exist('correctedFrames','var')
                                frameInd = 1;
                                correctedFrames(frameInd) = frameNum;
                            else
                                frameInd = frameInd + 1;
                                correctedFrames(frameInd) = frameNum;
                            end
                        else
                            correctedFlag = false;
                            lastFrameOK = frameNum;
                        end
                    end
                end
            else
                if frameNum > 1
                    % check for aberrant pupil diameter changes
                    fracChange = (areaFilled.FilledArea - pupilArea(frameNum - 1))/pupilArea(frameNum - 1);
                    volFlag = fracChange < -0.1;
                    if ~isnan(pupilArea(frameNum - 1))
                        if volFlag == true
                            if correctedFlag == false
                                lastFrameOK = frameNum - 1;
                                correctedFlag = true;
                            end
                            % correct aberrant diameters with previous pupil locations
                            pupilSweep = intensityThresh - (1:100); % adjust threshold of binary image instead of radon image KWG
                            for sweepNum = 1:size(pupilSweep,2)
                                if volFlag == true
                                    isoSweep = threshImg; % get video frame
                                    isoSweep(isoSweep < pupilSweep(sweepNum)) = 0; % set all pixel intensities below thresh to zero
                                    isoSweep(isoSweep >= pupilSweep(sweepNum)) = 1; % set all pixel intensities above thresh to one
                                    radSweep = radon(isoSweep); % take radon transform of binarized image
                                    minPupil = min(radSweep,[],1);
                                    minMat = repmat(minPupil,size(radSweep,1),1);
                                    maxMat = repmat(max((radSweep - minMat),[],1),size(radSweep,1),1);
                                    % normalize each projection angle to its min and max values. Each value should now be between [0 1]
                                    normSweep = (radSweep - minMat)./maxMat;
                                    threshSweep = normSweep;
                                    % binarize radon projection
                                    threshSweep(normSweep >= radonThresh) = 1;
                                    threshSweep(normSweep < radonThresh) = 0;
                                    % transform back to image space
                                    radonSweep = gather(iradon(double(threshSweep),theAngles,'linear','Hamming',size(workingImg,2)));
                                    sweepArea = [];
                                    % get image objects
                                    [~,sweepBoundaries] = bwboundaries(radonSweep > pupilSweep(sweepNum)*max(radonSweep(:)),8,'noholes');
                                    fillSweep = sweepBoundaries;
                                    fillSweep = imfill(fillSweep,8,'holes');
                                    % get object properties
                                    areaSweep = regionprops(fillSweep,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength');
                                    for num = 1:size(areaSweep,1)
                                        sweepArea(num) = areaSweep(num).FilledArea; %#ok<*AGROW>
                                    end
                                    maxSweep = max(sweepArea);
                                    sweepLogical = sweepArea == maxSweep;
                                    fracChange = (maxSweep - pupilArea(lastFrameOK))/pupilArea(lastFrameOK);
                                    volFlag = fracChange < -0.1;
                                end
                            end
                            % this can be used to insert NaN if the change is > 10%
                            if abs(fracChange) < 0.1
                                fillPupil = fillCorrection;
                                areaFilled = areaCorrect(sweepLogical);
                            end
                            if ~exist('correctedFrames','var')
                                frameInd = 1;
                                correctedFrames(frameInd)=frameNum;
                            else
                                frameInd = frameInd + 1;
                                correctedFrames(frameInd) = frameNum;
                            end
                        else
                            correctedFlag = false;
                            lastFrameOK = frameNum;
                        end
                    end
                end
            end
            pupilArea(frameNum) = areaFilled.FilledArea;
            pupilMajor(frameNum) = areaFilled.MajorAxisLength;
            pupilMinor(frameNum) = areaFilled.MinorAxisLength;
            pupilCentroid(frameNum,:) = areaFilled.Centroid;
            pupilBoundary(:,:,frameNum) = fillPupil;
            holdMat = labeloverlay(imageStack(:,:,frameNum),fillPupil,'Transparency',0.8); % this is the measured pupil overlayed on the movie frame
            if size(holdMat,3) == 1
                overlay(:,:,:,frameNum) = repmat(holdMat,1,1,3);
            else
                overlay(:,:,:,frameNum) = holdMat; % this is the measured pupil overlayed on the movie
            end
        else
            pupilArea(frameNum) = NaN;
            pupilMajor(frameNum) = NaN;
            pupilMinor(frameNum) = NaN;
            pupilCentroid(frameNum,:) = NaN;
            pupilBoundary(:,:,frameNum) = fillPupil;
            holdMat = labeloverlay(imageStack(:,:,frameNum),fillPupil); % this is the measured pupil overlayed on the movie frame
            if size(holdMat,3) == 1
                overlay(:,:,:,frameNum) = repmat(holdMat,1,1,3);
            else
                overlay(:,:,:,frameNum) = holdMat; % this is the measured pupil overlayed on the movie
            end
        end
    end
    % save data
    PupilStruct.pupilArea = pupilArea;
    PupilStruct.pupilMajor = pupilMajor;
    PupilStruct.pupilMinor = pupilMinor;
    PupilStruct.pupilCentroid = pupilCentroid;
    PupilStruct.eyeROI = eyeROI;
    PupilStruct.roiIntensity = gather(roiInt);
    PupilStruct.inspectFile = inspectFlag; % did the first frame include more than one ROI?
    blinks = find((abs(diff(PupilStruct.roiIntensity))./PupilStruct.roiIntensity(2:end)) >= blinkThresh) + 1;
    PupilStruct.blinkFrames = overlay(:,:,:,blinks);
    PupilStruct.blinkInds = blinks;
    pupilStop = toc(pupilStart);
    disp(['Whisker Tracking time was ' num2str(pupilStop) ' seconds.']); disp(' ')
else
    % save data
    PupilStruct.firstFrame = ones(200,200);
    PupilStruct.pupilArea = ones(1,str2double(trialData.pupilCamSamplingRate)*str2double(trialData.trialDuration_sec));
    PupilStruct.pupilMajor = ones(1,str2double(trialData.pupilCamSamplingRate)*str2double(trialData.trialDuration_sec));
    PupilStruct.pupilMinor = ones(1,str2double(trialData.pupilCamSamplingRate)*str2double(trialData.trialDuration_sec));
    PupilStruct.pupilCentroid = ones(str2double(trialData.pupilCamSamplingRate)*str2double(trialData.trialDuration_sec),2);
    PupilStruct.eyeROI = ones(200,200);
    PupilStruct.roiIntensity = ones(1,str2double(trialData.pupilCamSamplingRate)*str2double(trialData.trialDuration_sec));
    PupilStruct.inspectFile = 0;
    PupilStruct.blinkFrames = [];
    PupilStruct.blinkInds = [];
end