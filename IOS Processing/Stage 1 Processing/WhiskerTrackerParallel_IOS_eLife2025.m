function [angle] = WhiskerTrackerParallel_IOS_eLife2025(fileName,trialData)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
disp(['Running whisker tracker for: ' fileName]); disp(' ')
% angles used for radon
theta = -40:80;
% import whisker movie
imageHeight = str2double(trialData.whiskCamPixelHeight);
imageWidth = str2double(trialData.whiskCamPixelWidth);
% calculate pixels per frame for fread
pixelsPerFrame = imageWidth*imageHeight;
% open the file, get file size, back to the begining
fid = fopen(fileName);
fseek(fid,0,'eof');
fileSize = ftell(fid);
fseek(fid,0,'bof');
% identify the number of frames to read. Each frame has a previously defined width and height (as inputs), U8 has a depth of 1.
nFrameToRead = floor(fileSize/(pixelsPerFrame));
% pre-allocate
imageGrad = int8(zeros(imageWidth,imageHeight,nFrameToRead));
for n = 1:nFrameToRead
    z = fread(fid,pixelsPerFrame,'*uint8',0,'l');
    indImg = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
    imageGrad(:,:,n) = int8(gradient(double(indImg)));
end
fclose(fid);
% transfer the images to the GPU
gpuFrame = gpuArray(imageGrad);
% pre-allocate array of whisker angles, use NaN as a place holder
angle = NaN*ones(1,length(imageGrad));
radonTime1 = tic;
for aa = 1:(length(imageGrad) - 1)
    % radon on individual frame
    [R,~] = radon(gpuFrame(:,:,aa),theta);
    % get transformed image from GPU and calculate the variance
    colVar = var(gather(R));
    % sort the columns according to variance
    ordVar = sort(colVar);
    % choose the top 0.1*number columns which show the highest variance
    thresh = round(numel(ordVar)*0.9);
    sieve = gt(colVar,ordVar(thresh));
    % associate the columns with the corresponding whisker angle
    angles = nonzeros(theta.*sieve);
    % calculate the average of the whisker angles
    angle(aa) = mean(angles);
end
inds = isnan(angle);
angle(inds) = [];
radonTime = toc(radonTime1);
disp(['Whisker Tracking time was ' num2str(radonTime) ' seconds.']); disp(' ')