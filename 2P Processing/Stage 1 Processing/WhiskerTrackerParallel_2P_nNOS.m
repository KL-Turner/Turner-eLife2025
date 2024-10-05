function [angle] = WhiskerTrackerParallel_2P(fileName)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% Variable Setup
theta = -40:80;   % Angles used for radon
% Import whisker movie
importStart = tic;
baslerFrames = ReadBinFileU8MatrixGradient_2P([fileName '_WhiskerCam.bin'],350,30);
importTime = toc(importStart);
disp(['WhiskerTrackerParallel: Binary file import time was ' num2str(importTime) ' seconds.']); disp(' ')
% Transfer the images to the GPU
gpuTrans1 = tic;
gpuFrame = gpuArray(baslerFrames);
gpuTransfer = toc(gpuTrans1);
disp(['WhiskerTrackerParallel: GPU transfer time was ' num2str(gpuTransfer) ' seconds.']); disp(' ')
% PreAllocate array of whisker angles, use NaN as a place holder
angle = NaN*ones(1,length(baslerFrames));
radonTime1 = tic;
for a = 1:(length(baslerFrames) - 1)
    % Radon on individual frame
    [R,~] = radon(gpuFrame(:,:,a),theta);
    % Get transformed image from GPU and calculate the variance
    colVar = var(gather(R));
    % Sort the columns according to variance
    ordVar = sort(colVar);
    % Choose the top 0.1*number columns which show the highest variance
    thresh = round(numel(ordVar)*0.9);
    sieve = gt(colVar,ordVar(thresh));
    % Associate the columns with the corresponding whisker angle
    angles = nonzeros(theta.*sieve);
    % Calculate the average of the whisker angles
    angle(a) = mean(angles);
end
radonTime = toc(radonTime1);
disp(['WhiskerTrackerParallel: Whisker Tracking time was ' num2str(radonTime) ' seconds.']); disp(' ')
inds = isnan(angle) == 1;
angle(inds) = [];

end
