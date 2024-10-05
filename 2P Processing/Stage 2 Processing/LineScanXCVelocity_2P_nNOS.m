function [MScanData,xcFig] = LineScanXCVelocity_2P_nNOS(MScanData)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: uses the cross correlation method to measure velocity using the method of Kim et al, PLoS ONE 2012
%________________________________________________________________________________________________________________________

% maxspot-amplitude of the cross correlation peak.
% xcorr_image- power-spectra normalized cross correlation image
% xshift- shift of the peak away from the center, in pixels
xSpread = round(max(1,.5/MScanData.notes.xFactor)); % spatial gaussian with 0.5um std
tSpread = round(MScanData.notes.tFactor/10); % temporal gaussian with 10ms std
theKernel = gaussian2d_2P_nNOS(xSpread,tSpread,3*xSpread,3*tSpread);
theImage = double(MScanData.data.bloodFlow.Image'); % DOUBLE CHECK IS THIS REALLY THE RIGHT IMAGE????
nLines = size(theImage,2);
nPoints = size(theImage,1);
averageLine = mean(theImage,2);
xcImage = zeros(size(theImage));
for tt = 1:nLines
    theImage(:,tt) = theImage(:,tt) - averageLine; % subtract out the background
end
for tt = 3:nLines
    %take the convolution of the two line-scans normalized by the joint
    xcImage(:,tt) = ifft(fft(theImage(:,tt)).*conj(fft(theImage(:,tt - 1)))./sqrt(abs(fft(theImage(:,tt))).*abs(fft(theImage(:,tt - 1)))));
end
% convolve with a gaussian space-time kernel to average velocity
xcImage = conv2(fftshift(xcImage,1),theKernel,'same');% fftshift puts the lower absolute values for the velocities togetehr in the middel of the matrix,then we filter with a matched kernel
[~,maxSpot] = max(xcImage(:,:)); % find the peak in the cross correlation
xShift = (round(nPoints/2) - maxSpot);
[b,a] = butter(2,200/(2*MScanData.notes.tFactor)); % 200 Hz cutoff
xcVelocity = filtfilt(b,a,xShift*MScanData.notes.xFactor*MScanData.notes.tFactor); % calculate the velocity from the displacement and filter
MScanData.data.bloodFlow.xcVelocity = xcVelocity; % velocity obtained with the cross correlation method
params.Fs = MScanData.notes.tFactor;
params.tapers = [20,39];
[S_xc,f_xc] = mtspectrumc(xShift - mean(xShift(:)),params);
% figures
xcFig = figure;
sgtitle([MScanData.notes.animalID ' ' MScanData.notes.date ' ' MScanData.notes.imageID])
subplot(1,2,1)
hold off
imagesc((1:length(xcImage))/MScanData.notes.tFactor,(size(xcImage,1)/2:-1:-size(xcImage,1)/2)*MScanData.notes.tFactor*MScanData.notes.xFactor/1000,xcImage)
hold on
plot((1:length(MScanData.data.bloodFlow.xcVelocity))/MScanData.notes.tFactor,MScanData.data.bloodFlow.xcVelocity/1000,'w')
xlabel('Time (s)')
ylabel('velocity, mm/sec')
title('XC Velocity');
axis xy
axis square
set(gca,'box','off')
subplot(1,2,2)
loglog(f_xc,S_xc)
hold on
x1 = xline(10,'r');
title('Power Spectrum')
xlabel('Freq (Hz)')
ylabel('Power (a.u.)')
legend(x1,'Heart rate peak')
axis tight
axis square
set(gca,'box','off')
hold on
% save path
[pathstr,~,~] = fileparts(cd);
dirpath = [pathstr '/Figures/XC Velocity/'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(xcFig,[dirpath MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_XC_Velocity']);
close(xcFig)

end
