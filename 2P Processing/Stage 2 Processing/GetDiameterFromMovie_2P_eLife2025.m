function [MScanData] = GetDiameterFromMovie_2P_nNOS(MScanData,fileID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab and Yurong Gao
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________

MScanData.notes.firstFrame = imread(fileID,'TIFF','Index',1);
fftFirstFrame = fft2(double(MScanData.notes.firstFrame));
X = repmat(1:MScanData.notes.xSize,MScanData.notes.ySize,1);
Y = repmat((1:MScanData.notes.ySize)',1,MScanData.notes.xSize);
MScanData.notes.vesselROI.projectionAngle = atand(diff(MScanData.notes.vesselROI.vesselLine.position.xy(:,1))/diff(MScanData.notes.vesselROI.vesselLine.position.xy(:,2)));
for theFrame = MScanData.notes.startframe:MScanData.notes.endframe
    rawFrame = imread(fileID,'TIFF','Index',theFrame);
    fftRawFrame = fft2(double(rawFrame));
    [MScanData.notes.pixelShift(:,theFrame),~] = DftRegistration_2P_nNOS(fftFirstFrame,fftRawFrame,1);
    inpolyFrame = inpolygon(X + MScanData.notes.pixelShift(3,theFrame),Y + MScanData.notes.pixelShift(4,theFrame),MScanData.notes.vesselROI.boxPosition.xy(:,1),MScanData.notes.vesselROI.boxPosition.xy(:,2));
    boundedrawFrame = rawFrame.*uint16(inpolyFrame);
    MScanData.notes.vesselROI.projection(theFrame,:) = radon(boundedrawFrame,MScanData.notes.vesselROI.projectionAngle);
end

end
