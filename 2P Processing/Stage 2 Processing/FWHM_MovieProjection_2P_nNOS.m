function [MScanData] = FWHM_MovieProjection_2P(MScanData,theFrames)
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

for f = min(theFrames):max(theFrames)
    % Add in a 5 pixel median filter
    MScanData.data.rawVesselDiameter(f) = CalcFWHM_2P(medfilt1(MScanData.notes.vesselROI.projection(f,:),5));
end
MScanData.data.tempVesselDiameter = MScanData.data.rawVesselDiameter*MScanData.notes.xFactor;
[holdHist,d] = hist(MScanData.data.tempVesselDiameter,0:.25:100); %#ok<HIST>
[~,maxD] = max(holdHist);
MScanData.notes.vesselROI.modalFixedDiameter = d(maxD);

end