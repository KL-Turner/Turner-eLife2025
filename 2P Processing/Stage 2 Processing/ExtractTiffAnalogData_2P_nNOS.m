function [MScanData] = ExtractTiffAnalogData_2P_nNOS(MScanData, fileID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab and Yurong Gao
%________________________________________________________________________________________________________________________
%
%   Purpose: Uses the TIFF image stack and txt files from MScan to extract the vessel diameter and the analog data.
%________________________________________________________________________________________________________________________

% takes a tiff file (movie) and an analog ascii file, and extracts the diameters
MScan_analogData = [fileID '.TXT'];
disp(['Loading MScan file: ' MScan_analogData '...']); disp(' ');
analogData = load(MScan_analogData,'-ascii');
MScanData.data.corticalNeural = analogData(:,2);
MScanData.data.EMG = analogData(:,3);
MScanData.data.hippocampalNeural = analogData(:,4);
MScanData.data.forceSensor = analogData(:,5);
MScanData.notes.analogSamplingRate = 20000;
disp('Analyzing vessel projections from defined polygons...'); disp(' ');
[MScanData] = GetDiameterFromMovie_2P_nNOS(MScanData,fileID);
try
    [MScanData] = FWHM_MovieProjection_2P_nNOS(MScanData,[MScanData.notes.startframe,MScanData.notes.endframe]);
catch
    disp([MScanData.notes.imageID ' FWHM calculation failed!']); disp(' ')
    keyboard
end
try
    % 1 dural/vein, >40% changes spline, artery: >60% spline
    % 2 dural/vein, >30% changes interpolate, artery: >50% interpolate
    if strcmp(MScanData.notes.vesselType,'D') || strcmp(MScanData.notes.vesselType,'V')
        MScanData.data.vesselDiameter = RemoveMotion_2P_nNOS(MScanData.data.tempVesselDiameter,MScanData.notes.vesselROI.modalFixedDiameter,2,0.3);
    else
        MScanData.data.vesselDiameter = RemoveMotion_2P_nNOS(MScanData.data.tempVesselDiameter,MScanData.notes.vesselROI.modalFixedDiameter,2,0.5);
    end
    [diamPerc,S,f] = DiamPercPower_2P_nNOS(MScanData.data.vesselDiameter,MScanData.notes.vesselROI.modalFixedDiameter,MScanData.notes.frameRate);
    MScanData.notes.vessel.diamPerc = diamPerc;
    MScanData.notes.vessel.power_f = f;
    MScanData.notes.vessel.power_S = S;
catch
    disp([MScanData.notes.imageID ' Diameter percentage analysis failed!']); disp(' ')
    keyboard
end

end
