function [MScanData] = GetVelocityLineScan_2P_eLife2025(MScanData,fileID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

% takes a tiff file (movie) and an analog ascii file, and extracts the line scan velocity
MScan_analogData = [fileID '.TXT'];
disp(['Loading MScan file: ' MScan_analogData '...']); disp(' ');
analogData = load(MScan_analogData,'-ascii');
MScanData.data.corticalNeural = analogData(:,2);
MScanData.data.EMG = analogData(:,3);
MScanData.data.hippocampalNeural = analogData(:,4);
MScanData.data.forceSensor = analogData(:,5);
MScanData.notes.analogSamplingRate = 20000;
maxVelocity = 10000; % maximum physiological velocity in micrometers/sec
angleSpan = 15; % number of degrees (+/-) around previous angle to look
angleResolution = .1; % how accurate to determine the angle
interleaved = 1;
disp(['Processing ' MScanData.notes.date '_' MScanData.notes.imageID '...']); disp(' ')
currentFilePath = mfilename('fullpath');
MScanData.notes.code = GetFunctionCode_2P_eLife2025(currentFilePath);
[MScanData] = GetVelocityRadon_2P_eLife2025(fileID,angleSpan,angleResolution,maxVelocity,interleaved,MScanData);
% calculate the velocity using the cross correlation method of Kim et al, 2012
[MScanData] = LineScanXCVelocity_2P_eLife2025(MScanData);

end
