%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Pull vessel notes from Excel sheet.
%            2) Analyze vessel diameter.
%            3) Add the LabVIEW filename to each MScan file.
%            4) Analyze neural bands, force sensors, and whisker motion.
%            5) Correct LabVIEW time offset.
%            6) Combine LabVIEWData and MScan Data files to ceate MergedData.
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
zap;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
msExcelFile = uigetfile('*.xlsx');
%% BLOCK PURPOSE: [1] Use ms Excel sheet to create MScanData.mat files with vessel information.
disp('Analyzing Block [1] Pulling vessel notes from Excel sheet.'); disp(' ')
Analyze2PDataNotes_2P_nNOS(msExcelFile);
%% BLOCK PURPOSE: [2] Analyze vessel diameter and add it to MScanData.mat.
disp('Analyzing Block [2] Analyzing vessel diameter.'); disp(' ')
mscanDirectory = dir('*_MScanData.mat');
mscanDataFiles = {mscanDirectory.name}';
mscanDataFiles = char(mscanDataFiles);
Analyze2PDiameter_2P_nNOS(mscanDataFiles);
%% BLOCK PURPOSE: [3] Add LabVIEW filename to MScan file
disp('Analyzing Block [3] Adding LabVIEW file IDs to MScan data files.'); disp(' ')
labviewDirectory = dir('*_LabVIEWData.mat');
labviewDataFiles = {labviewDirectory.name}';
labviewDataFiles = char(labviewDataFiles);
AddLabVIEWFileID_2P_nNOS(mscanDataFiles,labviewDataFiles,msExcelFile);
%% BLOCK PURPOSE: [4] Process neural, whiskers, and force sensor data.
disp('Analyzing Block [4] Analyzing neural bands, force sensors, and whiskers.'); disp(' ')
Process2PDataFiles_2P_nNOS(labviewDataFiles,mscanDataFiles)
%% BLOCK PURPOSE: [5] Correct the offset between the MScan and LabVIEW acquisiton.
disp('Analyzing Block [5] Correcting LabVIEW time offset.'); disp(' ')
trimTime = 15;   % sec
CorrectLabVIEWOffset_2P_nNOS(mscanDataFiles,trimTime)
%% BLOCK PURPOSE: [6] Combine the MScan and LabVIEW structures into one.
disp('Analyzing Block [6] Combing LabVIEWData and MScan Data files to create MergedData.'); disp(' ')
CombineLabVIEWMScanFiles_2P_nNOS(mscanDataFiles)

disp('Two Photon Stage Two Processing - Complete.'); disp(' ')
