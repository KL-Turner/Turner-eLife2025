function [] = Analyze2PDataNotes_2P(msExcelFile)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Extract the vessel notes and information from a MS Excel file, and create a struct titled '_MScanData.mat'.
%________________________________________________________________________________________________________________________

% Read the image info from the formated xls file
[~,~,alldata] = xlsread(msExcelFile);
for a = 2:size(alldata,1)   % Loop through all rows of the excel sheet except the first row
    clear MScanData
    %% notes
    MScanData.notes.date = num2str(alldata{a,1});
    MScanData.notes.animalID = alldata{a,2};
    MScanData.notes.imageID = alldata{a,3};
    MScanData.notes.movieType = alldata{a,4};
    MScanData.notes.laserPower = alldata{a,5};
    MScanData.notes.objectiveID = alldata{a,6};
    MScanData.notes.frameRate = alldata{a,7};
    MScanData.notes.numberOfFrames = alldata{a,8};
    MScanData.notes.vesselType = alldata{a,9};
    MScanData.notes.vesselDepth = alldata{a,10};
    MScanData.notes.comments = alldata{a,11};
    MScanData.notes.vesselID = alldata{a,12};
    MScanData.notes.drug = alldata{a,13};
    currentFileID = ([MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_' MScanData.notes.vesselID '_MScanData.mat']);
    if ~exist(currentFileID,'file')   % Only run analysis if the current file doesn't exist yet
        % Vessel diameter calculation    for surface vessels
        if strcmp(MScanData.notes.movieType,'MS') == true
            [MScanData] = DiamCalcSurfaceVessel_2P(MScanData,[MScanData.notes.date '_' MScanData.notes.imageID]);
            % Vessel diameter (area) calculation for penetrating vessels
        elseif strcmp(MScanData.notes.movieType,'MP') == true
            [MScanData] = AreaCalcPenVessel_2P(MScanData,[MScanData.notes.date '_' MScanData.notes.imageID]);
            % Line scan calculation for capillaries
        elseif strcmp(MScanData.notes.movieType,'MC') == true
            [MScanData] = CapLinesScan_2P(MScanData,[MScanData.notes.date '_' MScanData.notes.imageID]);
        end
        % Checklist for analysis steps - debugging purposes
        MScanData.notes.checklist.analyzeVessel = false;
        MScanData.notes.checklist.processData = false;
        MScanData.notes.checklist.offsetCorrect = false;
        % Save the RawData file for the current movie type
        disp(['File Created. Saving MScanData File ' num2str(a - 1) '...']); disp(' ')
        save([MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_' MScanData.notes.vesselID '_MScanData'],'MScanData')
        close all
    else
        disp([currentFileID ' already exists in the current directory. Continuing...']); disp(' ')
    end
end

end
