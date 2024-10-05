function [] = CombineLabVIEWMScanFiles_2P(mscanDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Combine the MScan and LabVIEW data structures into one.
%________________________________________________________________________________________________________________________

% Loop through each set of files and combine the processed and corrected LabVEIW and MScan data into one
for a = 1:size(mscanDataFiles,1)
    disp(['Combining the data from LabVIEW and MScan file(s) number ' num2str(a) ' of ' num2str(size(mscanDataFiles, 1)) '...']); disp(' ');
    mscanDataFile = mscanDataFiles(a,:);
    load(mscanDataFile);
    labviewDataFile = MScanData.notes.labviewFileID;
    load(labviewDataFile);  
    [animalID,hem,~,fileID] = GetFileInfo_2P(labviewDataFile);
    vesselID = MScanData.notes.vesselID;
    imageID = MScanData.notes.imageID;   
    % Pull the notes and data from LabVIEW
    MergedData.notes.LabVIEW = LabVIEWData.notes;
    MergedData.data.whiskerAngle = LabVIEWData.data.dsWhiskerAngle_trim;
    MergedData.data.rawWhiskerAngle = LabVIEWData.data.whiskerAngle_trim;
    MergedData.data.binWhiskerAngle = LabVIEWData.data.binWhiskerAngle_trim;
    MergedData.data.forceSensorL = LabVIEWData.data.dsForceSensorL_trim;
    MergedData.data.binForceSensorL = LabVIEWData.data.binForceSensorL_trim;
    % Save solenoid times (in seconds). Identify the solenoids by amplitude.
    MergedData.data.solenoids.LPadSol = find(diff(LabVIEWData.data.solenoids_trim) == 1)/LabVIEWData.notes.analogSamplingRate_Hz;
    MergedData.data.solenoids.RPadSol = find(diff(LabVIEWData.data.solenoids_trim) == 2)/LabVIEWData.notes.analogSamplingRate_Hz;
    MergedData.data.solenoids.AudSol = find(diff(LabVIEWData.data.solenoids_trim) == 3)/LabVIEWData.notes.analogSamplingRate_Hz;  
    % Pull the notes and data from MScan
    MergedData.notes.MScan = MScanData.notes;
    MergedData.data.rawCorticalNeural = MScanData.data.corticalNeural_trim;
    MergedData.data.rawHippocampalNeural = MScanData.data.hippocampalNeural_trim;
    MergedData.data.corticalNeural.muaPower = MScanData.data.cortical.muaPower_trim;
    MergedData.data.corticalNeural.gammaBandPower = MScanData.data.cortical.gammaBandPower_trim;
    MergedData.data.corticalNeural.betaBandPower = MScanData.data.cortical.betaBandPower_trim;
    MergedData.data.corticalNeural.alphaBandPower = MScanData.data.cortical.alphaBandPower_trim;
    MergedData.data.corticalNeural.thetaBandPower = MScanData.data.cortical.thetaBandPower_trim;
    MergedData.data.corticalNeural.deltaBandPower = MScanData.data.cortical.deltaBandPower_trim;
    MergedData.data.hippocampalNeural.muaPower = MScanData.data.hippocampal.muaPower_trim;
    MergedData.data.hippocampalNeural.gammaBandPower = MScanData.data.hippocampal.gammaBandPower_trim;
    MergedData.data.hippocampalNeural.betaBandPower = MScanData.data.hippocampal.betaBandPower_trim;
    MergedData.data.hippocampalNeural.alphaBandPower = MScanData.data.hippocampal.alphaBandPower_trim;
    MergedData.data.hippocampalNeural.thetaBandPower = MScanData.data.hippocampal.thetaBandPower_trim;
    MergedData.data.hippocampalNeural.deltaBandPower = MScanData.data.hippocampal.deltaBandPower_trim;
    MergedData.data.forceSensorM = MScanData.data.dsForceSensorM_trim;
    MergedData.data.EMG.data = MScanData.data.filtEMG_trim;
    MergedData.data.binForceSensorM = MScanData.data.binForceSensorM_trim;
    MergedData.data.vesselDiameter.data = MScanData.data.vesselDiameter_trim; 
    % Most useful notes to be referenced in future analysis
    MergedData.notes.trialDuration_Sec = LabVIEWData.notes.trialDuration_Seconds_trim;
    if strcmp(MScanData.notes.movieType,'MC') == true
        MergedData.notes.p2Fs = MScanData.notes.p2Fs;
    else
        MergedData.notes.p2Fs = MScanData.notes.frameRate;
    end
    MergedData.notes.dsFs = MScanData.notes.dsFs;
    MergedData.notes.anFs = LabVIEWData.notes.analogSamplingRate_Hz;   
    save([animalID '_' hem '_' fileID '_' imageID '_' vesselID '_MergedData'],'MergedData')
end

end
