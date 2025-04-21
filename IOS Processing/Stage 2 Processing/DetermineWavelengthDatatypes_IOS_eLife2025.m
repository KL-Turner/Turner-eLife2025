function [dataTypes] = DetermineWavelengthDatatypes_IOS_eLife2025(imagingWavelengths,iteration)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
if strcmp(imagingWavelengths,'Red, Green, & Blue') == true
    if iteration == 1
        dataTypes = {'green','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 2
        dataTypes = {'green','blue','red','HbT','HbO','HbR','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 3
        dataTypes = {'HbT','HbO','HbR','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil','heartRate','whiskerAngle','forceSensor'};
    end
elseif strcmp(imagingWavelengths,'Red, Lime, & Blue') == true
    if iteration == 1
        dataTypes = {'lime','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 2
        dataTypes = {'lime','blue','red','HbT','HbO','HbR','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 3
        dataTypes = {'HbT','HbO','HbR','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil','heartRate','whiskerAngle','forceSensor'};
    end
elseif strcmp(imagingWavelengths,'Green & Blue') == true
    if iteration == 1
        dataTypes = {'green','blue','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 2
        dataTypes = {'green','blue','HbT','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 3
        dataTypes = {'HbT','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil','heartRate','whiskerAngle','forceSensor'};
    end
elseif strcmp(imagingWavelengths,'Lime & Blue') == true
    if iteration == 1
        dataTypes = {'lime','blue','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 2
        dataTypes = {'lime','blue','HbT','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 3
        dataTypes = {'HbT','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil','heartRate','whiskerAngle','forceSensor'};
    end
elseif strcmp(imagingWavelengths,'Green') == true
    if iteration == 1
        dataTypes = {'green','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 2
        dataTypes = {'green','HbT','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 3
        dataTypes = {'HbT','cortical_LH','cortical_RH','hippocampus','EMG','pupil','heartRate','whiskerAngle','forceSensor'};
    end
elseif strcmp(imagingWavelengths,'Lime') == true
    if iteration == 1
        dataTypes = {'lime','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 2
        dataTypes = {'lime','HbT','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 3
        dataTypes = {'HbT','cortical_LH','cortical_RH','hippocampus','EMG','pupil','heartRate','whiskerAngle','forceSensor'};
    end
elseif strcmp(imagingWavelengths,'Blue') == true
    if iteration == 1
        dataTypes = {'blue','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 2
        dataTypes = {'blue','HbT','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    elseif iteration == 3
        dataTypes = {'HbT','cortical_LH','cortical_RH','hippocampus','EMG','pupil','heartRate','whiskerAngle','forceSensor'};
    end
end