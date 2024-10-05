function [newVesselDiameter] = RemoveMotion_2P(vesselDiameter,baseline,diffThresh,rawThresh)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Motion correct vessel diameter. Performs best in x-y plane, not z.
%________________________________________________________________________________________________________________________

% index points for motion removal
indx1 = find(diff(vesselDiameter) > diffThresh);
indx2 = find(abs((vesselDiameter - baseline)/baseline) > (rawThresh));
indx = union(indx1 + 1,indx2);   % indx: points need to be interpolated
indx0 = 1:length(vesselDiameter);
indx0(indx) = [];   % indx0: good points
count = 1;
if isempty(indx) ~= 1
    if indx0(1) ~= 1
        indx0 = [1:indx0(1) - 1,indx0];
    end
    for a = 1:length(indx0) - 1
        step = indx0(a + 1) - indx0(a);
        if step == 1
            newVesselDiameter(count) = vesselDiameter(count); %#ok<*AGROW>
        end
        if (step ~= 1)
            newVesselDiameter(count) = vesselDiameter(count);
            newVesselDiameter(count + 1:count + step - 1) = (vesselDiameter(indx0(a + 1)) + vesselDiameter(indx0(a)))/2;
        end
        count = count + step;
    end
    newVesselDiameter(count) = vesselDiameter(indx0(end));
    if indx(end) == length(vesselDiameter)
        newVesselDiameter(indx0(end) + 1:length(vesselDiameter)) = vesselDiameter(indx0(end));
    end
else
    newVesselDiameter = vesselDiameter;
end

end
