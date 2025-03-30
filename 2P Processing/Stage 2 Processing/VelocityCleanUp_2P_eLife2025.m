function [fixedVelocity]=VelocityCleanUp_2P_nNOS(velocity,threshold)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: finds points faster than a threshold velocity and interpolated between last good points
%________________________________________________________________________________________________________________________

badpoints = find(abs(velocity) > threshold);
goodpoints = union(setdiff(badpoints + 1,badpoints),setdiff(badpoints - 1,badpoints));
fixedVelocity = velocity;
for pp = 1:length(badpoints)
    lastGood = find(goodpoints < badpoints(pp),1,'last');
    nextGood = find(goodpoints > badpoints(pp),1);
    if ((goodpoints(1) > 0) && (goodpoints(end) <= length(velocity)))
        if((numel(lastGood) > 0) && (isempty(nextGood) == 0))
            fixedVelocity(badpoints(pp)) = .5*velocity(goodpoints(lastGood)) + .5*velocity(goodpoints(nextGood));
        elseif((numel(lastGood) > 0))
            fixedVelocity(badpoints(pp)) = velocity(goodpoints(nextGood));
        else
            fixedVelocity(badpoints(pp)) = velocity(goodpoints(lastGood));
        end
    elseif (goodpoints(1) == 0)
        fixedVelocity(badpoints(pp)) = velocity(goodpoints(nextGood));
    elseif (goodpoints(end) == length(velocity) + 1)
        fixedVelocity(badpoints(pp)) = velocity(goodpoints(lastGood));
    elseif numel(goodpoints) == 0
        fixedVelocity(badpoints(pp)) = -velocity;
    end
end

end
