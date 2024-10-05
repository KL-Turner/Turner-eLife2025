function [fixedVelocity] = VelocityCleanUpSTDOutliers_2P_nNOS(velocity,threshold)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: finds points outside a threshold velocity and interpolated between last
%            and next good points, unless at the begining or end, then puts in mean velocities
%________________________________________________________________________________________________________________________

upperBound = threshold*std(velocity) + mean(velocity);
lowerBound = -threshold*std(velocity) + mean(velocity);
badPoints = union(find(velocity > upperBound),find(velocity < lowerBound));
nOutliers = length(badPoints);
fixedVelocity = velocity;
while nOutliers > 0
    goodPoints = union(setdiff(badPoints + 1,badPoints),setdiff(badPoints - 1,badPoints));
    for pp = 1:length(badPoints)
        lastGood = find(goodPoints < badPoints(pp),1,'last');
        nextGood = find(goodPoints > badPoints(pp),1);
        if ((goodPoints(1) > 0) && (goodPoints(end) <= length(fixedVelocity)))
            if ((numel(lastGood) > 0) && (isempty(nextGood) == 0))
                fixedVelocity(badPoints(pp)) = .5*fixedVelocity(goodPoints(lastGood)) + .5*fixedVelocity(goodPoints(nextGood));
            elseif ((numel(lastGood) > 0))
                fixedVelocity(badPoints(pp)) = fixedVelocity(goodPoints(nextGood));
            else
                fixedVelocity(badPoints(pp))=fixedVelocity(goodPoints(lastGood));
            end
        elseif (goodPoints(1) == 0)
            fixedVelocity(badPoints(pp)) = mean(fixedVelocity);
        elseif (goodPoints(end) == length(fixedVelocity) + 1)
            fixedVelocity(badPoints(pp)) = mean(fixedVelocity);
        elseif numel(goodPoints) == 0
            fixedVelocity(badPoints(pp)) = -fixedVelocity;
        end      
    end
    upperBound = threshold*std(fixedVelocity) + mean(fixedVelocity);
    lowerBound = -threshold*std(fixedVelocity) + mean(fixedVelocity);
    badPoints = union(find(fixedVelocity > upperBound),find(fixedVelocity < lowerBound));
    nOutliers = length(badPoints);   
end

end
