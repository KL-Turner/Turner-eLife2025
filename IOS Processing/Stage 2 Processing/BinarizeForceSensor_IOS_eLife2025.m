function [binForceSensor] = BinarizeForceSensor_IOS_nNOS(forceSensor,thresh)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
y = hilbert(diff(forceSensor));
env = abs(y);
binForceSensor = gt(env,thresh);