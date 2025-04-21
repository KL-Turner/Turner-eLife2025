function [binForceSensor] = BinarizeForceSensor_IOS_eLife2025(forceSensor,thresh)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
y = hilbert(diff(forceSensor));
env = abs(y);
binForceSensor = gt(env,thresh);