function [thresh] = CreateForceSensorThreshold_2P_eLife2025(forceSensor)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: View the force sensor data and determine what a good value is to binarize movement.
%________________________________________________________________________________________________________________________

y = hilbert(diff(forceSensor));
force = abs(y);
forceThresh = figure;
isok = 'n';
while strcmp(isok,'y') == 0
    plot(force,'k');
    drawnow
    thresh = input('No Threshold to binarize pressure sensor found. Please enter a threshold: '); disp(' ')
    binForceSensor = BinarizeForceSensor_2P_eLife2025(forceSensor,thresh);
    subplot(3,1,1)
    plot(forceSensor,'k')
    axis tight
    set(gca,'box','off')
    subplot(3,1,2)
    plot(force,'k')
    axis tight
    set(gca,'box','off')
    subplot(3,1,3)
    plot(binForceSensor,'k')
    axis tight
    set(gca,'box','off')
    drawnow
    isok = input('Is this threshold okay? (y/n) ','s'); disp(' ')
end
close(forceThresh);

end
