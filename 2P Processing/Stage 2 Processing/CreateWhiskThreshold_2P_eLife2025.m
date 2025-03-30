function [thresh1,thresh2] = CreateWhiskThreshold_2P_nNOS(angle,fs)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Set the threshold for whisker acceleration to be used to binarize whisker movement.
%________________________________________________________________________________________________________________________

isok = 'n';
dd_wwf = abs((diff(angle, 2)))*fs^2;
whiskThresh = figure;
while strcmp(isok,'y') == 0
    plot(dd_wwf,'k');
    drawnow
    thresh2 = input('No Threshold for volitional whisks found. Please enter a threshold: '); disp(' ')
    thresh1 = input('No Threshold for resting behavior found. Please enter a threshold: '); disp(' ')
    bin_wwf = BinarizeWhiskers_2P_nNOS(angle,fs,thresh1,thresh2);
    ax1 = subplot(3,1,1);
    plot(angle,'k');
    axis tight
    set(gca,'box','off')
    ylabel('Angle')
    ax2 = subplot(3,1,2);
    plot(abs(diff(angle,2))*fs^2,'k');
    axis tight
    set(gca,'box','off')
    ylabel('Second Derivative')
    ax3 = subplot(3,1,3);
    plot(bin_wwf,'k');
    axis tight
    set(gca,'box','off')
    ylabel('Binarization')
    linkaxes([ax1,ax2,ax3],'x');
    drawnow
    isok = input('Is this threshold okay? (y/n) ','s'); disp(' ')
end
close(whiskThresh);

end

