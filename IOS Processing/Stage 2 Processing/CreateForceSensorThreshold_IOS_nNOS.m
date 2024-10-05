function [thresh] = CreateForceSensorThreshold_IOS_nNOS(forceSensor,fs,strDay)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
isok = 'n';
force = abs(hilbert(diff(forceSensor)));
forceThresh = figure;
while strcmp(isok,'y') == 0
    plot((1:length(force))/fs,force,'k');
    xlabel('Time (s)')
    ylabel('|diff(Force sensor)|')
    set(get(gca,'XAxis'),'TickLength',[0,0])
    set(gca,'box','off')
    drawnow
    thresh = input(['Please enter ' strDay ' resting threshold from body motion: ']); disp(' ')
    binForceSensor = BinarizeForceSensor_IOS_nNOS(forceSensor,thresh);
    % force sensor V
    ax1 = subplot(3,1,1);
    plot((1:length(forceSensor))/fs,forceSensor,'k')
    axis tight;
    xlabel('Time (s)')
    ylabel('Force sensor (V)')
    set(get(gca,'XAxis'),'TickLength',[0,0])
    set(gca,'box','off')
    % diff (force)
    ax2 = subplot(3,1,2);
    plot((1:length(force))/fs,force,'k')
    hold on
    y1 = yline(thresh,'r');
    axis tight;
    xlabel('Time (s)')
    ylabel('|diff(Force sensor)|')
    legend(y1,'Resting threshold')
    set(get(gca,'XAxis'),'TickLength',[0,0])
    set(gca,'box','off')
    % binarization
    ax3 = subplot(3,1,3);
    plot((1:length(binForceSensor))/fs,binForceSensor,'k')
    axis tight;
    xlabel('Time (s)')
    ylabel('Binarization')
    set(get(gca,'XAxis'),'TickLength',[0,0])
    set(gca,'box','off')
    linkaxes([ax1,ax2,ax3],'x');
    isok = input('Is this threshold okay? (y/n): ','s'); disp(' ')
end
close(forceThresh);