function [thresh1,thresh2] = CreateWhiskThreshold_IOS(angl,fs,strDay)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
isok = 'n';
dd_wwf = abs((diff(angl,2)))*fs^2;
whiskThresh = figure;
while strcmp(isok,'y') == 0
    plot((1:length(dd_wwf))/fs,dd_wwf,'k');
    xlabel('Time (s)')
    ylabel('|Acceleration*fs^2|')
    set(get(gca,'XAxis'),'TickLength',[0,0])
    set(gca,'box','off')
    drawnow
    thresh1 = input(['Please enter ' strDay ' resting threshold from whisker tracking: ']); disp(' ')
    thresh2 = input(['Please enter ' strDay ' whisking threshold from whisker tracking: ']); disp(' ')
    bin_wwf = BinarizeWhiskers_IOS(angl,fs,thresh1,thresh2);
    % angle
    ax1 = subplot(3,1,1);
    plot((1:length(angl))/fs,angl,'k');
    axis tight;
    xlabel('Time (s)')
    ylabel('Angle')
    set(get(gca,'XAxis'),'TickLength',[0,0])
    set(gca,'box','off')
    % abs velocity^2
    ax2 = subplot(3,1,2);
    plot(((1:length(abs(diff(angl,2))*fs^2)))/fs,abs(diff(angl,2))*fs^2,'k');
    hold on
    y1 = yline(thresh1,'r');
    y2 = yline(thresh2,'b');
    axis tight;
    xlabel('Time (s)')
    ylabel('|Acceleration*fs^2|')
    legend([y1,y2],'Resting threshold','Whisking threshold')
    set(get(gca,'XAxis'),'TickLength',[0,0])
    set(gca,'box','off')
    % binarization
    ax3 = subplot(3,1,3);
    plot((1:length(bin_wwf))/fs,bin_wwf,'k');
    axis tight;
    xlabel('Time (s)')
    ylabel('Binarization')
    set(get(gca,'XAxis'),'TickLength',[0,0])
    set(gca,'box','off')
    linkaxes([ax1,ax2,ax3],'x');
    isok = input('Is this threshold okay? (y/n): ','s'); disp(' ')
end
close(whiskThresh);