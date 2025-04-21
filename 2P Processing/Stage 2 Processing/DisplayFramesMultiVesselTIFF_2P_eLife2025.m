function [imagefig,xstart,xstop,frames_hold,nframes,scan_type] = DisplayFramesMultiVesselTIFF_2P_eLife2025(fname,MScanData)

    nframes = MScanData.notes.frameCount;
    MScanData.xsize = str2double(MScanData.notes.frameWidth);
    MScanData.ysize = str2double(MScanData.notes.frameHeight);
    frame_height = MScanData.ysize;
 
    frames_hold = LoadTiffConcatenate_2P_eLife2025(fname,[1,5]);%
    imagefig = figure(2);
    colormap gray
    %plot the data
 
    hold off
    imagesc(double(frames_hold))
    hold on
    %axis image
    axis off
    title(fname)
    %plot the linescan trajectory
    if length(mfname)>0 %plot the
        try
            load(mfname)
            line_boundaries=find(abs(diff(scanData.pathObjSubNum)));
            scan_velocity=scanData.scanVelocity;
            for ss=1:length(line_boundaries)
                plot(line_boundaries(ss)*[1 1],[0 1000],'r')
            end
            subplot(212)
            hold off
            imagesc(scanData.im)
            axis image
            hold on
            pathImCoords(:,1) = scanData.path(:,1) * (size(scanData.im,2)-1)/(abs(diff(scanData.axisLimCol))) + 1 - (size(scanData.im,2)-1)/(abs(diff(scanData.axisLimCol)))*min(scanData.axisLimCol);
            pathImCoords(:,2) = scanData.path(:,2) * (size(scanData.im,1)-1)/(abs(diff(scanData.axisLimRow))) + 1 - (size(scanData.im,1)-1)/(abs(diff(scanData.axisLimRow)))*min(scanData.axisLimRow);
            plot(pathImCoords(:,1),pathImCoords(:,2),'r')

        end
    end
    for vessel=1:nvessels
        disp('Select velocity calculation range.')
        [x,y]=ginput(2);
        xstart(vessel)=round(min(x));
        xstop(vessel)=round(max(x));
        scan_type(vessel)=1;
    end
end
