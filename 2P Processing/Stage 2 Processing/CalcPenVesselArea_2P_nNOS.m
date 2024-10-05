function [MScanData] = CalcPenVesselArea_2P(MScanData,fileID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Extract cross-sectional area from a penetrating vessel TIFF stack
%________________________________________________________________________________________________________________________

% takes a tiff file (movie) and an analog ascii file, and extracts the diameters
MScan_analogData = [fileID '.TXT'];
disp(['Loading MScan file: ' MScan_analogData '...']); disp(' ');
analogData = load(MScan_analogData, '-ascii');
MScanData.data.corticalNeural = analogData(:,2);
MScanData.data.EMG = analogData(:,3);
MScanData.data.hippocampalNeural = analogData(:,4);
MScanData.data.forceSensor = analogData(:,5);
MScanData.notes.analogSamplingRate = 20000;
MScanData.notes.firstFrame = double(imread(fileID,'TIFF','Index',1));
for n = 2:5
    MScanData.notes.firstFrame = MScanData.notes.firstFrame+double(imread(fileID,'tif','Index',n));
end
MScanData.notes.firstFrame = MScanData.notes.firstFrame/5;
fileInfo = imfinfo([fileID '.TIF']);
nframes = length(fileInfo);
angles = 1:1:180;
rtdThrehold = .5;   %0.5 'half max'
irtdThrehold = .2;   %0.2 gets rid of the noise from the inverse radon
area = -ones(min(nframes,MScanData.notes.numberOfFrames),1);
radoncontours = cell(min(nframes,MScanData.notes.numberOfFrames),1);
fftFirstFrame = fft2(MScanData.notes.firstFrame); % *
for f = 1:min(nframes,MScanData.notes.numberOfFrames)
    rawHoldImage = double(sum(imread(fileID,'tif','Index',f),3));
    fftRawHoldFrame = fft2(rawHoldImage); % *
    [MScanData.notes.pixelShift(:,f),~] = DftRegistration_2P(fftFirstFrame,fftRawHoldFrame,1);
    holdImage = rawHoldImage(round(MScanData.notes.vesselROI.boxPosition.xy(2):MScanData.notes.vesselROI.boxPosition.xy(2) + MScanData.notes.vesselROI.boxPosition.xy(4)),...
        round(MScanData.notes.vesselROI.boxPosition.xy(1):MScanData.notes.vesselROI.boxPosition.xy(1)+MScanData.notes.vesselROI.boxPosition.xy(3)));
    holdImage = holdImage-mean(holdImage(:));
    radonHoldImage = radon(holdImage,angles);
    for k = 1:length(angles)
        % normalize so that for each agle, the transfomed image is between
        % 0 and 1
        radonHoldImage(:,k) = radonHoldImage(:,k) - min(radonHoldImage(:,k));
        radonHoldImage(:,k) = radonHoldImage(:,k)/max(radonHoldImage(:,k));
        % find the peak of the projection at this angle
        [~,maxpointlocation(k)] = max(radonHoldImage(:,k));
        % threshold at half maximum
        try
            [~,minEdge(k)] = max(find(radonHoldImage(1:maxpointlocation(k),k) < rtdThrehold)); %#ok<*MXFND>
            [~,maxEdge(k)] = max(find(radonHoldImage(maxpointlocation(k)+1:end,k) > rtdThrehold));
        catch
            minEdge(k) = 0;
            maxEdge(k) = 0;
        end
        radonHoldImage(1:minEdge(k),k) = 0;
        radonHoldImage((maxpointlocation(k) + maxEdge(k)):end,k) = 0;
    end
    % transform the threhlded image back to an x,y image
    irtdNorm = iradon(double(radonHoldImage > rtdThrehold*max(radonHoldImage(:))),(angles),'linear','Hamming');
    %draw boundaries
    [cc,l] = bwboundaries(irtdNorm > irtdThrehold*max(irtdNorm(:)));
    %find the biggest area
    numPixels = cellfun(@length,cc);
    [~,idx] = max(numPixels);
    %find the nuber of pixes enclosed by the largest area
    areaFilled = regionprops(l,'FilledArea','Image','FilledImage');
    area(f) = length(find(areaFilled(idx).FilledImage));
    radoncontours{f} = contourc(double(irtdNorm(1:end,1:end)),double([irtdThrehold irtdThrehold]*max(irtdNorm(:))));
    % plot the raw image, the inverese of the radon transfoermed imaage, and the contours
    MScanData.notes.vesselROI.pixelArea = area;
    MScanData.notes.vesselROI.radoncontours = radoncontours;
    MScanData.notes.vesselROI.area = MScanData.notes.vesselROI.pixelArea*MScanData.notes.xFactor*MScanData.notes.xFactor;
    [holdHist,d] = hist(MScanData.notes.vesselROI.area,0:10:10000); %#ok<HIST>
    [~,maxD] = max(holdHist);
    MScanData.notes.vesselROI.modal_fixed_area = d(maxD);
    MScanData.data.vesselDiameter=2*sqrt(MScanData.notes.vesselROI.area/(pi))';
    MScanData.notes.vesselROI.modal_fixed_diameter=2*sqrt(MScanData.notes.vesselROI.modal_fixed_area/(pi));
end

end
