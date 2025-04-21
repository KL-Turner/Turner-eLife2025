function [MScanData] = DiamCalcSurfaceVessel_2P_eLife2025(MScanData,imageID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Draw an ROI for the surface vessel and fill in notes information.
%________________________________________________________________________________________________________________________

% take file info and extract magnification and frame rate
movieInfo = imfinfo([imageID '.TIF']);
MScanData.notes.fileName = movieInfo(1).Filename;
MScanData.notes.frameWidth = num2str(movieInfo(1).Width);
MScanData.notes.frameHeight = num2str(movieInfo(1).Height);
MScanData.notes.numberOfFrames = length(movieInfo);
MScanData.notes.xSize = str2double(MScanData.notes.frameWidth);
MScanData.notes.ySize = str2double(MScanData.notes.frameHeight);
% Read header and take further action based on header information
textHold = strread(movieInfo(1).ImageDescription,'%s','delimiter','\n'); %#ok<*DSTRRD>
mscanStrings = textscan(movieInfo(1).ImageDescription,'%s','Delimiter', ':');
MScanData.notes.scanMode = mscanStrings{1}{17};%
MScanData.notes.magnification = textHold{20}(strfind(textHold{20},': ') + 2:end);
MScanData.notes.rotation = textHold{19}(strfind(textHold{19},': ') + 2:end);
MScanData.notes.frameTime = (textHold{24}(strfind(textHold{24},': ') + 2:end - 3));
MScanData.notes.frameRate = 1/str2double(MScanData.notes.frameTime);
MScanData.notes.startframe = 1;
MScanData.notes.endframe = MScanData.notes.numberOfFrames;
if MScanData.notes.objectiveID == 1       %10X
    MScanData.notes.micronsPerPixel = 1.2953;
elseif MScanData.notes.objectiveID == 2   % Small 20X
    MScanData.notes.micronsPerPixel = 0.5595;
elseif MScanData.notes.objectiveID == 3   % Big 20X
    MScanData.notes.micronsPerPixel = 0.64;
elseif MScanData.notes.objectiveID == 4   % 40X
    MScanData.notes.micronsPerPixel = 0.3619;
elseif MScanData.notes.objectiveID == 5   % 16X
    MScanData.notes.micronsPerPixel = 0.825;
end
MScanData.notes.timePerLine = 1/(MScanData.notes.frameRate*str2double(MScanData.notes.frameHeight));
MScanData.notes.xFactor = MScanData.notes.micronsPerPixel/(str2double(MScanData.notes.magnification(1:end - 1)));
% generate figure for drawing ROIs
image = imread(imageID,'TIFF','Index',1);
yString = 'y';
theInput = 'n';
while strcmp(yString,theInput) == false
    vesROI = figure;
    imagesc(double(image))
    title([MScanData.notes.animalID ' ' MScanData.notes.date ' ' MScanData.notes.imageID])
    colormap('gray');
    axis image
    xlabel('pixels')
    ylabel('pixels')
    xSize = size(image,2);
    ySize = size(image,1);
    area = impoly(gca,[1,1;1,20;20,20;20,1]); %#ok<*IMPOLY>
    drawnow
    theInput = input('Is the diameter of the box ok? (y/n): ', 's'); disp(' ')
    getAPI = iptgetapi(area);
    MScanData.notes.vesselROI.boxPosition.xy = getAPI.getPosition();
    MScanData.notes.vesselROI.xSize = xSize;
    MScanData.notes.vesselROI.ySize = ySize;
end
yString = 'y';
theInput = 'n';
diamAxis = imline(gca,round(xSize*[.25,.75]),round(ySize*[.25,.75])); %#ok<*IMLINE>
while strcmp(yString,theInput) == false
    theInput = input('Is the line along the diameter axis ok? (y/n): ','s'); disp(' ')
    getAPI = iptgetapi(diamAxis);
    MScanData.notes.vesselROI.vesselLine.position.xy = getAPI.getPosition();
    if strcmp(yString,theInput) == true
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Vessel ROIs/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(vesROI,[dirpath MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_ROIs']);
    end
end
close(vesROI);

end
