function [MScanData] = CapLinesScan_2P_eLife2025(MScanData,imageID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Draw edge points for linescan velocity calculation
%________________________________________________________________________________________________________________________

% take file info and extract magnification and frame rate
movieInfo = imfinfo([imageID '.TIF']);
MScanData.notes.filename = movieInfo(1).Filename;
MScanData.notes.frameWidth = num2str(movieInfo(1).Width);
MScanData.notes.frameHeight = num2str(movieInfo(1).Height);
MScanData.notes.numberOfFrames = length(movieInfo);
MScanData.notes.xSize = str2double(MScanData.notes.frameWidth);
MScanData.notes.ySize = str2double(MScanData.notes.frameHeight);
% read header and take further action based on header information
textHold = strread(movieInfo(1).ImageDescription,'%s','delimiter','\n'); %#ok<DSTRRD>
mscanStrings = textscan(movieInfo(1).ImageDescription,'%s','Delimiter', ':');
MScanData.notes.scanMode = mscanStrings{1}{17};%
MScanData.notes.magnification = str2double(textHold{20}(strfind(textHold{20},': ') + 2:end - 1));
MScanData.notes.rotation = textHold{19}(strfind(textHold{19},': ')+2:end);
MScanData.notes.frameTime = str2double((textHold{24}(strfind(textHold{24},': ') + 2:end - 3)));
MScanData.notes.frameRate = 1/MScanData.notes.frameTime;
MScanData.notes.startFrame = 1;
MScanData.notes.endFrame = MScanData.notes.numberOfFrames;
% magnification uM per pixel
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
MScanData.notes.pixelClock = 1/((5/4)*MScanData.notes.frameRate*MScanData.notes.xSize*MScanData.notes.ySize*.05e-6);
MScanData.notes.timePerLine = 1/(MScanData.notes.frameRate*str2double(MScanData.notes.frameHeight));
MScanData.notes.lineRate = 1/MScanData.notes.timePerLine;
MScanData.notes.tFactor = MScanData.notes.lineRate;
MScanData.notes.xFactor = MScanData.notes.micronsPerPixel/MScanData.notes.magnification;
MScanData.notes.decimate = 1;
% generate figure for drawing ROIs
framesHold = LoadTiffConcatenate_2P_eLife2025([imageID '.TIF'],[1,5]);
yString = 'y';
theInput = 'n';
while strcmp(yString,theInput) == false
    lineScanBoundaries = figure;
    imagesc(double(framesHold))
    title([MScanData.notes.animalID ' ' MScanData.notes.date ' ' MScanData.notes.imageID])
    colormap('gray');
    axis image
    axis off
    % draw boundary lines for linescan velocity calculation
    disp('Select left and right boundary for velocity calculation range'); disp(' ')
    [x,~] = ginput(2);
    MScanData.notes.xStart = round(min(x));
    MScanData.notes.xStop = round(max(x));
    hold on;
    xline(MScanData.notes.xStart,'r')
    xline(MScanData.notes.xStop,'r')
    theInput = input('Are the boundaries okay? (y/n): ','s'); disp(' ')
    if strcmp(yString,theInput) == true
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Vessel ROIs/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(lineScanBoundaries,[dirpath MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_ROIs']);
    end
    close(lineScanBoundaries);
end

end
