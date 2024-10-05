function [images] = LoadTiffConcatenate_2P(filename,frameIndeces)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Pulls out frame number from tiff stack
%________________________________________________________________________________________________________________________

if isempty(frameIndeces) == true
    movieInfo = imfinfo(filename);
    startFrame = 1;
    endFrame = length(movieInfo);
    numnberOfFrames = endFrame - startFrame + 1;
else
    startFrame = min(frameIndeces);
    endFrame = max(frameIndeces);
    numnberOfFrames = endFrame - startFrame + 1;
end
imageFile = Tiff(filename);
firstFrame = imageFile.read();
imageHeight = size(firstFrame,1);
imageWidth = size(firstFrame,2);
images = zeros(imageHeight*numnberOfFrames,imageWidth);
for aa = 1:numnberOfFrames
    imageFile.setDirectory(aa);
    images((1 + (aa - 1)*imageHeight):(aa*imageHeight),:) = double(imageFile.read());
end
imageFile.close();

end
