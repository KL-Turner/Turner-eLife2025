function [ROIs] = CreateFreeHandROIs_IOS_eLife2025(img,ROIname,animalID,ROIs)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% create figure of the image frame
roiFig = figure;
imagesc(img)
colormap(gray)
axis image
xlabel('Caudal')
% draw ROI over the cement
disp(['Please select your region of interest for ' animalID ' ' ROIname '.']); disp(' ')
[~,xi,yi] = roipoly;
ROIs.(ROIname).xi = xi;
ROIs.(ROIname).yi = yi;
close(roiFig)