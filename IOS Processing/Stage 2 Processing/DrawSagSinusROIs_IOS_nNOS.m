function [ROIs] = DrawSagSinusROIs_IOS_nNOS(animalID,fileID,ROIs,imagingWavelengths)
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
% draw ROI for the mask over the entire windows
isok = false;
while isok == false
    disp('Draw an ROI over the superior sagital sinus'); disp(' ')
    [~,rect] = imcrop;
    hold on;
    ROIoutline = rectangle('Position',rect,'EdgeColor','r');
    checkMask = input('Is the ROI okay? (y/n): ','s'); disp(' ')
    if strcmp(checkMask,'y') == true
        isok = true;
        ROIs.(['SSS_' strDay]).rect = rect;
    end
    delete(ROIoutline);
end

width = rect(3)/2;
lrect = rect;
lrect(1) = lrect(1) - width;
rrect = rect;
rrect(1) = rrect(1) + width;

rectangle('Position',lrect,'EdgeColor','r');
ROIs.(['lSSS_' strDay]).xi = [lrect(1),lrect(1) + lrect(3),lrect(1) + lrect(3),lrect(1)];
ROIs.(['lSSS_' strDay]).yi = [lrect(2),lrect(2),lrect(2) + lrect(4),lrect(2) + lrect(4)] ;

rectangle('Position',rect,'EdgeColor','b');
ROIs.(['SSS_' strDay]).xi = [rect(1),rect(1) + rect(3),rect(1) + rect(3),rect(1)];
ROIs.(['SSS_' strDay]).yi = [rect(2),rect(2),rect(2) + rect(4),rect(2) + rect(4)] ;

rectangle('Position',rrect,'EdgeColor','g');
ROIs.(['rSSS_' strDay]).xi = [rrect(1),rrect(1) + rrect(3),rrect(1) + rrect(3),rrect(1)];
ROIs.(['rSSS_' strDay]).yi = [rrect(2),rrect(2),rrect(2) + rrect(4),rrect(2) + rrect(4)] ;

% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% go through each file and check the color o
for qq = 1:size(procDataFileIDs,1)
    load(procDataFileIDs(qq,:));
    ProcData.notes.blueFrames = 1;
    ProcData.notes.redFrames = 2;
    ProcData.notes.greenFrames = 3;
    save(procDataFileIDs(qq,:),'ProcData')
end
close(roiFig)