
% this function plots the body parts in heatmap, vector arrow, and
% population plot.
% go to folder named Population_compare/Session1 and select any of the
% animal named files i.e. SK01_Run1. you can do this for Session2 and 3
% also. once the data is loaded into the work directory run this code. this
% will generate the plot for that animal.
%         close all;
%% load a mouse and make the plots
%         clear all; close all; clc
function vectormap_plot(rootFolder,saveFigs,delim)
% %% first mouse
%         load('T261_Run1.mat');
%         animalID = OFTData.MouseID;
%         frameX = 1;
%         %%
%         % mouse body parts
%         Head = OFTData.TrackData(:,5:6);
%         mid = OFTData.TrackData(:,7:8);
%         Back = OFTData.TrackData(:,9:10);
%
%         tailTrunk = OFTData.TrackData(:,11:12);
%         tailMid = OFTData.TrackData(:,13:14);
%         tailEnd = OFTData.TrackData(:,15:16);
%
%         % outer
%         OTR = OFTData.TrackData(:,17:18);
%         OBR = OFTData.TrackData(:,19:20);
%         OTL = OFTData.TrackData(:,21:22);
%         OBL = OFTData.TrackData(:,23:24);
%         %%
%         Zodd = 1:2:9;
%         Zeven = 2:2:10;
%         outerpoints = [OTR,OTL,OBL,OBR,OTR];
%         %% Plot vector arrow
%         figure;
%         plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'*-r')
%         hold on;
%
%         A(1) =  mid(1,1);
%         A(2) =  mid(1,2);
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         p1 = plot(A(1),A(2),'ko','MarkerSize',10);
%         for i=1:1:1200
%         B(1) =  mid(i,1);
%         B(2) =  mid(i,2);
%
%         A(1) =  tailTrunk(i,1);
%         A(2) =  tailTrunk(i,2);
%         hold on;
%         vectarrow(A,B)
%         end
%         hold on
%         A(1) =  mid(1201,1);
%         A(2) =  mid(1201,2);
%         p2 = plot(A(1),A(2),'g*','MarkerSize',10);
%         axis equal
%         set(gca,'xticklabel',{[]},'yticklabel',{[]})
%         legend([p1,p2],'start','end')
%         ylabel('30 cm wide')
%         xlabel('60 cm long')
%         axis([-50 1000 -50 600])
%         xticks([])
%         yticks([])
%         title([animalID 'Vector Plot'])

%% plot another mouse from injected group
load('T267_Run1.mat');
animalID = OFTData.MouseID;
frameX = 1;
%%
% mouse body parts
Head = OFTData.TrackData(:,5:6);
mid = OFTData.TrackData(:,7:8);
Back = OFTData.TrackData(:,9:10);

tailTrunk = OFTData.TrackData(:,11:12);
tailMid = OFTData.TrackData(:,13:14);
tailEnd = OFTData.TrackData(:,15:16);

% outer
OTR = OFTData.TrackData(:,17:18);
OBR = OFTData.TrackData(:,19:20);
OTL = OFTData.TrackData(:,21:22);
OBL = OFTData.TrackData(:,23:24);
%%
Zodd = 1:2:9;
Zeven = 2:2:10;
outerpoints = [OTR,OTL,OBL,OBR,OTR];
%% Plot vector arrow
summaryFigure = figure;
plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'*-r')
hold on;

A(1) =  mid(1,1);
A(2) =  mid(1,2);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
p1 = plot(A(1),A(2),'ko','MarkerSize',10);
for i=1:1:1200
    B(1) =  mid(i,1);
    B(2) =  mid(i,2);

    A(1) =  tailTrunk(i,1);
    A(2) =  tailTrunk(i,2);
    hold on;
    vectarrow(A,B)
end
hold on
A(1) =  mid(1201,1);
A(2) =  mid(1201,2);
p2 = plot(A(1),A(2),'g*','MarkerSize',10);
axis equal
set(gca,'xticklabel',{[]},'yticklabel',{[]})
legend([p1,p2],'start','end')
ylabel('30 cm wide')
xlabel('60 cm long')
axis([-50 1000 -50 600])
xticks([])
yticks([])
title([animalID 'Vector Plot'])

if saveFigs == true
    % statistical diary
    dirpath = [rootFolder delim 'Summary Figures' delim 'Behavior' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'OpenField_VectorMapExample']);
end