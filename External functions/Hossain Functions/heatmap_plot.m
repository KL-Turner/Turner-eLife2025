        
% this function plots the body parts in heatmap. 

        %% load a mouse and make the plots
%         clear all; close all; clc
        function heatmap_plot()
        load('T261_Run1.mat');
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
        %% plot heatmap
        resolution = 1; % adjust to change n x n binning
        cameraX = 953/resolution;
        cameraY = 551/resolution;
        popMat = zeros(cameraY,cameraX);
        headX = round(mid(:,1)/resolution);
        headY = round(mid(:,2)/resolution);
        % populate popMat matrix
        for aa = 1:length(headX)
        headX_index = headX(aa,1);
        headY_index = headY(aa,1);
        popMat(headY_index,headX_index) = popMat(headY_index,headX_index) + 1; 
        end
        figure; 
        imagesc(popMat);
        hold on;
        % plot outline of maze
        plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'-r','LineWidth',2)
        
        axis image
        axis xy
        axis off
        caxis([0,1]) % adjust bin counts
        ylabel('30 cm wide')
        xlabel('60 cm long')
        axis([-50 1000 -50 600])
        xticks([])
        yticks([])
        title([animalID 'heatmap'])
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
        %% plot heatmap
       resolution = 1; % adjust to change n x n binning
        cameraX = 953/resolution;
        cameraY = 551/resolution;
        popMat = zeros(cameraY,cameraX);
        headX = round(mid(:,1)/resolution);
        headY = round(mid(:,2)/resolution);
        % populate popMat matrix
        for aa = 1:length(headX)
        headX_index = headX(aa,1);
        headY_index = headY(aa,1);
        popMat(headY_index,headX_index) = popMat(headY_index,headX_index) + 1; 
        end
        figure; 
        imagesc(popMat);
        hold on;
        % plot outline of maze
        plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'-r','LineWidth',2)
        
        axis image
        axis xy
        axis off
        caxis([0,1]) % adjust bin counts
        ylabel('30 cm wide')
        xlabel('60 cm long')
        axis([-50 1000 -50 600])
        xticks([])
        yticks([])
        title([animalID 'heatmap'])