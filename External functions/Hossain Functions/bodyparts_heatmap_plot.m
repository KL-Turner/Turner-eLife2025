        
% this function plots the body parts in heatmap, vector arrow, and
% population plot. 
% go to folder named Population_compare/Session1 and select any of the
% animal named files i.e. SK01_Run1. you can do this for Session2 and 3
% also. once the data is loaded into the work directory run this code. this
% will generate the plot for that animal. 
%         close all;
        %% load a mouse and make the plots
        clear all; close all; clc
% function bodyparts_heatmap_plot()
%% first mouse        
        load('T247_Run1.mat');
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
        figure;       
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
        
        %% plot dotmap mid
%         figure;       
% 
%         plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'*-k')
%         hold on;
%         axis equal
%         set(gca,'xticklabel',{[]},'yticklabel',{[]})
%         xlabel({[]})
%         ylabel({[]})
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         plot(mid(:,1),mid(:,2),'.r')
%         axis([-50 1000 -50 600])
%         title([animalID 'dot map'])
        %% plot heatmap
       
%         resolution = 1; % adjust to change n x n binning
%         cameraX = round(953/resolution);
%         cameraY = round(551/resolution);
%         popMat = zeros(cameraY,cameraX);
%         headX = round(mid(:,1)/resolution);
%         headY = round(mid(:,2)/resolution);
%         % populate popMat matrix
%         for aa = 1:1:length(headX)
%         headX_index = headX(aa,1);
%         headY_index = headY(aa,1);
%         popMat(headY_index,headX_index) = popMat(headY_index,headX_index) + 1; 
%         end
%         figure; 
%         imagesc(popMat);
%         hold on;
%         % plot outline of maze
%         plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'-r','LineWidth',2)
%         
%         axis image
% 
%         ylabel('30 cm wide')
%         xlabel('60 cm long')
%         axis([-50 1000 -50 600])
%         xticks([])
%         yticks([])
%         axis xy
%         axis off
%         caxis([0,2]) % adjust bin counts
% %         axis([-50 1000 -50 600])
%         title([animalID 'heatmap'])
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
        figure;       
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
        %% plot dotmap mid
%         figure;       
% 
%         plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'*-k')
%         hold on;
%         axis equal
%         set(gca,'xticklabel',{[]},'yticklabel',{[]})
%         xlabel({[]})
%         ylabel({[]})
%         set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%         plot(mid(:,1),mid(:,2),'.r')
%         axis([-50 1000 -50 600])
%         title([animalID 'dot map'])
        %% plot heatmap
%         resolution = 1; % adjust to change n x n binning
%         cameraX = (953/resolution);
%         cameraY = (551/resolution);
%         popMat = zeros(cameraY,cameraX);
%         headX = round(mid(:,1)/resolution);
%         headY = round(mid(:,2)/resolution);
%         % populate popMat matrix
%         for aa = 1:length(headX)
%         headX_index = headX(aa,1);
%         headY_index = headY(aa,1);
%         popMat(headY_index,headX_index) = popMat(headY_index,headX_index) + 1; 
%         end
%         figure; 
%         imagesc(popMat);
%         hold on;
%         % plot outline of maze
%         plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'-r','LineWidth',2)
%         
%         axis image
%         % title('Heatmap')
%         ylabel('pixels')
%         xlabel('pixels')
%         axis xy
%         axis off
%         caxis([0,2]) % adjust bin counts
%         axis([-50 1000 -50 600])
%         title([animalID 'heatmap'])