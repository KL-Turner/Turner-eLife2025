clc; clear all; close all
delete(('OFTTable.mat'))
% function CorrectedOFTAnalysis(subjectPath)
% cd(subjectPath)
%% load the data
MatFiles = ls('*.mat'); % list the files in the folder
for csvList = 1:1:size(MatFiles,1)
        fileName = MatFiles(csvList,:);
%         extInd = strfind(fileName(csvList,:),'.'); 
%         matfileName = [fileName(1:extInd) 'mat']; 
        load(fileName);

        OFTData = rmfield(OFTData,'TrackAnalysis');
        %% Assign values
        
        % Short name alphabets
        % O --> Outer
        % C --> Center
        % T --> Top
        % B --> Bottom
        % L --> Left
        % R --> Right
        
        % mouse body parts
        Head = OFTData.TrackData(:,5:6);
        Mid = OFTData.TrackData(:,7:8);
        Back = OFTData.TrackData(:,9:10);

        tailTrunk = OFTData.TrackData(:,11:12);
        tailMid = OFTData.TrackData(:,13:14);
        tailEnd = OFTData.TrackData(:,15:16);
        
        % non mouse bodyparts
        % outer 
        OTR = OFTData.TrackData(:,17:18);
        OBR = OFTData.TrackData(:,19:20);
        OTL = OFTData.TrackData(:,21:22);
        OBL = OFTData.TrackData(:,23:24);
         % ears
        LEar = OFTData.TrackData(:,1:2);
        REar = OFTData.TrackData(:,3:4);
        %% calculate frame rate
        OFTData.TrackAnalysis.totalframes = length( LEar);
        totalTime = 10*60; % 10 minutes
        OFTData.TrackAnalysis.FrameRate = OFTData.TrackAnalysis.totalframes/totalTime;
        %% measure the pixel size
        OFTData.TrackAnalysis.PixelsSize = sqrt((OTL(1,1)-OTR(1,1)).^2 + (OTL(1,2)-OTR(1,2)).^2); % distance in pixels
        % the actual distance is 60cm.
        OFTData.TrackAnalysis.PixelstoCm = OFTData.TrackAnalysis.PixelsSize/60; 
        %% calculate the inner points
        CBR(:,1) = OBR(:,1) - (5*OFTData.TrackAnalysis.PixelstoCm);
        CBR(:,2) = OBR(:,2) - (5*OFTData.TrackAnalysis.PixelstoCm);
        CBL(:,1) = OBL(:,1) + (5*OFTData.TrackAnalysis.PixelstoCm);
        CBL(:,2) = OBL(:,2) - (5*OFTData.TrackAnalysis.PixelstoCm);
        CTL(:,1) = OTL(:,1) + (5*OFTData.TrackAnalysis.PixelstoCm);
        CTL(:,2) = OTL(:,2) + (5*OFTData.TrackAnalysis.PixelstoCm);
        CTR(:,1) = OTR(:,1) - (5*OFTData.TrackAnalysis.PixelstoCm);
        CTR(:,2) = OTR(:,2) + (5*OFTData.TrackAnalysis.PixelstoCm);
        %% plot the data
        if sum(strcmp(fileName, ["T208_Run1.mat" , "T217_Run1.mat" , "T269_Run1.mat" ,"T273_Run1.mat"]))
            frameX = 137;
        else
            frameX = 59;
        end
        %{
        Zodd = 1:2:9;
        Zeven = 2:2:10;
        outerpoints = [OTR,OTL,OBL,OBR,OTR];
        innerpoints = [CTR,CTL,CBL,CBR,CTR];
        figure;       
        plot(outerpoints(frameX,Zodd),outerpoints(frameX,Zeven),'*-r')
        hold on
        plot(innerpoints(frameX,Zodd),innerpoints(frameX,Zeven),'o-k')
        plot(Mid(:,1),Mid(:,2),'.r')

        axis([-50 1000 -50 600])
        set(gca,'xticklabel',{[]},'yticklabel',{[]})
        xlabel({[]})
        ylabel({[]})
        str1 = strfind(fileName,'.');
        FigName = [fileName(1:str1(1)) 'png'];
        saveas(gcf,FigName)
        close all
        %}
        %% Calculate Displacement
        for mm = 2:1:size(Head,1)
            OFTData.TrackAnalysis.displ(mm-1,1) = sqrt((Head(mm,1)-Head(mm-1,1)).^2 + (Head(mm,2)-Head(mm-1,2)).^2);
            OFTData.TrackAnalysis.displ(mm-1,2) = sqrt((Mid(mm,1)-Mid(mm-1,1)).^2 + (Mid(mm,2)-Mid(mm-1,2)).^2);
            OFTData.TrackAnalysis.displ(mm-1,3) = sqrt((tailTrunk(mm,1)-tailTrunk(mm-1,1)).^2 + (tailTrunk(mm,2)-tailTrunk(mm-1,2)).^2);
        end
        OFTData.TrackAnalysis.disCm = OFTData.TrackAnalysis.displ./OFTData.TrackAnalysis.PixelstoCm; % convert to cm
        OFTData.TrackAnalysis.avg_dis_cms =  OFTData.TrackAnalysis.FrameRate*mean(OFTData.TrackAnalysis.disCm); % average displacement in cm per seconds      
        %% dintance binned in minutes
        DataSplit = floor(length(OFTData.TrackAnalysis.disCm)/(30));
        
            for isplit = 1:1:10
                startInd = DataSplit*(isplit-1)+1;
                endInd = DataSplit*isplit;
                OFTData.TrackAnalysis.tot_dis_binned(isplit,:) =  sum(OFTData.TrackAnalysis.disCm(startInd:endInd,:)); % distance per minute in cm
            end
        OFTData.TrackAnalysis.tot_dis =  sum(OFTData.TrackAnalysis.tot_dis_binned,1); % total displacement in cm
        OFTData.TrackAnalysis.tot_dis_5min =  sum(OFTData.TrackAnalysis.tot_dis_binned((1:end/2),2)); % total displacement in cm
        %% Calcule center time     
        center_poly = [CTL,CTR,CBR,CBL];
        [OFTData.TrackAnalysis.center_in,center_on] = inpolygon(Mid(:,1),Mid(:,2),center_poly(frameX,[1,3,5,7]),center_poly(frameX,[2,4,6,8]));
        OFTData.TrackAnalysis.center_time = sum(OFTData.TrackAnalysis.center_in)/OFTData.TrackAnalysis.FrameRate;
        OFTData.TrackAnalysis.center_time_percentage = sum(OFTData.TrackAnalysis.center_in)/length(OFTData.TrackAnalysis.center_in);

        OFTData.TrackAnalysis.center_time_5min = sum(OFTData.TrackAnalysis.center_in(1:end/2))/OFTData.TrackAnalysis.FrameRate;
        OFTData.TrackAnalysis.center_time_percentage_5min = sum(OFTData.TrackAnalysis.center_in(1:end/2))/length(OFTData.TrackAnalysis.center_in(1:end/2));
        %% calculate time spend per minutes and total
%         DataSplit = floor(length(OFTData.TrackAnalysis.center_in)/(10));
%         
%             for isplit = 1:1:10
%                 startInd = DataSplit*(isplit-1)+1;
%                 endInd = DataSplit*isplit;
%                 OFTData.TrackAnalysis.timespend_center_binned(isplit) =  sum(double(OFTData.TrackAnalysis.center_in(startInd:endInd)))/OFTData.TrackAnalysis.FrameRate; % distance per minute in cm
%                 OFTData.TrackAnalysis.timespend_outside_binned(isplit) =  sum(double(OFTData.TrackAnalysis.outside_in(startInd:endInd)))/OFTData.TrackAnalysis.FrameRate; % distance per minute in cm
%             end
%         OFTData.TrackAnalysis.timespend_center_total =  sum(OFTData.TrackAnalysis.timespend_center_binned); % total displacement in cm
%         OFTData.TrackAnalysis.timespend_outside_total =  sum(OFTData.TrackAnalysis.timespend_outside_binned); % total displacement in cm  
        %% Calculate latency to first center entries
% 
%         [OFTData.TrackAnalysis.CenterEntires, OFTData.TrackAnalysis.FramesCenterEntries] = findpeaks(double(OFTData.TrackAnalysis.center_in));
%         if size(OFTData.TrackAnalysis.FramesCenterEntries,2)>0
%             if OFTData.TrackAnalysis.FramesCenterEntries(1) == 1 % exclude if the mouse start from center
%                 OFTData.TrackAnalysis.Latencyto_1stCenterEntry = OFTData.TrackAnalysis.FramesCenterEntries(2) / OFTData.TrackAnalysis.FrameRate;
%             else
%                 OFTData.TrackAnalysis.Latencyto_1stCenterEntry = OFTData.TrackAnalysis.FramesCenterEntries(1) / OFTData.TrackAnalysis.FrameRate;
%             end
%         elseif size(OFTData.TrackAnalysis.FramesCenterEntries,2)==0
%                 OFTData.TrackAnalysis.Latencyto_1stCenterEntry = length(OFTData.TrackAnalysis.center_in)/OFTData.TrackAnalysis.FrameRate;
%         end
        %%
        clearvars -EXCEPT OFTData.TrackAnalysis OFTData fileName csvList MatFiles 
        
        save(fileName,'OFTData')        
end