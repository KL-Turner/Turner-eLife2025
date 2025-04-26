% this is the final function of the data analysis for the behavioral
% experiment. 
% function FinalOFTAnalysis()
delete(('OFTTable.mat'))
clc; clear all; close all
AnimalTable_OFT = table(); 

%% load the data
MatFiles = ls('*.mat'); % list the files in the folder
for csvList = 1:1:size(MatFiles,1)
        fileName = MatFiles(csvList,:);
        extInd = strfind(fileName,'_'); 
        load(fileName);
        
        tempTable = table(); 
        tempTable.AnimalID = fileName(1:extInd(1)-1);
        %% Calculate population average 
        % I calculated a variety of measurements for this tracking;
        % however, we decided to drop all of them after last presentation
        % and discussion. So we only have the total distance travelled
        % calculated for the final version.
        % distance traveled over 10 min (total and binned by minute)
        tempTable.total_distance_binned = OFTData.TrackAnalysis.tot_dis_binned(:,2)';
        

        % distance
        tempTable.avg_distance =  OFTData.TrackAnalysis.avg_dis_cms(2); % 
        tempTable.total_distance = OFTData.TrackAnalysis.tot_dis(:,2);
        % centerTime 
        tempTable.center_time =  OFTData.TrackAnalysis.center_time; % 
        % center percentage
        tempTable.center_percentage =  OFTData.TrackAnalysis.center_time_percentage*100; % 
% 5 min
        tempTable.total_distance_5min = OFTData.TrackAnalysis.tot_dis_5min;
        % centerTime 
        tempTable.center_time_5min =  OFTData.TrackAnalysis.center_time_5min; % 
        % center percentage
        tempTable.center_percentage_5min =  OFTData.TrackAnalysis.center_time_percentage_5min*100; % 

        % add data to main table
        AnimalTable_OFT = [AnimalTable_OFT;tempTable];

        clearvars -EXCEPT AnimalTable_OFT fileName csvList MatFiles         
end

    BlindKey = 'ExperimentRecords_Keys.xlsx'; % this is the blind key 

    T = readtable(BlindKey);

    BlindTable = table();

    BlindTable.Sex = categorical(T.Sex);
    BlindTable.DrugGroup = categorical(T.SxDescription);

    AnimalTable_OFT = [AnimalTable_OFT,BlindTable]; % now the table has information about the drug and sec group
    clearvars -EXCEPT AnimalTable_OFT
    save('OFTTable.mat','AnimalTable_OFT') 

    GLME_STATS = fitglme(AnimalTable_OFT,"center_percentage~1+Sex+DrugGroup") % we are running a GLME stats
    GLME_STATS = fitglme(AnimalTable_OFT,"total_distance_5min~1+Sex+DrugGroup") % we are running a GLME stats

    
figure;
    boxplot(AnimalTable_OFT.total_distance,AnimalTable_OFT.DrugGroup);
    xlabel('Experimental Groups');
    ylabel('Distance Travelled (cm)')
    title('Population comparison of total distance travelled')

figure;
    boxplot(AnimalTable_OFT.center_percentage,AnimalTable_OFT.DrugGroup);
    xlabel('Experimental Groups');
    ylabel('Pertentage time spent')
    title('Pertentage of time spent in center')

figure; % plot mean distance
    [G,GN] = grp2idx(AnimalTable_OFT.DrugGroup);
    Cidx = find(G==2);
    Bidx = find(G==1);
    Sidx = find(G==3);
    Distance_Control = AnimalTable_OFT.total_distance(Cidx,:);
    Distance_Blank = AnimalTable_OFT.total_distance(Bidx,:);
    Distance_SSP = AnimalTable_OFT.total_distance(Sidx,:);

    mean_distance_Control = mean(Distance_Control);
    mean_distance_Blank = mean(Distance_Blank);
    mean_distance_SSP = mean(Distance_SSP);
    serr_distance_control = std(Distance_Control)/sqrt(length(Distance_Control));
    serr_distance_blank = std(Distance_Blank)/sqrt(length(Distance_Blank));
    serr_distance_ssp = std(Distance_SSP)/sqrt(length(Distance_SSP));

    cplot=bar(1,mean_distance_Control,0.1); 
    cplot.FaceColor = [0.9290 0.6940 0.1250];
    hold on
    plot(1,Distance_Control,'ob','MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','k');
    errorbar(1,mean_distance_Control,serr_distance_control,'-k','CapSize',18,...
    'LineWidth',3);
   
    kplot= bar(1.25,mean_distance_Blank,0.1);
    kplot.FaceColor = [0.6350 0.0780 0.1840];
    plot(1.25,Distance_Blank,'go','MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','k');
    errorbar(1.25,mean_distance_Blank,serr_distance_blank,'-k','CapSize',18,...
    'LineWidth',3);
    
    nplot = bar(1.5,mean_distance_SSP,0.1);
    nplot.FaceColor = [0.4660 0.6740 0.1880];
    plot(1.5,Distance_SSP,'or','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k');
    errorbar(1.5,mean_distance_SSP,serr_distance_ssp,'-k','CapSize',18,...
    'LineWidth',3);
    xlabel('Experimental Groups');
    ylabel('Distance Travelled (cm)')
    title('Population comparison of total distance travelled')
    xlim([0.75 1.75])
    xticks([1 1.25 1.5]);
    xticklabels({'WT Control','Blank-SAP','SSP-SAP'})

figure; % plot center percentage
    [G,GN] = grp2idx(AnimalTable_OFT.DrugGroup);
    Cidx = find(G==2);
    Bidx = find(G==1);
    Sidx = find(G==3);
    Center_Control = AnimalTable_OFT.center_percentage(Cidx,:);
    Center_Blank = AnimalTable_OFT.center_percentage(Bidx,:);
    Center_SSP = AnimalTable_OFT.center_percentage(Sidx,:);

    mean_center_Control = mean(Center_Control);
    mean_center_Blank = mean(Center_Blank);
    mean_center_SSP = mean(Center_SSP);
    serr_center_control = std(Center_Control)/sqrt(length(Center_Control));
    serr_center_blank = std(Center_Blank)/sqrt(length(Center_Blank));
    serr_center_ssp = std(Center_SSP)/sqrt(length(Center_SSP));

    cplot=bar(1,mean_center_Control,0.1); 
    cplot.FaceColor = [0.9290 0.6940 0.1250];
    hold on
    plot(1,Center_Control,'ob','MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','k');
    errorbar(1,mean_center_Control,serr_center_control,'-k','CapSize',18,...
    'LineWidth',3);
   
    kplot= bar(1.25,mean_center_Blank,0.1);
    kplot.FaceColor = [0.6350 0.0780 0.1840];
    plot(1.25,Center_Blank,'go','MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','k');
    errorbar(1.25,mean_center_Blank,serr_center_blank,'-k','CapSize',18,...
    'LineWidth',3);
    
    nplot = bar(1.5,mean_center_SSP,0.1);
    nplot.FaceColor = [0.4660 0.6740 0.1880];
    plot(1.5,Center_SSP,'or','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k');
    errorbar(1.5,mean_center_SSP,serr_center_ssp,'-k','CapSize',18,...
    'LineWidth',3);
    xlabel('Experimental Groups');
    ylabel('Pertentage time spent')
    title('Pertentage of time spent in center')
    xlim([0.75 1.75])
    xticks([1 1.25 1.5]);
    xticklabels({'WT Control','Blank-SAP','SSP-SAP'})

    %% 5 min

    figure;
    boxplot(AnimalTable_OFT.total_distance_5min,AnimalTable_OFT.DrugGroup);
    xlabel('Experimental Groups');
    ylabel('Distance Travelled (cm)')
    title('Population comparison of total distance travelled')

figure;
    boxplot(AnimalTable_OFT.center_percentage_5min,AnimalTable_OFT.DrugGroup);
    xlabel('Experimental Groups');
    ylabel('Pertentage time spent')
    title('Pertentage of time spent in center')

figure; % plot mean distance
    [G,GN] = grp2idx(AnimalTable_OFT.DrugGroup);
    Cidx = find(G==2);
    Bidx = find(G==1);
    Sidx = find(G==3);
    Distance_Control = AnimalTable_OFT.total_distance_5min(Cidx,:);
    Distance_Blank = AnimalTable_OFT.total_distance_5min(Bidx,:);
    Distance_SSP = AnimalTable_OFT.total_distance_5min(Sidx,:);

    mean_distance_Control = mean(Distance_Control);
    mean_distance_Blank = mean(Distance_Blank);
    mean_distance_SSP = mean(Distance_SSP);
    serr_distance_control = std(Distance_Control)/sqrt(length(Distance_Control));
    serr_distance_blank = std(Distance_Blank)/sqrt(length(Distance_Blank));
    serr_distance_ssp = std(Distance_SSP)/sqrt(length(Distance_SSP));

    cplot=bar(1,mean_distance_Control,0.1); 
    cplot.FaceColor = [0.9290 0.6940 0.1250];
    hold on
    plot(1,Distance_Control,'ob','MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','k');
    errorbar(1,mean_distance_Control,serr_distance_control,'-k','CapSize',18,...
    'LineWidth',3);
   
    kplot= bar(1.25,mean_distance_Blank,0.1);
    kplot.FaceColor = [0.6350 0.0780 0.1840];
    plot(1.25,Distance_Blank,'go','MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','k');
    errorbar(1.25,mean_distance_Blank,serr_distance_blank,'-k','CapSize',18,...
    'LineWidth',3);
    
    nplot = bar(1.5,mean_distance_SSP,0.1);
    nplot.FaceColor = [0.4660 0.6740 0.1880];
    plot(1.5,Distance_SSP,'or','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k');
    errorbar(1.5,mean_distance_SSP,serr_distance_ssp,'-k','CapSize',18,...
    'LineWidth',3);
    xlabel('Experimental Groups');
    ylabel('Distance Travelled (cm)')
    title('Population comparison of total distance travelled')
    xlim([0.75 1.75])
    xticks([1 1.25 1.5]);
    xticklabels({'WT Control','Blank-SAP','SSP-SAP'})

figure; % plot center percentage
    [G,GN] = grp2idx(AnimalTable_OFT.DrugGroup);
    Cidx = find(G==2);
    Bidx = find(G==1);
    Sidx = find(G==3);
    Center_Control = AnimalTable_OFT.center_percentage_5min(Cidx,:);
    Center_Blank = AnimalTable_OFT.center_percentage_5min(Bidx,:);
    Center_SSP = AnimalTable_OFT.center_percentage_5min(Sidx,:);

    mean_center_Control = mean(Center_Control);
    mean_center_Blank = mean(Center_Blank);
    mean_center_SSP = mean(Center_SSP);
    serr_center_control = std(Center_Control)/sqrt(length(Center_Control));
    serr_center_blank = std(Center_Blank)/sqrt(length(Center_Blank));
    serr_center_ssp = std(Center_SSP)/sqrt(length(Center_SSP));

    cplot=bar(1,mean_center_Control,0.1); 
    cplot.FaceColor = [0.9290 0.6940 0.1250];
    hold on
    plot(1,Center_Control,'ob','MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','k');
    errorbar(1,mean_center_Control,serr_center_control,'-k','CapSize',18,...
    'LineWidth',3);
   
    kplot= bar(1.25,mean_center_Blank,0.1);
    kplot.FaceColor = [0.6350 0.0780 0.1840];
    plot(1.25,Center_Blank,'go','MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','k');
    errorbar(1.25,mean_center_Blank,serr_center_blank,'-k','CapSize',18,...
    'LineWidth',3);
    
    nplot = bar(1.5,mean_center_SSP,0.1);
    nplot.FaceColor = [0.4660 0.6740 0.1880];
    plot(1.5,Center_SSP,'or','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k');
    errorbar(1.5,mean_center_SSP,serr_center_ssp,'-k','CapSize',18,...
    'LineWidth',3);
    xlabel('Experimental Groups');
    ylabel('Pertentage time spent')
    title('Pertentage of time spent in center')
    xlim([0.75 1.75])
    xticks([1 1.25 1.5]);
    xticklabels({'WT Control','Blank-SAP','SSP-SAP'})