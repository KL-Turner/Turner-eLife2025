function plotStatistics_first5min(rootFolder,saveFigs,delim)
% clear all; close all; clc;
load('OFTTable.mat')
% this will plot the stats for first 5 minutes of the experiment.
% you just need to run this code if you are only interested to see the
% final results and you have the other results. You need to be in the
% Session1 folder inside population_compare.
% run the stats
if saveFigs == true
    % statistical diary
    dirpath = [rootFolder delim 'Summary Figures' delim 'Behavior' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    diaryFile = [dirpath 'OpenField_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    [GLME_STATS_distance_5min] = fitglme(AnimalTable_OFT,"total_distance_5min~1+Sex+DrugGroup"); % we are running a GLME stats
    disp('Distance Traveled'); disp(' '); disp(GLME_STATS_distance_5min);
    [GLME_STATS_centertime_5min] = fitglme(AnimalTable_OFT,"center_time_5min~1+Sex+DrugGroup"); % we are running a GLME stats
    disp('Time in Center'); disp(' '); disp(GLME_STATS_centertime_5min);
    diary off
end
WTChoice = 'y';

if WTChoice == 'n'
    %% distance
    %     figure;
    %     boxplot(AnimalTable_OFT.total_distance(13:end),AnimalTable_OFT.DrugGroup(13:end));
    %     xlabel('Experimental Groups');
    %     ylabel('Distance Travelled (cm)')
    %     title('Population comparison of total distance travelled (boxplot)')

    figure; % plot mean distance
    [G,~] = grp2idx(AnimalTable_OFT.DrugGroup);
    Bidx = find(G==1);
    Sidx = find(G==3);
    Distance_Blank = AnimalTable_OFT.total_distance_5min(Bidx,:);
    Distance_SSP = AnimalTable_OFT.total_distance_5min(Sidx,:);

    mean_distance_Blank = mean(Distance_Blank);
    mean_distance_SSP = mean(Distance_SSP);
    serr_distance_blank = std(Distance_Blank)/sqrt(length(Distance_Blank));
    serr_distance_ssp = std(Distance_SSP)/sqrt(length(Distance_SSP));


    kplot= bar(1.25,mean_distance_Blank,0.1);
    kplot.FaceColor = [0.6350 0.0780 0.1840];
    hold on
    plot(1.25,Distance_Blank,'go','MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','k');
    errorbar(1.25,mean_distance_Blank,serr_distance_blank,'-k','CapSize',18,...
        'LineWidth',3);

    nplot = bar(1.5,mean_distance_SSP,0.1);
    nplot.FaceColor = [0.4660 0.6740 0.1880];
    hold on
    plot(1.5,Distance_SSP,'or','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k');
    errorbar(1.5,mean_distance_SSP,serr_distance_ssp,'-k','CapSize',18,...
        'LineWidth',3);

    xlabel('Experimental Groups');
    ylabel('Distance Travelled (cm)')
    title('Population comparison of total distance travelled (mean+-stderror) First 5 minutes')
    xlim([1 1.75])
    xticks([1.25 1.5]);
    xticklabels({'Blank-SAP','SSP-SAP'})

    %%
    figure; % plot mean CenterTime
    [G,~] = grp2idx(AnimalTable_OFT.DrugGroup);
    Bidx = find(G==1);
    Sidx = find(G==3);
    CenterTime_Blank = AnimalTable_OFT.center_time_5min(Bidx,:);
    CenterTime_SSP = AnimalTable_OFT.center_time_5min(Sidx,:);

    mean_CenterTime_Blank = mean(CenterTime_Blank);
    mean_CenterTime_SSP = mean(CenterTime_SSP);
    serr_CenterTime_blank = std(CenterTime_Blank)/sqrt(length(CenterTime_Blank));
    serr_CenterTime_ssp = std(CenterTime_SSP)/sqrt(length(CenterTime_SSP));


    kplot= bar(1.25,mean_CenterTime_Blank,0.1);
    kplot.FaceColor = [0.6350 0.0780 0.1840];
    hold on
    plot(1.25,CenterTime_Blank,'go','MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','k');
    errorbar(1.25,mean_CenterTime_Blank,serr_CenterTime_blank,'-k','CapSize',18,...
        'LineWidth',3);

    nplot = bar(1.5,mean_CenterTime_SSP,0.1);
    nplot.FaceColor = [0.4660 0.6740 0.1880];
    hold on
    plot(1.5,CenterTime_SSP,'or','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k');
    errorbar(1.5,mean_CenterTime_SSP,serr_CenterTime_ssp,'-k','CapSize',18,...
        'LineWidth',3);

    xlabel('Experimental Groups');
    ylabel('Center Time  Percentage')
    title('Population comparison of percentage of time spent in the center (mean+-stderror). First 5 minutes')
    xlim([1 1.75])
    xticks([1.25 1.5]);
    xticklabels({'Blank-SAP','SSP-SAP'})

elseif WTChoice == 'y'
    %     figure;
    %     boxplot(AnimalTable_OFT.total_distance,AnimalTable_OFT.DrugGroup);
    %     xlabel('Experimental Groups');
    %     ylabel('Distance Travelled (cm)')
    %     title('Population comparison of total distance travelled (boxplot)')
    %% distance
    summaryFigure = figure; % plot mean distance
    [G,~] = grp2idx(AnimalTable_OFT.DrugGroup);
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
    cplot.FaceColor = colors('sapphire');
    hold on
    plot(1,Distance_Control,'ob','MarkerSize',8,'MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','k');
    errorbar(1,mean_distance_Control,serr_distance_control,'-k','CapSize',18,...
        'LineWidth',3);

    kplot= bar(1.25,mean_distance_Blank,0.1);
    kplot.FaceColor = colors('north texas green');
    plot(1.25,Distance_Blank,'go','MarkerSize',8,'MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','k');
    errorbar(1.25,mean_distance_Blank,serr_distance_blank,'-k','CapSize',18,...
        'LineWidth',3);

    nplot = bar(1.5,mean_distance_SSP,0.1);
    nplot.FaceColor = colors('electric purple');
    plot(1.5,Distance_SSP,'or','MarkerSize',8,'MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','k');
    errorbar(1.5,mean_distance_SSP,serr_distance_ssp,'-k','CapSize',18,...
        'LineWidth',3);
    xlabel('Experimental Groups');
    ylabel('Distance Travelled (cm)')
    title('Population comparison of total distance travelled (mean+-stderror). First 5 minutes')
    xlim([0.75 1.75])
    xticks([1 1.25 1.5]);
    xticklabels({'Naive','Blank-SAP','SSP-SAP'})
    % save figure(s)
    if saveFigs == true
        savefig(summaryFigure,[dirpath 'OpenField_DistanceTraveled']);
    end
    %%
    summaryFigure = figure; % plot mean CenterTime
    [G,~] = grp2idx(AnimalTable_OFT.DrugGroup);
    Cidx = find(G==2);
    Bidx = find(G==1);
    Sidx = find(G==3);
    CenterTime_Control = AnimalTable_OFT.center_time_5min(Cidx,:);
    CenterTime_Blank = AnimalTable_OFT.center_time_5min(Bidx,:);
    CenterTime_SSP = AnimalTable_OFT.center_time_5min(Sidx,:);

    mean_CenterTime_Control = mean(CenterTime_Control);
    mean_CenterTime_Blank = mean(CenterTime_Blank);
    mean_CenterTime_SSP = mean(CenterTime_SSP);
    serr_CenterTime_control = std(CenterTime_Control)/sqrt(length(CenterTime_Control));
    serr_CenterTime_blank = std(CenterTime_Blank)/sqrt(length(CenterTime_Blank));
    serr_CenterTime_ssp = std(CenterTime_SSP)/sqrt(length(CenterTime_SSP));

    cplot=bar(1,mean_CenterTime_Control,0.1);
    cplot.FaceColor = colors('sapphire');
    hold on
    plot(1,CenterTime_Control,'ob','MarkerSize',8,'MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','k');
    errorbar(1,mean_CenterTime_Control,serr_CenterTime_control,'-k','CapSize',18,...
        'LineWidth',3);

    kplot= bar(1.25,mean_CenterTime_Blank,0.1);
    kplot.FaceColor = colors('north texas green');
    plot(1.25,CenterTime_Blank,'go','MarkerSize',8,'MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','k');
    errorbar(1.25,mean_CenterTime_Blank,serr_CenterTime_blank,'-k','CapSize',18,...
        'LineWidth',3);

    nplot = bar(1.5,mean_CenterTime_SSP,0.1);
    nplot.FaceColor = colors('electric purple');
    plot(1.5,CenterTime_SSP,'or','MarkerSize',8,'MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','k');
    errorbar(1.5,mean_CenterTime_SSP,serr_CenterTime_ssp,'-k','CapSize',18,...
        'LineWidth',3);
    xlabel('Experimental Groups');
    ylabel('Percentage Center Time')
    title('Population comparison of percentage of time spent in the center (mean+-stderror). First 5 minutes.')
    xlim([0.75 1.75])
    xticks([1 1.25 1.5]);
    xticklabels({'Naive','Blank-SAP','SSP-SAP'})
    % save figure(s)
    if saveFigs == true
        savefig(summaryFigure,[dirpath 'OpenField_TimeInCenter']);
    end
end