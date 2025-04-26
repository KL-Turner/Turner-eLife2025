% function BehaviorMainScript()
% just run this one
        clear all; close all; clc
        plotStatistics() % plot the statistics and figures showing stats
        plotStatistics_first5min(); % plot the stats for the first 5 minutes
        heatmap_plot() % plot the heatmap
        vectormap_plot() % plot the vector map
