% function figure_LTA(control_ave, drug1_ave, drug2_ave, control_sem, drug1_sem, drug2_sem)
function [summaryFigure] = figure_LTA(varargin)
if nargin == 2
    control_ave = varargin{1};
    control_sem = varargin{2};
    time = -89/30:1/30:5;
    % mean and std
    summaryFigure = figure;
    box_line(time,control_ave,control_sem,colors('north texas green')); 
    hold on;
    plot(time,control_ave,'color',colors('north texas green'),'LineWidth',2);
    xlim([-3,5]);
    axis square
    set(gca,'box','off')
elseif nargin == 4
    control_ave = varargin{1};
    drug1_ave = varargin{2};
    control_sem = varargin{3};
    drug1_sem = varargin{4};
    time = -89/30:1/30:5;
    % mean and std
    summaryFigure = figure;
    box_line(time,control_ave,control_sem,colors('north texas green')); 
    hold on;
    box_line(time,drug1_ave,drug1_sem,colors('electric purple'));
    plot(time, control_ave,'color',colors('north texas green'),'LineWidth',2);
    plot(time, drug1_ave,'color',colors('electric purple'),'LineWidth',2);
    xlim([-3,5]);
    axis square
    set(gca,'box','off')
elseif nargin == 6
    control_ave = varargin{1};
    drug1_ave = varargin{2};
    drug2_ave = varargin{3};
    control_sem = varargin{4};
    drug1_sem = varargin{5};
    drug2_sem = varargin{6};
    time = -89/30:1/30:5;
    % mean and std
    summaryFigure = figure;
    box_line(time,control_ave,control_sem,colors('north texas green')); 
    hold on;
    box_line(time,drug1_ave,drug1_sem,colors('electric purple'));
    box_line(time,drug2_ave,drug2_sem,'b');
    plot(time, control_ave,'color',colors('north texas green'),'LineWidth',2);
    plot(time, drug1_ave,'color',colors('electric purple'),'LineWidth',2);
    plot(time, drug2_ave,'b');
    xlim([-3,5]);
    axis square
    set(gca,'box','off')
end